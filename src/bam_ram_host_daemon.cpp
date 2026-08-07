/*
 * bam_ram_host_daemon.cpp   (version V1_R4; see revision history)
 *
 *
 * BAM RAM-host daemon (C++ replacement for bam_ram_host_daemon.py).
 *
 * Stages a set of BAM files into a node-local RAM-backed filesystem (/dev/shm,
 * tmpfs) as real files, so every job on that node reads them from one explicitly
 * owned resident copy at a known path, opened normally by htslib/samtools. This
 * is explicit residency under our control, not the page cache: files stay until
 * we tear them down.
 *
 * Unlike vcf_loader_daemon (which parses VCFs into a POSIX shared-memory data
 * structure that demux_parallel attaches to), BAMs are consumed as files by
 * path, so this daemon copies the files into tmpfs and publishes a manifest of
 * orig -> hosted paths. A client (the bounded sampler in
 * synth_composition_injection.py) resolves a wanted BAM by reading
 * <shm-dir>/hosted_manifest.tsv: hosted path if present, else the original path
 * (or, under SYNTH_BAM_HOST_REQUIRE, a hard error).
 *
 * Each BAM is copied as-is ("compressed"); its .csi sidecar is colocated when
 * present so the hosted copy is a complete indexed BAM (reads are full scans, so
 * the index is not required, but colocating it is cheap and drop-in safe). The
 * daemon hosts files in manifest order while they fit under the cap and the
 * tmpfs has room (minus a reserve left for other /dev/shm users such as the
 * per-task RAMDISKs), and leaves anything that does not fit unhosted; clients
 * fall back to the original path unless SYNTH_BAM_HOST_REQUIRE forbids it.
 *
 * Modes:
 *   (default) host : clean the shm-dir, stage files, write hosted_manifest.tsv,
 *                    then hold (foreground or daemonized) until SIGTERM/SIGINT,
 *                    on which it removes the shm-dir (clean teardown).
 *   --oneshot      : stage and exit, leaving files resident (for srun staging
 *                    and week-long holds); teardown is then --destroy.
 *   --destroy      : remove a previously hosted shm-dir and exit.
 *   --status       : print the hosted manifest and exit.
 *
 * Manifest format (tab-separated, one header line):
 *   orig_path<TAB>hosted_path<TAB>form<TAB>bytes
 * orig_path is the realpath of the source BAM so it matches a client that
 * canonicalizes the path it looks up.
 *
 * Build (added to the Makefile parallel_tools target):
 *   g++ -std=c++11 -O3 -o bam_ram_host_daemon src/bam_ram_host_daemon.cpp -lpthread
 *
 * Revision history:
 *   V1_R1: Initial C++ daemon. tmpfs file staging with statvfs-bounded cap and a
 *          reserve, parallel copy across --threads workers, atomic manifest
 *          write, foreground/daemonized hold with signal teardown, --oneshot
 *          persist, and --destroy/--status modes. Replaces the Python daemon.
 *   V1_R2: Diagnostics only; hosting behavior byte-identical. Copy failures now
 *          report the failing call's errno (number + strerror) and bytes written
 *          so far. The threaded stager additionally logs the shm-dir state at the
 *          moment of failure (shmdir_exists / shmdir_writable) to separate
 *          "dir vanished/perms" from "the write itself failed". The manifest open
 *          failure appends errno plus the ofstream state bits. Staging start logs
 *          the resolved shm-dir, daemon uid/gid, dir_preexists, and the mkdir
 *          rc/errno. errno is captured into a local immediately after each failing
 *          call (before any logging call can clobber it) and strerror is read via
 *          a thread-safe strerror_r wrapper. Internal version string bumped to
 *          v1.1.0 so the rebuilt binary is identifiable in logs.
 *   V1_R3: The explicit --cap-gb is now authoritative for the staging budget instead
 *          of only ever lowering the 10% auto-reserve. budget = min(cap_bytes,
 *          avail - floor_reserve) when a cap is given, so a cap sized to the known
 *          cache footprint raises the budget above the conservative auto-reserve and
 *          binds, while a cap can never exceed what tmpfs physically holds. Fixes the
 *          holder over-cap skip (B27): at an 850 GiB tmpfs the auto budget was 765 GiB,
 *          below the 766.47 GiB cache, so one cache (lib40) was skipped and the
 *          all-or-nothing residency gate failed; the orchestrator's --cap-gb 783 was
 *          inert because it could only lower 765, never raise it. Now 783 binds and all
 *          40 caches stage. Budget/sizing only: staged bytes are byte-identical copies;
 *          which files are RAM-resident vs read from BeeGFS is the only thing that
 *          changes, and the manifest and teardown paths are untouched.
 *   V1_R4: Version string only; hosting/budget behavior byte-identical to V1_R3.
 *          The V1_R3 cap-authoritative change did not bump DAEMON_VERSION, so a
 *          freshly built V1_R3 binary and a stale V1_R2-era binary both logged
 *          "...v1.1.0" at startup, leaving the only runtime discriminator of the
 *          cap fix the budget value in the staging log. DAEMON_VERSION is bumped to
 *          v1.2.0 so a binary carrying the cap-authoritative budget logic is
 *          identifiable from its first startup line; no code path, sizing, manifest,
 *          or teardown behavior changes.
 */

#include <getopt.h>
#include <unistd.h>
#include <dirent.h>
#include <signal.h>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <sys/stat.h>
#include <sys/statvfs.h>
#include <sys/types.h>

using namespace std;

static const char* DAEMON_VERSION = "bam_ram_host_daemon_v1.2.0";
static const char* HOSTED_MANIFEST_NAME = "hosted_manifest.tsv";

static volatile sig_atomic_t keep_running = 1;

static void signal_handler(int) {
    keep_running = 0;
}

// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

static void log_line(const string& msg) {
    char ts[32];
    time_t now = time(nullptr);
    struct tm tmv;
    localtime_r(&now, &tmv);
    strftime(ts, sizeof(ts), "%H:%M:%S", &tmv);
    fprintf(stderr, "[bam_ram_host %s] %s\n", ts, msg.c_str());
}

[[noreturn]] static void die(const string& msg) {
    log_line("ERROR: " + msg);
    exit(1);
}

// Thread-safe strerror. The stager logs from multiple worker threads, so plain
// strerror (shared static buffer) would race; strerror_r is the safe form. The
// two xsi_select overloads pick the right return-type handling at compile time
// for either the GNU (char*) or XSI/POSIX (int) strerror_r signature.
static inline const char* xsi_select(int ret, const char* buf) {
    return (ret == 0) ? buf : "unknown error";
}
static inline const char* xsi_select(char* ret, const char* /*buf*/) {
    return ret;
}
static string errno_str(int e) {
    char buf[256];
    buf[0] = '\0';
    return string(xsi_select(strerror_r(e, buf, sizeof(buf)), buf));
}

static string basename_of(const string& path) {
    size_t pos = path.find_last_of('/');
    return (pos == string::npos) ? path : path.substr(pos + 1);
}

static bool is_regular_nonempty(const string& path, uint64_t* size_out) {
    struct stat st;
    if (stat(path.c_str(), &st) != 0) return false;
    if (!S_ISREG(st.st_mode)) return false;
    if (st.st_size <= 0) return false;
    if (size_out) *size_out = static_cast<uint64_t>(st.st_size);
    return true;
}

static string realpath_or_self(const string& path) {
    char buf[PATH_MAX];
    if (realpath(path.c_str(), buf) != nullptr) return string(buf);
    return path;  // fall back to the given path if it cannot be resolved
}

static int mkdir_p(const string& dir) {
    string acc;
    size_t i = 0;
    if (!dir.empty() && dir[0] == '/') { acc = "/"; i = 1; }
    while (i <= dir.size()) {
        if (i == dir.size() || dir[i] == '/') {
            if (!acc.empty() && acc != "/") {
                if (mkdir(acc.c_str(), 0755) != 0 && errno != EEXIST) return -1;
            }
            if (i < dir.size()) acc += '/';
        } else {
            acc += dir[i];
        }
        ++i;
    }
    return 0;
}

// Recursive remove of a directory tree (the shm-dir is flat, but recurse to be safe).
static void remove_tree(const string& path) {
    struct stat st;
    if (lstat(path.c_str(), &st) != 0) return;
    if (S_ISDIR(st.st_mode)) {
        DIR* d = opendir(path.c_str());
        if (d != nullptr) {
            struct dirent* ent;
            while ((ent = readdir(d)) != nullptr) {
                string name = ent->d_name;
                if (name == "." || name == "..") continue;
                remove_tree(path + "/" + name);
            }
            closedir(d);
        }
        rmdir(path.c_str());
    } else {
        unlink(path.c_str());
    }
}

// Buffered copy of src -> dst. Returns bytes copied, or -1 on error.
// Copy mechanism: C stdio FILE* streams (fopen/fread/fwrite/fflush/fclose);
// glibc sets errno on these, so errno is the correct failure channel here.
// On error, *err_out receives the errno captured at the FIRST failing call
// (captured into a local immediately, before any later call can clobber it) and
// *written_out receives the bytes successfully written so far. Copy logic, order,
// and control flow are unchanged from V1_R1; only failure diagnostics are added.
static int64_t copy_file(const string& src, const string& dst,
                         int* err_out, int64_t* written_out) {
    if (err_out) *err_out = 0;
    if (written_out) *written_out = 0;
    FILE* in = fopen(src.c_str(), "rb");
    if (in == nullptr) {
        int e = errno;                       // capture before anything else
        if (err_out) *err_out = e;
        return -1;
    }
    FILE* out = fopen(dst.c_str(), "wb");
    if (out == nullptr) {
        int e = errno;                       // capture before fclose clobbers it
        if (err_out) *err_out = e;
        fclose(in);
        return -1;
    }
    const size_t BUF = 16 * 1024 * 1024;  // 16 MiB
    vector<char> buf(BUF);
    int64_t total = 0;
    bool ok = true;
    int first_errno = 0;                     // errno of the first failing op
    while (true) {
        size_t n = fread(buf.data(), 1, BUF, in);
        if (n == 0) {
            if (ferror(in)) { if (first_errno == 0) first_errno = errno; ok = false; }
            break;
        }
        if (fwrite(buf.data(), 1, n, out) != n) { if (first_errno == 0) first_errno = errno; ok = false; break; }
        total += static_cast<int64_t>(n);
    }
    if (fflush(out) != 0) { if (first_errno == 0) first_errno = errno; ok = false; }
    fclose(in);
    if (fclose(out) != 0) { if (first_errno == 0) first_errno = errno; ok = false; }
    if (written_out) *written_out = total;
    if (!ok) {
        if (err_out) *err_out = first_errno;  // store before unlink() can change errno
        unlink(dst.c_str());
        return -1;
    }
    return total;
}

static uint64_t fs_available_bytes(const string& dir) {
    struct statvfs vfs;
    if (statvfs(dir.c_str(), &vfs) != 0) return 0;
    return static_cast<uint64_t>(vfs.f_bavail) * static_cast<uint64_t>(vfs.f_frsize);
}

// ---------------------------------------------------------------------------
// manifest
// ---------------------------------------------------------------------------

struct HostedRow {
    string orig;    // realpath of the source BAM
    string hosted;  // tmpfs path
    string form;    // "compressed"
    uint64_t bytes; // hosted BAM bytes
};

static void write_manifest(const string& shm_dir, const vector<HostedRow>& rows) {
    string final_path = shm_dir + "/" + HOSTED_MANIFEST_NAME;
    string tmp_path = final_path + ".tmp";
    errno = 0;
    ofstream out(tmp_path.c_str(), ios::out | ios::trunc);
    int open_errno = errno;   // capture immediately after the stream op
    if (!out.is_open()) {
        die("cannot open manifest for writing: " + tmp_path +
            ": errno=" + to_string(open_errno) + " (" + errno_str(open_errno) + ")" +
            " [stream fail=" + to_string(out.fail() ? 1 : 0) +
            " bad=" + to_string(out.bad() ? 1 : 0) +
            " eof=" + to_string(out.eof() ? 1 : 0) + "]");
    }
    out << "orig_path\thosted_path\tform\tbytes\n";
    for (const HostedRow& r : rows) {
        out << r.orig << '\t' << r.hosted << '\t' << r.form << '\t' << r.bytes << '\n';
    }
    out.flush();
    if (!out.good()) { out.close(); unlink(tmp_path.c_str()); die("failed writing manifest: " + tmp_path); }
    out.close();
    if (rename(tmp_path.c_str(), final_path.c_str()) != 0) {
        unlink(tmp_path.c_str());
        die("failed to install manifest " + final_path + ": " + strerror(errno));
    }
}

static void print_manifest(const string& shm_dir) {
    string path = shm_dir + "/" + HOSTED_MANIFEST_NAME;
    ifstream in(path.c_str());
    if (!in.is_open()) { log_line("no hosted manifest at " + path); return; }
    string line;
    while (getline(in, line)) cout << line << "\n";
}

// ---------------------------------------------------------------------------
// staging
// ---------------------------------------------------------------------------

static vector<string> read_manifest_paths(const string& manifest) {
    vector<string> out;
    ifstream in(manifest.c_str());
    if (!in.is_open()) die("cannot open --manifest: " + manifest);
    string line;
    while (getline(in, line)) {
        size_t b = line.find_first_not_of(" \t\r\n");
        if (b == string::npos) continue;
        size_t e = line.find_last_not_of(" \t\r\n");
        string p = line.substr(b, e - b + 1);
        if (p.empty() || p[0] == '#') continue;
        out.push_back(p);
    }
    return out;
}

struct StageState {
    vector<string> srcs;          // source BAM paths (manifest order)
    string shm_dir;
    uint64_t remaining_budget;    // bytes still allowed to host
    size_t next_index;            // next file to claim
    vector<HostedRow> rows;       // successful hosts
    size_t skipped;               // over-budget skips
    bool had_error;               // a copy failed
    mutex mtx;
};

static void stage_worker(StageState* s) {
    while (true) {
        size_t idx;
        string src;
        uint64_t need = 0;
        string csi;
        uint64_t csi_bytes = 0;
        bool have_csi = false;
        // claim a file and reserve its budget under the lock
        {
            lock_guard<mutex> lk(s->mtx);
            if (s->next_index >= s->srcs.size()) return;
            idx = s->next_index++;
            src = s->srcs[idx];
            uint64_t bam_bytes = 0;
            if (!is_regular_nonempty(src, &bam_bytes)) {
                s->had_error = true;
                log_line("source BAM missing or empty: " + src);
                continue;
            }
            csi = src + ".csi";
            if (is_regular_nonempty(csi, &csi_bytes)) have_csi = true;
            need = bam_bytes + (have_csi ? csi_bytes : 0);
            if (need > s->remaining_budget) {
                s->skipped++;
                log_line("over cap, NOT hosting (clients fall back): " + src);
                continue;
            }
            s->remaining_budget -= need;
        }
        // copy outside the lock
        string base = basename_of(src);
        string dst = s->shm_dir + "/" + base;
        int copy_errno = 0;
        int64_t copy_wrote = 0;
        int64_t copied = copy_file(src, dst, &copy_errno, &copy_wrote);
        if (copied < 0) {
            // shm-dir state at the moment of failure (separates "dir vanished/perms"
            // from "the write itself failed"); copy_errno is already saved, so the
            // stat/access errno does not matter here.
            struct stat dst_dir_st;
            int shmdir_exists = (stat(s->shm_dir.c_str(), &dst_dir_st) == 0) ? 1 : 0;
            int shmdir_writable = (access(s->shm_dir.c_str(), W_OK) == 0) ? 1 : 0;
            lock_guard<mutex> lk(s->mtx);
            s->had_error = true;
            log_line("failed to copy " + src + " -> " + dst +
                     ": errno=" + to_string(copy_errno) +
                     " (" + errno_str(copy_errno) + ")" +
                     " [wrote " + to_string(copy_wrote) + " bytes]" +
                     " shmdir_exists=" + to_string(shmdir_exists) +
                     " shmdir_writable=" + to_string(shmdir_writable));
            continue;
        }
        if (have_csi) {
            string dst_csi = dst + ".csi";
            int csi_errno = 0;
            int64_t csi_wrote = 0;
            if (copy_file(csi, dst_csi, &csi_errno, &csi_wrote) < 0) {
                struct stat dst_dir_st;
                int shmdir_exists = (stat(s->shm_dir.c_str(), &dst_dir_st) == 0) ? 1 : 0;
                int shmdir_writable = (access(s->shm_dir.c_str(), W_OK) == 0) ? 1 : 0;
                lock_guard<mutex> lk(s->mtx);
                s->had_error = true;
                log_line("failed to copy index " + csi + " -> " + dst_csi +
                         ": errno=" + to_string(csi_errno) +
                         " (" + errno_str(csi_errno) + ")" +
                         " [wrote " + to_string(csi_wrote) + " bytes]" +
                         " shmdir_exists=" + to_string(shmdir_exists) +
                         " shmdir_writable=" + to_string(shmdir_writable));
                unlink(dst.c_str());
                continue;
            }
        }
        HostedRow row;
        row.orig = realpath_or_self(src);
        row.hosted = dst;
        row.form = "compressed";
        row.bytes = static_cast<uint64_t>(copied);
        {
            lock_guard<mutex> lk(s->mtx);
            s->rows.push_back(row);
            log_line("hosted " + base + " (" +
                     to_string(copied / (1024 * 1024)) + " MiB) at " + dst);
        }
    }
}

static int cmd_host(const string& shm_dir, const string& manifest, int threads,
                    const string& cap_gb, bool oneshot, bool foreground);

static int cmd_host(const string& shm_dir, const string& manifest, int threads,
                    const string& cap_gb, bool oneshot, bool foreground) {
    vector<string> srcs = read_manifest_paths(manifest);
    if (srcs.empty()) die("no BAM paths in --manifest: " + manifest);

    // Does the shm-dir already exist before we (re)create it? (read-only probe)
    struct stat pre_st;
    int dir_preexists = (stat(shm_dir.c_str(), &pre_st) == 0) ? 1 : 0;

    // Clean any stale residency, then create the dir fresh.
    remove_tree(shm_dir);
    errno = 0;
    int mkdir_rc = mkdir_p(shm_dir);
    int mkdir_errno = errno;   // capture immediately after mkdir_p
    log_line("staging into " + shm_dir +
             " (uid=" + to_string(static_cast<unsigned>(getuid())) +
             " gid=" + to_string(static_cast<unsigned>(getgid())) + ")" +
             " dir_preexists=" + to_string(dir_preexists) +
             " mkdir_rc=" + to_string(mkdir_rc) +
             " errno=" + to_string(mkdir_rc == 0 ? 0 : mkdir_errno));
    if (mkdir_rc != 0) die("failed to create shm-dir " + shm_dir + ": " + errno_str(mkdir_errno));

    // Budget: leave a reserve for other /dev/shm users (per-task RAMDISKs, etc.).
    // Reserve scales with what is actually available so it can never exceed it.
    uint64_t avail = fs_available_bytes(shm_dir);
    uint64_t reserve = avail / 10;                       // ~10% of available
    uint64_t floor_reserve = static_cast<uint64_t>(8) * 1024 * 1024 * 1024;  // >= 8 GiB
    if (reserve < floor_reserve) reserve = floor_reserve;
    uint64_t budget = (avail > reserve) ? (avail - reserve) : 0;
    if (cap_gb != "auto") {
        char* end = nullptr;
        double g = strtod(cap_gb.c_str(), &end);
        if (end == cap_gb.c_str() || g <= 0) die("--cap-gb must be 'auto' or a positive number: " + cap_gb);
        uint64_t cap_bytes = static_cast<uint64_t>(g * 1024.0 * 1024.0 * 1024.0);
        // An explicit --cap-gb is authoritative: it sets the staging budget directly,
        // raising it above the conservative 10% auto-reserve when the caller has sized
        // the cap to the known cache footprint, or lowering it below the auto-reserve
        // when asked. It is bounded only by what the device can physically hold
        // (avail - floor_reserve), so a too-large cap can never exceed the tmpfs.
        // Previously the cap could only LOWER the auto budget, so a cap set just above
        // the cache (e.g. 783 GiB) was inert against the 0.9*avail auto-reserve (765 GiB
        // at an 850 GiB tmpfs) and one file was skipped over-cap even though it fit.
        uint64_t hard_ceiling = (avail > floor_reserve) ? (avail - floor_reserve) : 0;
        budget = (cap_bytes < hard_ceiling) ? cap_bytes : hard_ceiling;
    }
    log_line("staging " + to_string(srcs.size()) + " BAM(s) into " + shm_dir +
             "; tmpfs avail " + to_string(avail / (1024ULL * 1024 * 1024)) + " GiB, budget " +
             to_string(budget / (1024ULL * 1024 * 1024)) + " GiB");

    StageState s;
    s.srcs = srcs;
    s.shm_dir = shm_dir;
    s.remaining_budget = budget;
    s.next_index = 0;
    s.skipped = 0;
    s.had_error = false;

    int nthreads = threads;
    if (nthreads < 1) nthreads = 1;
    if (static_cast<size_t>(nthreads) > srcs.size()) nthreads = static_cast<int>(srcs.size());
    vector<thread> pool;
    for (int i = 0; i < nthreads; ++i) pool.emplace_back(stage_worker, &s);
    for (thread& t : pool) t.join();

    write_manifest(shm_dir, s.rows);
    log_line("hosted " + to_string(s.rows.size()) + "/" + to_string(srcs.size()) +
             " file(s) in " + shm_dir + " (skipped over-cap: " + to_string(s.skipped) + ")");

    if (s.had_error) die("one or more files failed to stage; manifest reflects only successful hosts");

    if (oneshot) {
        log_line("oneshot: files left resident; exiting");
        return 0;
    }

    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    if (!foreground) {
        pid_t pid = fork();
        if (pid < 0) { remove_tree(shm_dir); die(string("fork failed: ") + strerror(errno)); }
        if (pid > 0) {
            fprintf(stderr, "\nDaemon started with PID %d\n", pid);
            fprintf(stderr, "Hosted %zu file(s) in %s\n", s.rows.size(), shm_dir.c_str());
            fprintf(stderr, "To destroy: bam_ram_host_daemon --destroy --shm-dir %s\n", shm_dir.c_str());
            exit(0);
        }
        setsid();
        close(STDIN_FILENO);
        close(STDOUT_FILENO);
        close(STDERR_FILENO);
    } else {
        log_line("holding resident (foreground); SIGTERM/SIGINT to release");
    }

    while (keep_running) sleep(1);

    if (foreground) log_line("signal received; tearing down " + shm_dir);
    remove_tree(shm_dir);
    return 0;
}

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

static void help(int code) {
    fprintf(stderr, "bam_ram_host_daemon [OPTIONS]\n");
    fprintf(stderr, "Stage BAM files into node-local tmpfs (/dev/shm) for fast shared reads.\n\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED (host mode) =====\n");
    fprintf(stderr, "    --manifest -m   File with one absolute BAM path per line\n");
    fprintf(stderr, "    --shm-dir -s    Absolute node-local tmpfs dir (e.g. /dev/shm/bam_host_<name>)\n");
    fprintf(stderr, "                    (or pass --name to derive /dev/shm/bam_host_<name>)\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --name -n       Host name; shm-dir becomes /dev/shm/bam_host_<name>\n");
    fprintf(stderr, "    --threads -t    Parallel copy workers (default: 8)\n");
    fprintf(stderr, "    --cap-gb -c     Per-node cap in GB, or 'auto' (default: auto = tmpfs avail minus reserve)\n");
    fprintf(stderr, "    --oneshot -1    Stage and exit, leaving files resident\n");
    fprintf(stderr, "    --foreground -f Hold in the foreground instead of daemonizing\n");
    fprintf(stderr, "    --destroy -d    Remove the shm-dir and exit\n");
    fprintf(stderr, "    --status        Print the hosted manifest and exit\n");
    fprintf(stderr, "    --help -h       Display this message and exit\n");
    exit(code);
}

int main(int argc, char* argv[]) {
    static struct option long_options[] = {
        {"manifest",   required_argument, 0, 'm'},
        {"shm-dir",    required_argument, 0, 's'},
        {"name",       required_argument, 0, 'n'},
        {"threads",    required_argument, 0, 't'},
        {"cap-gb",     required_argument, 0, 'c'},
        {"oneshot",    no_argument,       0, '1'},
        {"foreground", no_argument,       0, 'f'},
        {"destroy",    no_argument,       0, 'd'},
        {"status",     no_argument,       0, 1001},
        {"help",       no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

    string manifest;
    string shm_dir;
    string name;
    int threads = 8;
    string cap_gb = "auto";
    bool oneshot = false;
    bool foreground = false;
    bool destroy_only = false;
    bool status_only = false;

    if (argc == 1) help(0);

    int option_index = 0;
    int ch;
    while ((ch = getopt_long(argc, argv, "m:s:n:t:c:1fdh", long_options, &option_index)) != -1) {
        switch (ch) {
            case 'm': manifest = optarg; break;
            case 's': shm_dir = optarg; break;
            case 'n': name = optarg; break;
            case 't': threads = atoi(optarg); break;
            case 'c': cap_gb = optarg; break;
            case '1': oneshot = true; break;
            case 'f': foreground = true; break;
            case 'd': destroy_only = true; break;
            case 1001: status_only = true; break;
            case 'h': help(0); break;
            default: help(1);
        }
    }

    if (shm_dir.empty() && !name.empty()) shm_dir = "/dev/shm/bam_host_" + name;
    if (shm_dir.empty()) die("--shm-dir (or --name) is required");
    if (shm_dir[0] != '/') die("--shm-dir must be an absolute path: " + shm_dir);

    if (destroy_only) {
        log_line("destroying " + shm_dir);
        remove_tree(shm_dir);
        return 0;
    }
    if (status_only) {
        print_manifest(shm_dir);
        return 0;
    }

    if (manifest.empty()) die("--manifest is required for host mode");
    if (oneshot && foreground) die("--oneshot and --foreground are mutually exclusive");

    log_line(string(DAEMON_VERSION) + " host -> " + shm_dir);
    return cmd_host(shm_dir, manifest, threads, cap_gb, oneshot, foreground);
}
