/*
 * bam_cb_cache_extract.cpp
 *
 * Build one per-library cell-barcode BAM cache from a coordinate-sorted GEX BAM.
 *
 * Contract:
 *   - One invocation handles one source library.
 *   - Input barcodes come from an allowlist file; no internal library/barcode discovery.
 *   - Reads are emitted only when:
 *       (flag & 0xD04) == 0
 *       AND the requested string tag, normally CB, exists and is in the allowlist.
 *   - No reads are modified. All BAM fields/tags are passed through unchanged.
 *   - Source BAM record order is preserved because this is a single sequential scan.
 *   - Output is a durable BAM + CSI + .meta.json sidecar.
 *   - Existing output is reused only when the sidecar cache key matches current inputs.
 *
 * Build:
 *   g++ -O3 -std=c++17 -Wall -Wextra -pedantic -o bam_cb_cache_extract \
 *       bam_cb_cache_extract.cpp -lhts -lz -lpthread
 *
 * Typical cluster module environment:
 *   module purge
 *   module load htslib/1.20
 */

#include <htslib/sam.h>
#include <htslib/hts.h>

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

static const char* EXTRACTOR_VERSION = "bam_cb_cache_extract_v1.0.0";
static constexpr uint32_t READ_FILTER_MASK = 0xD04;  // unmapped, secondary, duplicate, supplementary
static constexpr int CSI_MIN_SHIFT = 14;

struct Options {
    std::string library_id;
    std::string bam_path;
    std::string allowlist_path;
    std::string tag = "CB";
    int threads = 1;
    std::string output_dir;
    bool force = false;
    bool verbose = false;
};

struct FileSignature {
    uint64_t size = 0;
    int64_t mtime_sec = 0;
    int64_t mtime_nsec = 0;
};

struct AllowlistInfo {
    std::unordered_set<std::string> barcodes;
    std::unordered_map<uint64_t, std::vector<std::string>> lookup_by_hash;
    std::vector<std::string> sorted_barcodes;
    std::string hash_hex;
    size_t raw_nonempty_lines = 0;
    size_t duplicate_lines = 0;
    bool skipped_header = false;
};

struct Counters {
    uint64_t input_records = 0;
    uint64_t flag_pass_records = 0;
    uint64_t missing_tag_records = 0;
    uint64_t malformed_tag_records = 0;
    uint64_t not_in_allowlist_records = 0;
    uint64_t written_records = 0;
};

[[noreturn]] static void die(const std::string& msg) {
    std::cerr << "ERROR: " << msg << "\n";
    std::exit(1);
}

static void log_msg(const std::string& msg) {
    std::cerr << msg << "\n";
}

static std::string shell_quote_for_meta(int argc, char** argv) {
    std::ostringstream oss;
    for (int i = 0; i < argc; ++i) {
        if (i) oss << ' ';
        std::string s(argv[i]);
        oss << '\'';
        for (char c : s) {
            if (c == '\'') oss << "'\\''";
            else oss << c;
        }
        oss << '\'';
    }
    return oss.str();
}

static std::string trim_ascii(const std::string& in) {
    size_t b = 0;
    while (b < in.size() && (in[b] == ' ' || in[b] == '\t' || in[b] == '\r' || in[b] == '\n')) ++b;
    size_t e = in.size();
    while (e > b && (in[e-1] == ' ' || in[e-1] == '\t' || in[e-1] == '\r' || in[e-1] == '\n')) --e;
    return in.substr(b, e - b);
}

static std::string lower_ascii(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return s;
}

static std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    size_t start = 0;
    while (true) {
        size_t pos = line.find('\t', start);
        if (pos == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

static bool looks_like_allowlist_header(const std::string& first_field, const std::string& tag) {
    std::string x = lower_ascii(trim_ascii(first_field));
    std::string t = lower_ascii(tag);
    return x == t || x == "barcode" || x == "barcodes" || x == "cell_barcode" ||
           x == "cellbarcode" || x == "cell" || x == "cb_tag" || x == "cell_id";
}

static std::string json_escape(const std::string& s) {
    std::ostringstream oss;
    for (unsigned char c : s) {
        switch (c) {
            case '\\': oss << "\\\\"; break;
            case '"':  oss << "\\\""; break;
            case '\b': oss << "\\b"; break;
            case '\f': oss << "\\f"; break;
            case '\n': oss << "\\n"; break;
            case '\r': oss << "\\r"; break;
            case '\t': oss << "\\t"; break;
            default:
                if (c < 0x20) {
                    oss << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c)
                        << std::dec << std::setfill(' ');
                } else {
                    oss << static_cast<char>(c);
                }
        }
    }
    return oss.str();
}

static std::string now_iso8601_utc() {
    std::time_t now = std::time(nullptr);
    std::tm tm{};
    gmtime_r(&now, &tm);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return std::string(buf);
}

static std::string hostname_string() {
    char buf[256];
    if (gethostname(buf, sizeof(buf)) == 0) {
        buf[sizeof(buf)-1] = '\0';
        return std::string(buf);
    }
    return "unknown";
}

static bool file_exists_nonempty(const std::string& path) {
    std::error_code ec;
    return fs::is_regular_file(path, ec) && fs::file_size(path, ec) > 0;
}

static void ensure_absolute_path(const std::string& path, const std::string& label) {
    if (path.empty()) die(label + " is empty");
    if (!fs::path(path).is_absolute()) {
        die(label + " must be an absolute path: " + path);
    }
}

static FileSignature stat_file_signature(const std::string& path) {
    struct stat st{};
    if (stat(path.c_str(), &st) != 0) {
        die("cannot stat file " + path + ": " + std::strerror(errno));
    }
    if (!S_ISREG(st.st_mode)) {
        die("path is not a regular file: " + path);
    }
    if (st.st_size <= 0) {
        die("file is empty: " + path);
    }
    FileSignature sig;
    sig.size = static_cast<uint64_t>(st.st_size);
#if defined(__APPLE__)
    sig.mtime_sec = static_cast<int64_t>(st.st_mtimespec.tv_sec);
    sig.mtime_nsec = static_cast<int64_t>(st.st_mtimespec.tv_nsec);
#else
    sig.mtime_sec = static_cast<int64_t>(st.st_mtim.tv_sec);
    sig.mtime_nsec = static_cast<int64_t>(st.st_mtim.tv_nsec);
#endif
    return sig;
}

// SplitMix64/FNV-derived deterministic combiner. This is a content fingerprint, not a cryptographic hash.
static uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static void hash_update_u64(uint64_t& h1, uint64_t& h2, uint64_t x) {
    h1 ^= mix64(x + h2);
    h1 *= 0x100000001b3ULL;
    h2 ^= mix64(x + h1 + 0x517cc1b727220a95ULL);
    h2 *= 0x9e3779b185ebca87ULL;
}

static void hash_update_string(uint64_t& h1, uint64_t& h2, const std::string& s) {
    for (unsigned char c : s) {
        hash_update_u64(h1, h2, static_cast<uint64_t>(c));
    }
    hash_update_u64(h1, h2, 0x0aULL); // delimiter
}

static std::string hex128(uint64_t h1, uint64_t h2) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0') << std::setw(16) << h1 << std::setw(16) << h2;
    return oss.str();
}


static uint64_t barcode_lookup_hash(const std::string& s) {
    uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) {
        h ^= static_cast<uint64_t>(c);
        h *= 1099511628211ULL;
    }
    return h;
}

static uint64_t barcode_lookup_hash_cstr(const char* s) {
    uint64_t h = 1469598103934665603ULL;
    const unsigned char* p = reinterpret_cast<const unsigned char*>(s);
    while (*p) {
        h ^= static_cast<uint64_t>(*p++);
        h *= 1099511628211ULL;
    }
    return h;
}

static bool barcode_allowed(const AllowlistInfo& allow, const char* cb_value) {
    uint64_t h = barcode_lookup_hash_cstr(cb_value);
    auto it = allow.lookup_by_hash.find(h);
    if (it == allow.lookup_by_hash.end()) return false;
    for (const std::string& candidate : it->second) {
        if (candidate == cb_value) return true;
    }
    return false;
}

static std::string fingerprint_strings_sorted(const std::vector<std::string>& items) {
    uint64_t h1 = 0xcbf29ce484222325ULL;
    uint64_t h2 = 0x84222325cbf29ce4ULL;
    hash_update_u64(h1, h2, static_cast<uint64_t>(items.size()));
    for (const auto& s : items) hash_update_string(h1, h2, s);
    return hex128(h1, h2);
}

static std::string fingerprint_cache_key(const Options& opt,
                                         const FileSignature& bam_sig,
                                         const AllowlistInfo& allow) {
    uint64_t h1 = 0x6a09e667f3bcc908ULL;
    uint64_t h2 = 0xbb67ae8584caa73bULL;
    hash_update_string(h1, h2, EXTRACTOR_VERSION);
    hash_update_string(h1, h2, opt.library_id);
    hash_update_string(h1, h2, opt.bam_path);
    hash_update_u64(h1, h2, bam_sig.size);
    hash_update_u64(h1, h2, static_cast<uint64_t>(bam_sig.mtime_sec));
    hash_update_u64(h1, h2, static_cast<uint64_t>(bam_sig.mtime_nsec));
    hash_update_string(h1, h2, opt.allowlist_path);
    hash_update_string(h1, h2, allow.hash_hex);
    hash_update_u64(h1, h2, static_cast<uint64_t>(allow.sorted_barcodes.size()));
    hash_update_string(h1, h2, opt.tag);
    hash_update_u64(h1, h2, READ_FILTER_MASK);
    hash_update_u64(h1, h2, CSI_MIN_SHIFT);
    return hex128(h1, h2);
}

static std::string normalize_library_key(const std::string& library_id) {
    if (library_id.empty()) die("library id is empty");
    std::string key = library_id;
    if (!(key.size() >= 3 && key[0] == 'l' && key[1] == 'i' && key[2] == 'b')) {
        bool all_digits = !key.empty() && std::all_of(key.begin(), key.end(), [](unsigned char c){ return std::isdigit(c); });
        if (all_digits) key = "lib" + key;
    }
    for (char& c : key) {
        unsigned char uc = static_cast<unsigned char>(c);
        if (!(std::isalnum(uc) || c == '_' || c == '-' || c == '.')) c = '_';
    }
    return key;
}

static AllowlistInfo read_allowlist(const std::string& path, const std::string& tag) {
    std::ifstream in(path);
    if (!in) die("failed to open allowlist: " + path);

    AllowlistInfo info;
    std::string line;
    bool saw_first_nonempty = false;
    while (std::getline(in, line)) {
        std::string trimmed = trim_ascii(line);
        if (trimmed.empty()) continue;
        auto fields = split_tab(trimmed);
        std::string bc = trim_ascii(fields.empty() ? trimmed : fields[0]);
        if (bc.empty()) continue;

        if (!saw_first_nonempty) {
            saw_first_nonempty = true;
            if (looks_like_allowlist_header(bc, tag)) {
                info.skipped_header = true;
                continue;
            }
        }

        ++info.raw_nonempty_lines;
        auto inserted = info.barcodes.insert(bc);
        if (!inserted.second) ++info.duplicate_lines;
    }

    if (!in.eof()) die("error while reading allowlist: " + path);
    if (info.barcodes.empty()) die("allowlist is empty after parsing: " + path);

    info.sorted_barcodes.assign(info.barcodes.begin(), info.barcodes.end());
    std::sort(info.sorted_barcodes.begin(), info.sorted_barcodes.end());
    info.hash_hex = fingerprint_strings_sorted(info.sorted_barcodes);
    info.lookup_by_hash.reserve(info.sorted_barcodes.size() * 2 + 1);
    for (const auto& bc : info.sorted_barcodes) {
        info.lookup_by_hash[barcode_lookup_hash(bc)].push_back(bc);
    }
    return info;
}

static void validate_tag_name(const std::string& tag) {
    if (tag.size() != 2) die("tag name must be exactly two characters, got: " + tag);
    for (char c : tag) {
        if (!std::isalnum(static_cast<unsigned char>(c))) {
            die("tag name must be alphanumeric, got: " + tag);
        }
    }
}

static bool get_string_aux_tag(const bam1_t* rec,
                               const std::string& tag,
                               const char*& value,
                               bool& malformed) {
    malformed = false;
    uint8_t* aux = bam_aux_get(rec, tag.c_str());
    if (aux == nullptr) return false;
    char type = static_cast<char>(*aux);
    if (type == 'Z' || type == 'H') {
        value = reinterpret_cast<const char*>(aux + 1);
        return true;
    }
    malformed = true;
    value = nullptr;
    return false;
}

static bool meta_has_cache_key(const std::string& meta_path, const std::string& cache_key) {
    std::ifstream in(meta_path);
    if (!in) return false;
    std::ostringstream ss;
    ss << in.rdbuf();
    std::string txt = ss.str();
    std::string needle = "\"cache_key\": \"" + cache_key + "\"";
    return txt.find(needle) != std::string::npos;
}

static bool validate_bam_and_csi(const std::string& bam_path, const std::string& csi_path) {
    if (!file_exists_nonempty(bam_path) || !file_exists_nonempty(csi_path)) return false;

    samFile* fp = sam_open(bam_path.c_str(), "r");
    if (fp == nullptr) return false;
    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (hdr == nullptr) {
        sam_close(fp);
        return false;
    }
    hts_idx_t* idx = sam_index_load2(fp, bam_path.c_str(), csi_path.c_str());
    bool ok = (idx != nullptr);
    if (idx) hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return ok;
}

static void write_meta_json(const std::string& path,
                            const Options& opt,
                            const std::string& library_key,
                            const std::string& bam_out,
                            const std::string& csi_out,
                            const FileSignature& bam_sig,
                            const AllowlistInfo& allow,
                            const Counters& c,
                            const std::string& cache_key,
                            const std::string& command_line) {
    std::ofstream out(path);
    if (!out) die("failed to open temporary meta for writing: " + path);

    out << "{\n";
    out << "  \"extractor_version\": \"" << json_escape(EXTRACTOR_VERSION) << "\",\n";
    out << "  \"cache_key\": \"" << cache_key << "\",\n";
    out << "  \"created_at_utc\": \"" << now_iso8601_utc() << "\",\n";
    out << "  \"hostname\": \"" << json_escape(hostname_string()) << "\",\n";
    out << "  \"pid\": " << static_cast<long long>(getpid()) << ",\n";
    out << "  \"library_id\": \"" << json_escape(opt.library_id) << "\",\n";
    out << "  \"library_key\": \"" << json_escape(library_key) << "\",\n";
    out << "  \"source_bam\": {\n";
    out << "    \"path\": \"" << json_escape(opt.bam_path) << "\",\n";
    out << "    \"size_bytes\": " << bam_sig.size << ",\n";
    out << "    \"mtime_sec\": " << bam_sig.mtime_sec << ",\n";
    out << "    \"mtime_nsec\": " << bam_sig.mtime_nsec << "\n";
    out << "  },\n";
    out << "  \"allowlist\": {\n";
    out << "    \"path\": \"" << json_escape(opt.allowlist_path) << "\",\n";
    out << "    \"unique_barcodes\": " << allow.sorted_barcodes.size() << ",\n";
    out << "    \"raw_nonempty_data_lines\": " << allow.raw_nonempty_lines << ",\n";
    out << "    \"duplicate_data_lines\": " << allow.duplicate_lines << ",\n";
    out << "    \"skipped_header\": " << (allow.skipped_header ? "true" : "false") << ",\n";
    out << "    \"sorted_hash_type\": \"stable_hash128_sorted_unique_lines_v1\",\n";
    out << "    \"sorted_hash\": \"" << allow.hash_hex << "\"\n";
    out << "  },\n";
    out << "  \"read_filter\": {\n";
    out << "    \"drop_flag_mask_hex\": \"0xD04\",\n";
    out << "    \"drop_unmapped\": true,\n";
    out << "    \"drop_secondary\": true,\n";
    out << "    \"drop_duplicate\": true,\n";
    out << "    \"drop_supplementary\": true,\n";
    out << "    \"mapq_cutoff\": null,\n";
    out << "    \"additional_filters\": []\n";
    out << "  },\n";
    out << "  \"barcode_tag\": \"" << json_escape(opt.tag) << "\",\n";
    out << "  \"threads_requested\": " << opt.threads << ",\n";
    out << "  \"order_preserved\": true,\n";
    out << "  \"reads_modified\": false,\n";
    out << "  \"header_modified\": false,\n";
    out << "  \"output\": {\n";
    out << "    \"bam\": \"" << json_escape(bam_out) << "\",\n";
    out << "    \"csi\": \"" << json_escape(csi_out) << "\",\n";
    out << "    \"read_count\": " << c.written_records << "\n";
    out << "  },\n";
    out << "  \"counters\": {\n";
    out << "    \"input_records\": " << c.input_records << ",\n";
    out << "    \"flag_pass_records\": " << c.flag_pass_records << ",\n";
    out << "    \"missing_tag_records\": " << c.missing_tag_records << ",\n";
    out << "    \"malformed_tag_records\": " << c.malformed_tag_records << ",\n";
    out << "    \"not_in_allowlist_records\": " << c.not_in_allowlist_records << ",\n";
    out << "    \"written_records\": " << c.written_records << "\n";
    out << "  },\n";
    out << "  \"command_line\": \"" << json_escape(command_line) << "\"\n";
    out << "}\n";

    if (!out) die("failed while writing meta json: " + path);
}

static void usage(int code) {
    std::FILE* f = (code == 0) ? stdout : stderr;
    std::fprintf(f,
        "bam_cb_cache_extract %s\n"
        "\n"
        "Build/reuse one per-library CB-allowlist BAM cache.\n"
        "\n"
        "Required:\n"
        "  --library-id, -L   Source library id. Numeric ids are written under lib<N>.\n"
        "  --bam, -b          Absolute path to source coordinate-sorted gex.bam.\n"
        "  --allowlist, -a    Absolute path to barcode allowlist; one CB per line or one-column TSV with header.\n"
        "  --outdir, -o       Absolute cache root. Output goes to <outdir>/cell_bams/<lib>/\n"
        "\n"
        "Optional:\n"
        "  --tag, -t          Barcode tag name. Default: CB.\n"
        "  --threads, -@      htslib BGZF/index threads. Default: 1.\n"
        "  --force, -f        Rebuild even if matching cache exists.\n"
        "  --verbose, -v      Print progress summary.\n"
        "  --version          Print version and exit.\n"
        "  --help, -h         Show this help and exit.\n"
        "\n"
        "Invariant filter: writes reads with (FLAG & 0xD04) == 0 and %s in allowlist.\n",
        EXTRACTOR_VERSION, "TAG");
    std::exit(code);
}

static Options parse_args(int argc, char** argv) {
    Options opt;
    static struct option long_options[] = {
        {"library-id", required_argument, nullptr, 'L'},
        {"bam", required_argument, nullptr, 'b'},
        {"allowlist", required_argument, nullptr, 'a'},
        {"tag", required_argument, nullptr, 't'},
        {"threads", required_argument, nullptr, '@'},
        {"outdir", required_argument, nullptr, 'o'},
        {"force", no_argument, nullptr, 'f'},
        {"verbose", no_argument, nullptr, 'v'},
        {"version", no_argument, nullptr, 1000},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, 0, nullptr, 0}
    };

    if (argc == 1) usage(1);
    int option_index = 0;
    while (true) {
        int c = getopt_long(argc, argv, "L:b:a:t:@:o:fvh", long_options, &option_index);
        if (c == -1) break;
        switch (c) {
            case 'L': opt.library_id = optarg; break;
            case 'b': opt.bam_path = optarg; break;
            case 'a': opt.allowlist_path = optarg; break;
            case 't': opt.tag = optarg; break;
            case '@': opt.threads = std::atoi(optarg); break;
            case 'o': opt.output_dir = optarg; break;
            case 'f': opt.force = true; break;
            case 'v': opt.verbose = true; break;
            case 'h': usage(0); break;
            case 1000:
                std::cout << EXTRACTOR_VERSION << "\n";
                std::exit(0);
            default: usage(1);
        }
    }
    if (optind != argc) usage(1);

    if (opt.library_id.empty()) die("--library-id is required");
    if (opt.bam_path.empty()) die("--bam is required");
    if (opt.allowlist_path.empty()) die("--allowlist is required");
    if (opt.output_dir.empty()) die("--outdir is required");
    if (opt.threads < 1) die("--threads must be >= 1");
    validate_tag_name(opt.tag);
    ensure_absolute_path(opt.bam_path, "--bam");
    ensure_absolute_path(opt.allowlist_path, "--allowlist");
    ensure_absolute_path(opt.output_dir, "--outdir");
    return opt;
}

static Counters extract_bam(const Options& opt,
                            const AllowlistInfo& allow,
                            const std::string& tmp_bam_path) {
    Counters c;
    samFile* in = sam_open(opt.bam_path.c_str(), "r");
    if (in == nullptr) die("failed to open source BAM: " + opt.bam_path);
    if (opt.threads > 1) {
        if (hts_set_threads(in, opt.threads) != 0) {
            sam_close(in);
            die("failed to set input htslib threads");
        }
    }

    bam_hdr_t* hdr = sam_hdr_read(in);
    if (hdr == nullptr) {
        sam_close(in);
        die("failed to read source BAM header: " + opt.bam_path);
    }

    samFile* out = sam_open(tmp_bam_path.c_str(), "wb");
    if (out == nullptr) {
        bam_hdr_destroy(hdr);
        sam_close(in);
        die("failed to open temporary BAM for writing: " + tmp_bam_path);
    }
    if (opt.threads > 1) {
        if (hts_set_threads(out, opt.threads) != 0) {
            sam_close(out);
            bam_hdr_destroy(hdr);
            sam_close(in);
            die("failed to set output htslib threads");
        }
    }

    if (sam_hdr_write(out, hdr) < 0) {
        sam_close(out);
        bam_hdr_destroy(hdr);
        sam_close(in);
        die("failed to write BAM header to: " + tmp_bam_path);
    }

    bam1_t* rec = bam_init1();
    if (rec == nullptr) {
        sam_close(out);
        bam_hdr_destroy(hdr);
        sam_close(in);
        die("failed to allocate BAM record");
    }

    int ret = 0;
    while ((ret = sam_read1(in, hdr, rec)) >= 0) {
        ++c.input_records;
        if ((rec->core.flag & READ_FILTER_MASK) != 0) continue;
        ++c.flag_pass_records;

        const char* cb_value = nullptr;
        bool malformed = false;
        bool has_tag = get_string_aux_tag(rec, opt.tag, cb_value, malformed);
        if (malformed) {
            ++c.malformed_tag_records;
            continue;
        }
        if (!has_tag) {
            ++c.missing_tag_records;
            continue;
        }
        if (!barcode_allowed(allow, cb_value)) {
            ++c.not_in_allowlist_records;
            continue;
        }

        if (sam_write1(out, hdr, rec) < 0) {
            bam_destroy1(rec);
            sam_close(out);
            bam_hdr_destroy(hdr);
            sam_close(in);
            die("failed while writing record to temporary BAM: " + tmp_bam_path);
        }
        ++c.written_records;
    }

    if (ret < -1) {
        bam_destroy1(rec);
        sam_close(out);
        bam_hdr_destroy(hdr);
        sam_close(in);
        die("error while reading source BAM: " + opt.bam_path);
    }

    bam_destroy1(rec);
    bam_hdr_destroy(hdr);

    if (sam_close(out) != 0) {
        sam_close(in);
        die("failed to close temporary output BAM cleanly: " + tmp_bam_path);
    }
    if (sam_close(in) != 0) {
        die("failed to close source BAM cleanly: " + opt.bam_path);
    }

    if (c.malformed_tag_records > 0) {
        die("encountered " + std::to_string(c.malformed_tag_records) +
            " records where tag " + opt.tag + " exists but is not a string tag; refusing cache");
    }
    if (c.written_records == 0) {
        die("output would contain zero reads after -F 0xD04 and " + opt.tag +
            " allowlist filtering; refusing to declare a valid cache");
    }
    return c;
}

static void build_csi_index(const std::string& bam_path, const std::string& csi_path, int threads) {
    if (sam_index_build3(bam_path.c_str(), csi_path.c_str(), CSI_MIN_SHIFT, threads) != 0) {
        die("failed to build CSI index for temporary BAM: " + bam_path);
    }
    if (!file_exists_nonempty(csi_path)) {
        die("CSI index was not created or is empty: " + csi_path);
    }
}

static void atomic_rename(const std::string& src, const std::string& dst) {
    std::error_code ec;
    fs::rename(src, dst, ec);
    if (ec) {
        // Some filesystems/libraries refuse replacing existing targets through std::filesystem.
        std::error_code rm_ec;
        fs::remove(dst, rm_ec);
        ec.clear();
        fs::rename(src, dst, ec);
        if (ec) die("failed to rename " + src + " -> " + dst + ": " + ec.message());
    }
}

int main(int argc, char** argv) {
    const std::string command_line = shell_quote_for_meta(argc, argv);
    Options opt = parse_args(argc, argv);

    FileSignature bam_sig = stat_file_signature(opt.bam_path);
    (void)stat_file_signature(opt.allowlist_path);
    AllowlistInfo allow = read_allowlist(opt.allowlist_path, opt.tag);
    std::string cache_key = fingerprint_cache_key(opt, bam_sig, allow);

    std::string library_key = normalize_library_key(opt.library_id);
    fs::path lib_dir = fs::path(opt.output_dir) / "cell_bams" / library_key;
    fs::path final_bam = lib_dir / (library_key + ".cells.bam");
    fs::path final_csi = lib_dir / (library_key + ".cells.bam.csi");
    fs::path final_meta = lib_dir / (library_key + ".cells.bam.meta.json");

    std::error_code ec;
    fs::create_directories(lib_dir, ec);
    if (ec) die("failed to create output directory " + lib_dir.string() + ": " + ec.message());

    if (!opt.force &&
        file_exists_nonempty(final_bam.string()) &&
        file_exists_nonempty(final_csi.string()) &&
        file_exists_nonempty(final_meta.string()) &&
        meta_has_cache_key(final_meta.string(), cache_key) &&
        validate_bam_and_csi(final_bam.string(), final_csi.string())) {
        std::cout << "REUSED\t" << opt.library_id << "\t" << final_bam.string()
                  << "\t" << final_csi.string() << "\t" << final_meta.string() << "\n";
        return 0;
    }

    const auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    fs::path tmp_dir = lib_dir / (".tmp." + library_key + "." + std::to_string(getpid()) + "." + std::to_string(now_ns));
    fs::create_directories(tmp_dir, ec);
    if (ec) die("failed to create temporary directory " + tmp_dir.string() + ": " + ec.message());

    fs::path tmp_bam = tmp_dir / (library_key + ".cells.bam.tmp");
    fs::path tmp_csi = tmp_dir / (library_key + ".cells.bam.tmp.csi");
    fs::path tmp_meta = tmp_dir / (library_key + ".cells.bam.meta.json.tmp");

    if (opt.verbose) {
        log_msg("Building cache for library_id=" + opt.library_id);
        log_msg("  source BAM: " + opt.bam_path);
        log_msg("  allowlist unique barcodes: " + std::to_string(allow.sorted_barcodes.size()));
        log_msg("  output BAM: " + final_bam.string());
    }

    Counters counters = extract_bam(opt, allow, tmp_bam.string());
    if (!file_exists_nonempty(tmp_bam.string())) die("temporary BAM is empty: " + tmp_bam.string());

    build_csi_index(tmp_bam.string(), tmp_csi.string(), opt.threads);
    if (!validate_bam_and_csi(tmp_bam.string(), tmp_csi.string())) {
        die("temporary BAM/CSI failed validation: " + tmp_bam.string());
    }

    write_meta_json(tmp_meta.string(), opt, library_key, final_bam.string(), final_csi.string(),
                    bam_sig, allow, counters, cache_key, command_line);
    if (!file_exists_nonempty(tmp_meta.string())) die("temporary meta is empty: " + tmp_meta.string());

    atomic_rename(tmp_bam.string(), final_bam.string());
    atomic_rename(tmp_csi.string(), final_csi.string());
    atomic_rename(tmp_meta.string(), final_meta.string());

    fs::remove_all(tmp_dir, ec); // best-effort cleanup

    if (!validate_bam_and_csi(final_bam.string(), final_csi.string())) {
        die("final BAM/CSI failed validation after rename: " + final_bam.string());
    }

    std::cout << "BUILT\t" << opt.library_id << "\t" << final_bam.string()
              << "\t" << final_csi.string() << "\t" << final_meta.string()
              << "\treads=" << counters.written_records << "\n";
    if (opt.verbose) {
        log_msg("  input_records=" + std::to_string(counters.input_records));
        log_msg("  flag_pass_records=" + std::to_string(counters.flag_pass_records));
        log_msg("  missing_tag_records=" + std::to_string(counters.missing_tag_records));
        log_msg("  not_in_allowlist_records=" + std::to_string(counters.not_in_allowlist_records));
        log_msg("  written_records=" + std::to_string(counters.written_records));
    }
    return 0;
}
