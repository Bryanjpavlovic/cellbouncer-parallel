#include <getopt.h>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <vector>
#include <signal.h>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <htswrapper/bam.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_parallel_hts.h"

using namespace std;

static volatile sig_atomic_t keep_running = 1;

// Track all created segments for cleanup
static vector<string> active_segments;

void signal_handler(int sig) {
    keep_running = 0;
}

// Intersect BAM chromosomes with a VCF file, optionally filtered by a user chromosome list.
// Returns the set of chromosomes present in both BAM and VCF (and in user list if provided).
set<string> get_shared_chroms(const string& vcf_file,
                              const set<string>& chroms_bam,
                              const string& chroms_file) {
    set<string> chroms_vcf;
    string vcf_mut = vcf_file;
    get_vcf_chroms(vcf_mut, chroms_vcf);
    
    set<string> shared;
    for (auto& c : chroms_bam) {
        if (chroms_vcf.count(c)) {
            shared.insert(c);
        }
    }
    
    if (!chroms_file.empty()) {
        set<string> user_chroms;
        ifstream cf(chroms_file);
        string line;
        while (getline(cf, line)) {
            if (!line.empty()) {
                user_chroms.insert(line);
            }
        }
        set<string> filtered;
        for (auto& c : shared) {
            if (user_chroms.count(c)) {
                filtered.insert(c);
            }
        }
        return filtered;
    }
    return shared;
}

// Load a VCF into shared memory with the given segment name.
// Returns true if segment was created successfully.
// On success, adds shm_name to active_segments for cleanup tracking.
bool load_vcf_segment(const string& vcf_file,
                      const string& shm_name,
                      const string& label,
                      const set<string>& chroms_bam,
                      const string& chroms_file,
                      map<string, int>& seq2tid,
                      int min_vq) {
    
    set<string> chroms = get_shared_chroms(vcf_file, chroms_bam, chroms_file);
    
    fprintf(stderr, "%s chromosomes to load: %lu\n", label.c_str(), chroms.size());
    
    if (chroms.empty()) {
        fprintf(stderr, "WARNING: No shared chromosomes found for %s, skipping\n", label.c_str());
        return false;
    }
    
    fprintf(stderr, "Creating shared memory segment for %s: %s\n", label.c_str(), shm_name.c_str());
    
    if (!create_shared_vcf(vcf_file, shm_name, chroms, seq2tid, min_vq)) {
        fprintf(stderr, "ERROR: Failed to create shared memory segment for %s\n", label.c_str());
        return false;
    }
    
    active_segments.push_back(shm_name);
    return true;
}

// Destroy all active segments
void cleanup_all_segments() {
    for (auto& seg : active_segments) {
        destroy_shared_vcf(seg);
    }
    active_segments.clear();
}

void help(int code) {
    fprintf(stderr, "vcf_loader_daemon [OPTIONS]\n");
    fprintf(stderr, "Load VCF(s) into shared memory for use by multiple demux_parallel instances.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --vcf -v VCF/BCF file to load (main demux VCF)\n");
    fprintf(stderr, "    --bam -b BAM file (for chromosome/TID mapping)\n");
    fprintf(stderr, "    --name -n Shared memory segment name (base name for all segments)\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --het_vcf -H      Het VCF for ploidy detection (from downsample_vcf_parallel)\n");
    fprintf(stderr, "    --atac_vcf -A     ATAC-demux VCF panel\n");
    fprintf(stderr, "    --atac_het_vcf    ATAC het VCF for ATAC het balance diagnostics\n");
    fprintf(stderr, "    --species_vcf     Species-discrimination VCF panel\n");
    fprintf(stderr, "    --qual -q         Minimum variant quality (default: 50)\n");
    fprintf(stderr, "    --chroms -c       File listing chromosomes to include\n");
    fprintf(stderr, "    --foreground -f   Run in foreground (don't daemonize)\n");
    fprintf(stderr, "    --destroy -d      Destroy existing shared memory segment(s) and exit\n");
    fprintf(stderr, "    --help -h         Display this message and exit\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Shared memory segment naming:\n");
    fprintf(stderr, "  Given --name /myvcf, segments are created as:\n");
    fprintf(stderr, "    Main VCF:      /myvcf\n");
    fprintf(stderr, "    Het VCF:       /myvcf_het\n");
    fprintf(stderr, "    ATAC VCF:      /myvcf_atac\n");
    fprintf(stderr, "    ATAC het VCF:  /myvcf_atac_het\n");
    fprintf(stderr, "    Species VCF:   /myvcf_species\n");
    exit(code);
}

int main(int argc, char* argv[]) {
    
    static struct option long_options[] = {
        {"vcf",            required_argument, 0, 'v'},
        {"het_vcf",        required_argument, 0, 'H'},
        {"atac_vcf",       required_argument, 0, 'A'},
        {"atac_het_vcf",   required_argument, 0, 1001},
        {"species_vcf",    required_argument, 0, 1002},
        {"bam",            required_argument, 0, 'b'},
        {"name",           required_argument, 0, 'n'},
        {"qual",           required_argument, 0, 'q'},
        {"chroms",         required_argument, 0, 'c'},
        {"foreground",     no_argument,       0, 'f'},
        {"destroy",        no_argument,       0, 'd'},
        {"help",           no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };
    
    string vcf_file = "";
    string het_vcf_file = "";
    string atac_vcf_file = "";
    string atac_het_vcf_file = "";
    string species_vcf_file = "";
    string bam_file = "";
    string shm_name = "";
    string chroms_file = "";
    int min_vq = 50;
    bool foreground = false;
    bool destroy_only = false;
    
    int option_index = 0;
    int ch;
    
    if (argc == 1) {
        help(0);
    }
    
    while ((ch = getopt_long(argc, argv, "v:H:A:b:n:q:c:fdh", long_options, &option_index)) != -1) {
        switch (ch) {
            case 'v':
                vcf_file = optarg;
                break;
            case 'H':
                het_vcf_file = optarg;
                break;
            case 'A':
                atac_vcf_file = optarg;
                break;
            case 1001:
                atac_het_vcf_file = optarg;
                break;
            case 1002:
                species_vcf_file = optarg;
                break;
            case 'b':
                bam_file = optarg;
                break;
            case 'n':
                shm_name = optarg;
                break;
            case 'q':
                min_vq = atoi(optarg);
                break;
            case 'c':
                chroms_file = optarg;
                break;
            case 'f':
                foreground = true;
                break;
            case 'd':
                destroy_only = true;
                break;
            case 'h':
                help(0);
                break;
            default:
                help(1);
        }
    }
    
    if (shm_name.empty()) {
        fprintf(stderr, "ERROR: --name/-n is required\n");
        exit(1);
    }
    
    // Ensure name starts with /
    if (shm_name[0] != '/') {
        shm_name = "/" + shm_name;
    }
    
    // Derive all segment names from the base name
    string shm_het         = shm_name + "_het";
    string shm_atac        = shm_name + "_atac";
    string shm_atac_het    = shm_name + "_atac_het";
    string shm_species     = shm_name + "_species";
    
    if (destroy_only) {
        fprintf(stderr, "Destroying shared memory segments...\n");
        // Try to destroy all possible segments (non-existent ones just print a warning)
        destroy_shared_vcf(shm_name);
        destroy_shared_vcf(shm_het);
        destroy_shared_vcf(shm_atac);
        destroy_shared_vcf(shm_atac_het);
        destroy_shared_vcf(shm_species);
        return 0;
    }
    
    if (vcf_file.empty()) {
        fprintf(stderr, "ERROR: --vcf/-v is required\n");
        exit(1);
    }
    
    if (bam_file.empty()) {
        fprintf(stderr, "ERROR: --bam/-b is required for TID mapping\n");
        exit(1);
    }
    
    // Validate: --atac_het_vcf without --atac_vcf is likely a mistake
    if (!atac_het_vcf_file.empty() && atac_vcf_file.empty()) {
        fprintf(stderr, "WARNING: --atac_het_vcf provided without --atac_vcf\n");
    }
    
    // Get chromosome info from BAM
    bam_reader reader;
    reader.set_file(bam_file);
    map<string, int> seq2tid = reader.get_seq2tid();
    
    set<string> chroms_bam;
    for (auto& kv : seq2tid) {
        chroms_bam.insert(kv.first);
    }
    
    // ========================================================================
    // Load main VCF (required)
    // ========================================================================
    if (!load_vcf_segment(vcf_file, shm_name, "Main VCF",
                          chroms_bam, chroms_file, seq2tid, min_vq)) {
        fprintf(stderr, "ERROR: Main VCF loading failed\n");
        exit(1);
    }
    
    // ========================================================================
    // Load optional VCFs. On failure, clean up everything and exit.
    // ========================================================================
    
    struct OptionalVCF {
        const string& file;
        const string& seg_name;
        const char* label;
    };
    
    vector<OptionalVCF> optional_vcfs = {
        { het_vcf_file,      shm_het,      "Het VCF"      },
        { atac_vcf_file,     shm_atac,     "ATAC VCF"     },
        { atac_het_vcf_file, shm_atac_het, "ATAC het VCF" },
        { species_vcf_file,  shm_species,  "Species VCF"  },
    };
    
    for (auto& opt : optional_vcfs) {
        if (opt.file.empty()) continue;
        
        fprintf(stderr, "\n");
        if (!load_vcf_segment(opt.file, opt.seg_name, opt.label,
                              chroms_bam, chroms_file, seq2tid, min_vq)) {
            // load_vcf_segment returns false for both "no shared chroms" (warning)
            // and "create_shared_vcf failed" (error). The warning case already
            // printed, and we can continue. The error case also printed. We treat
            // create_shared_vcf failure as fatal.
            // Check if the VCF has any shared chroms at all to distinguish.
            set<string> check = get_shared_chroms(opt.file, chroms_bam, chroms_file);
            if (!check.empty()) {
                // Had chroms but create_shared_vcf failed: fatal
                fprintf(stderr, "ERROR: %s shared memory creation failed, aborting\n", opt.label);
                cleanup_all_segments();
                exit(1);
            }
            // Otherwise: no shared chroms, already warned, continue
        }
    }
    
    // ========================================================================
    // Report loaded segments
    // ========================================================================
    fprintf(stderr, "\n");
    fprintf(stderr, "Loaded %lu shared memory segment(s):\n", active_segments.size());
    for (auto& seg : active_segments) {
        fprintf(stderr, "  %s\n", seg.c_str());
    }
    
    // ========================================================================
    // Build demux_parallel usage hint
    // ========================================================================
    string usage_hint = "  demux_parallel --shared_vcf " + shm_name;
    for (auto& seg : active_segments) {
        if (seg == shm_name) continue;
        if (seg == shm_het)      usage_hint += " --shared_het_vcf " + seg;
        if (seg == shm_atac)     usage_hint += " --atac_shared_vcf " + seg;
        if (seg == shm_atac_het) usage_hint += " --atac_shared_het_vcf " + seg;
        if (seg == shm_species)  usage_hint += " --species_shared_vcf " + seg;
    }
    usage_hint += " ...";
    
    // Set up signal handlers
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);
    
    if (!foreground) {
        // Daemonize
        pid_t pid = fork();
        if (pid < 0) {
            perror("fork");
            cleanup_all_segments();
            exit(1);
        }
        if (pid > 0) {
            // Parent: print info and exit
            fprintf(stderr, "\nDaemon started with PID %d\n", pid);
            fprintf(stderr, "\nTo use with demux_parallel:\n");
            fprintf(stderr, "%s\n", usage_hint.c_str());
            fprintf(stderr, "\nTo destroy: vcf_loader_daemon --destroy --name %s\n", shm_name.c_str());
            exit(0);
        }
        
        // Child: become daemon
        setsid();
        close(STDIN_FILENO);
        close(STDOUT_FILENO);
        close(STDERR_FILENO);
    }
    else {
        fprintf(stderr, "\nRunning in foreground. Press Ctrl+C to stop.\n");
        fprintf(stderr, "\nTo use with demux_parallel:\n");
        fprintf(stderr, "%s\n", usage_hint.c_str());
    }
    
    // Keep running until signaled
    while (keep_running) {
        sleep(1);
    }
    
    // Cleanup
    if (foreground) {
        fprintf(stderr, "\nShutting down...\n");
    }
    cleanup_all_segments();
    
    return 0;
}
