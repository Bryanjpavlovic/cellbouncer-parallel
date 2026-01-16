#include <getopt.h>
#include <string>
#include <iostream>
#include <fstream>
#include <set>
#include <map>
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
static string shm_name_global;
static string shm_name_het_global;

void signal_handler(int sig) {
    keep_running = 0;
}

void help(int code) {
    fprintf(stderr, "vcf_loader_daemon [OPTIONS]\n");
    fprintf(stderr, "Load VCF(s) into shared memory for use by multiple demux_parallel instances.\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "[OPTIONS]:\n");
    fprintf(stderr, "===== REQUIRED =====\n");
    fprintf(stderr, "    --vcf -v VCF/BCF file to load (main demux VCF)\n");
    fprintf(stderr, "    --bam -b BAM file (for chromosome/TID mapping)\n");
    fprintf(stderr, "    --name -n Shared memory segment name\n");
    fprintf(stderr, "===== OPTIONAL =====\n");
    fprintf(stderr, "    --het_vcf -H Het VCF file for ploidy detection (from downsample_vcf_parallel)\n");
    fprintf(stderr, "    --qual -q Minimum variant quality (default: 50)\n");
    fprintf(stderr, "    --chroms -c File listing chromosomes to include\n");
    fprintf(stderr, "    --foreground -f Run in foreground (don't daemonize)\n");
    fprintf(stderr, "    --destroy -d Destroy existing shared memory segment(s) and exit\n");
    fprintf(stderr, "    --help -h Display this message and exit\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "When --het_vcf is provided, a second shared memory segment is created\n");
    fprintf(stderr, "with the suffix '_het' (e.g., /myvcf becomes /myvcf and /myvcf_het).\n");
    exit(code);
}

int main(int argc, char* argv[]) {
    
    static struct option long_options[] = {
        {"vcf", required_argument, 0, 'v'},
        {"het_vcf", required_argument, 0, 'H'},
        {"bam", required_argument, 0, 'b'},
        {"name", required_argument, 0, 'n'},
        {"qual", required_argument, 0, 'q'},
        {"chroms", required_argument, 0, 'c'},
        {"foreground", no_argument, 0, 'f'},
        {"destroy", no_argument, 0, 'd'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    string vcf_file = "";
    string het_vcf_file = "";
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
    
    while ((ch = getopt_long(argc, argv, "v:H:b:n:q:c:fdh", long_options, &option_index)) != -1) {
        switch (ch) {
            case 'v':
                vcf_file = optarg;
                break;
            case 'H':
                het_vcf_file = optarg;
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
    
    shm_name_global = shm_name;
    string shm_name_het = shm_name + "_het";
    shm_name_het_global = shm_name_het;
    
    if (destroy_only) {
        fprintf(stderr, "Destroying shared memory segments...\n");
        destroy_shared_vcf(shm_name);
        if (!het_vcf_file.empty()) {
            destroy_shared_vcf(shm_name_het);
        }
        // Also try to destroy het segment even if not specified (cleanup)
        // This will just print an error if it doesn't exist, which is fine
        destroy_shared_vcf(shm_name_het);
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
    
    // Get chromosome info from BAM
    bam_reader reader;
    reader.set_file(bam_file);
    map<string, int> seq2tid = reader.get_seq2tid();
    
    set<string> chroms_bam;
    for (auto& kv : seq2tid) {
        chroms_bam.insert(kv.first);
    }
    
    // Get chromosomes from VCF
    set<string> chroms_vcf;
    get_vcf_chroms(vcf_file, chroms_vcf);
    
    // Find shared chromosomes
    set<string> chroms_to_include;
    for (auto& c : chroms_bam) {
        if (chroms_vcf.find(c) != chroms_vcf.end()) {
            chroms_to_include.insert(c);
        }
    }
    
    // Apply user filter if provided
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
        for (auto& c : chroms_to_include) {
            if (user_chroms.find(c) != user_chroms.end()) {
                filtered.insert(c);
            }
        }
        chroms_to_include = filtered;
    }
    
    fprintf(stderr, "Chromosomes to load: %lu\n", chroms_to_include.size());
    
    if (chroms_to_include.empty()) {
        fprintf(stderr, "ERROR: No shared chromosomes found between BAM and VCF\n");
        exit(1);
    }
    
    // Create shared memory segment for main VCF
    fprintf(stderr, "Creating shared memory segment for main VCF: %s\n", shm_name.c_str());
    
    if (!create_shared_vcf(vcf_file, shm_name, chroms_to_include, seq2tid, min_vq)) {
        fprintf(stderr, "ERROR: Failed to create shared memory segment for main VCF\n");
        exit(1);
    }
    
    // Create shared memory segment for het VCF if provided
    bool het_loaded = false;
    if (!het_vcf_file.empty()) {
        fprintf(stderr, "\nCreating shared memory segment for het VCF: %s\n", shm_name_het.c_str());
        
        // Get chromosomes from het VCF
        set<string> chroms_het_vcf;
        get_vcf_chroms(het_vcf_file, chroms_het_vcf);
        
        // Find shared chromosomes (intersection of BAM and het VCF)
        set<string> chroms_to_include_het;
        for (auto& c : chroms_bam) {
            if (chroms_het_vcf.find(c) != chroms_het_vcf.end()) {
                chroms_to_include_het.insert(c);
            }
        }
        
        // Apply same user filter if provided
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
            for (auto& c : chroms_to_include_het) {
                if (user_chroms.find(c) != user_chroms.end()) {
                    filtered.insert(c);
                }
            }
            chroms_to_include_het = filtered;
        }
        
        fprintf(stderr, "Het VCF chromosomes to load: %lu\n", chroms_to_include_het.size());
        
        if (chroms_to_include_het.empty()) {
            fprintf(stderr, "WARNING: No shared chromosomes found for het VCF, skipping\n");
        }
        else {
            if (!create_shared_vcf(het_vcf_file, shm_name_het, chroms_to_include_het, seq2tid, min_vq)) {
                fprintf(stderr, "ERROR: Failed to create shared memory segment for het VCF\n");
                // Clean up main VCF segment
                destroy_shared_vcf(shm_name);
                exit(1);
            }
            het_loaded = true;
        }
    }
    
    // Set up signal handlers
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);
    
    if (!foreground) {
        // Daemonize
        pid_t pid = fork();
        if (pid < 0) {
            perror("fork");
            destroy_shared_vcf(shm_name);
            if (het_loaded) {
                destroy_shared_vcf(shm_name_het);
            }
            exit(1);
        }
        if (pid > 0) {
            // Parent - print info and exit
            fprintf(stderr, "\nDaemon started with PID %d\n", pid);
            fprintf(stderr, "Shared memory available as:\n");
            fprintf(stderr, "  Main VCF: %s\n", shm_name.c_str());
            if (het_loaded) {
                fprintf(stderr, "  Het VCF:  %s\n", shm_name_het.c_str());
            }
            fprintf(stderr, "\nTo use with demux_parallel:\n");
            fprintf(stderr, "  demux_parallel --shared_vcf %s", shm_name.c_str());
            if (het_loaded) {
                fprintf(stderr, " --shared_het_vcf %s", shm_name_het.c_str());
            }
            fprintf(stderr, " ...\n");
            fprintf(stderr, "\nTo destroy: vcf_loader_daemon --destroy --name %s\n", shm_name.c_str());
            exit(0);
        }
        
        // Child - become daemon
        setsid();
        
        // Close standard file descriptors
        close(STDIN_FILENO);
        close(STDOUT_FILENO);
        close(STDERR_FILENO);
    }
    else {
        fprintf(stderr, "\nRunning in foreground. Press Ctrl+C to stop.\n");
        fprintf(stderr, "Shared memory available as:\n");
        fprintf(stderr, "  Main VCF: %s\n", shm_name.c_str());
        if (het_loaded) {
            fprintf(stderr, "  Het VCF:  %s\n", shm_name_het.c_str());
        }
    }
    
    // Keep running until signaled
    while (keep_running) {
        sleep(1);
    }
    
    // Cleanup
    if (foreground) {
        fprintf(stderr, "\nShutting down...\n");
    }
    destroy_shared_vcf(shm_name);
    if (het_loaded) {
        destroy_shared_vcf(shm_name_het);
    }
    
    return 0;
}
