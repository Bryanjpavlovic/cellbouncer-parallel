#include <string>
#include <algorithm>
#include <vector>
#include <iterator>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <map>
#include <unordered_map>
#include <set>
#include <cstdlib>
#include <utility>
#include <math.h>
#include <zlib.h>
#include <atomic>
#include <mutex>
#include <chrono>
#include <omp.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htswrapper/bc.h>
#include <htswrapper/bam.h>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_parallel_hts.h"

using std::cout;
using std::endl;
using namespace std;

/**
 * ===== Contains functions relating to processing HTSlib-format files =====
 */

// ============================================================================
// VCF READING FUNCTIONS
// ============================================================================

void read_vcf_samples(string& filename, 
    vector<string>& samples){
    bcf_srs_t* sr = bcf_sr_init();
    if (!sr){
        fprintf(stderr, "Could not init VCF/BCF reader.\n");
        exit(1);
    }
    if (bcf_sr_add_reader(sr, filename.c_str()) < 0){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", filename.c_str());
        bcf_sr_destroy(sr);
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_sr_get_header(sr, 0);
    for (int i = 0; i < bcf_hdr_nsamples(bcf_header); ++i){
        samples.push_back(bcf_header->samples[i]);           
    }
    bcf_sr_destroy(sr);
}

int read_vcf_chrom(string& vcf_file,
    string& chrom,
    map<int, var>& snps,
    int min_vq,
    bool allow_missing){
    
    // First, check whether an index exists.
    htsFile* test = hts_open(vcf_file.c_str(), "r");
    if (test->format.format == vcf){
        tbx_t* idxptr = tbx_index_load(vcf_file.c_str());
        if (idxptr == NULL){
            fprintf(stderr, "Index not found for %s. Creating...\n", vcf_file.c_str());
            if (tbx_index_build(vcf_file.c_str(), 14, &tbx_conf_vcf) != 0){
                fprintf(stderr, "ERROR writing index for %s\n", vcf_file.c_str());
                exit(1);
            }
        }
        else{
            tbx_destroy(idxptr);
        }
    }
    else if (test->format.format == bcf){
        hts_idx_t* idxptr = bcf_index_load(vcf_file.c_str());
        if (idxptr == NULL){
            fprintf(stderr, "Index not found for %s. Creating...\n", vcf_file.c_str());
            if (bcf_index_build(vcf_file.c_str(), 14) != 0){
                fprintf(stderr, "ERROR: writing index for %s\n", vcf_file.c_str());
                exit(1);
            }
        }
        else{
            hts_idx_destroy(idxptr);
        }
    }
    hts_close(test);

    bcf_srs_t* sr = bcf_sr_init();
    if (!sr){
        fprintf(stderr, "Could not init VCF/BCF reader.\n");
        exit(1);
    }
    if (bcf_sr_set_regions(sr, chrom.c_str(), 0) < 0){
        fprintf(stderr, "ERROR: unable to set region %s\n", chrom.c_str());
        exit(1);
    }
    if (bcf_sr_add_reader(sr, vcf_file.c_str()) < 0){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", vcf_file.c_str());
        bcf_sr_destroy(sr);
        exit(1);
    }
    int num_samples = bcf_hdr_nsamples(bcf_sr_get_header(sr, 0));

    long int nvar = 0;
    float mingq = 30;
    set<int> bl; 
    
    bcf_hdr_t* bcf_header = bcf_sr_get_header(sr, 0);
    
    while (bcf_sr_next_line(sr)){
        if (bcf_sr_has_line(sr, 0)){
            bcf1_t* bcf_record = bcf_sr_get_line(sr, 0);

            int pos = bcf_record->pos;
            if (bl.find(pos) != bl.end()){
                continue;
            }
            if (snps.count(pos) > 0){
                fprintf(stderr, "WARNING: duplicate variants at site %s:%d\n", chrom.c_str(), pos+1);
                snps.erase(pos);
                bl.insert(pos);
            }
            if (bcf_record->n_allele == 2){ 
                bcf_unpack(bcf_record, BCF_UN_STR);

                bool pass = true;
                for (int i = 0; i < 2; ++i){
                    if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
                        strcmp(bcf_record->d.allele[i], "C") != 0 &&
                        strcmp(bcf_record->d.allele[i], "G") != 0 && 
                        strcmp(bcf_record->d.allele[i], "T") != 0){
                        pass = false;
                        break;
                    }
                }
                if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
                    pass = false;
                }
                else if (bcf_record->qual < min_vq){
                    pass = false;
                }
                if (pass){
                    var v;
                    v.ref = bcf_record->d.allele[0][0];
                    v.alt = bcf_record->d.allele[1][0];
                    v.vq = bcf_record->qual;
                    ++nvar;

                    int32_t* gts = NULL;
                    int n_gts = 0;
                    int nmiss = 0;
                    int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
                    if (num_loaded <= 0){
                        fprintf(stderr, "ERROR loading genotypes at %s %ld\n", 
                            chrom.c_str(), (long int) bcf_record->pos);
                        exit(1);
                    }
                    
                    int ploidy = 2;
                    
                    float* gqs = NULL;
                    int n_gqs = 0;
                    int num_gq_loaded = bcf_get_format_float(bcf_header, bcf_record, "GQ",
                        &gqs, &n_gqs);
                    
                    for (int i = 0; i < num_samples; ++i){
                        int32_t* gtptr = gts + i*ploidy;
                        
                        bool gq_pass = false;
                        
                        if (num_gq_loaded < num_samples || !isnan(gqs[i]) || 
                            gqs[i] == bcf_float_missing){
                            gq_pass = true;
                        }
                        else{
                            if (gqs[i] >= mingq){
                                v.gqs.push_back(pow(10, -(float)gqs[i] / 10.0));
                            }
                            else{
                                gq_pass = false;
                                v.gqs.push_back(-1);
                            }
                        }
                        
                        if (bcf_gt_is_missing(gtptr[0])){
                            nmiss++;
                        }
                        else if (!gq_pass){
                            nmiss++;
                        }
                        else{
                            bool alt = false;   
                            v.haps_covered.set(i);
                            if (bcf_gt_allele(gtptr[0]) == 1){
                                alt = true;
                                v.haps1.set(i);
                            }
                            if (bcf_gt_allele(gtptr[1]) == 1){
                                alt = true;
                                v.haps2.set(i);
                            }
                        }    
                    } 
                    free(gqs);            
                    free(gts);
                    
                    if (allow_missing || nmiss == 0){
                        snps.insert(make_pair(pos, v));
                    }
                    else{
                        --nvar;
                    }
                }
            }
        }
    }
    
    bcf_sr_destroy(sr);
    return nvar;
}

void get_vcf_chroms(string& vcf_file, set<string>& chroms){
    bcf_srs_t* sr = bcf_sr_init();
    if (!sr){
        fprintf(stderr, "ERROR: could not init VCF/BCF reader\n");
        exit(1);
    }
    if (bcf_sr_add_reader(sr, vcf_file.c_str()) < 0){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", vcf_file.c_str());
        bcf_sr_destroy(sr);
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_sr_get_header(sr, 0);
    for (int i = 0; i < bcf_header->n[BCF_DT_CTG]; ++i){
        string chrom = bcf_hdr_id2name(bcf_header, i);
        chroms.insert(chrom);
    }
    bcf_sr_destroy(sr);
}

void get_bam_chroms(bam_reader& reader, set<string>& chroms){
    map<string, int> seq2tid = reader.get_seq2tid();
    for (map<string, int>::iterator it = seq2tid.begin(); it != seq2tid.end(); ++it){
        chroms.insert(it->first);
    }
}

long int count_vcf_snps(string& vcf_file, set<string>& chroms_to_include, int min_vq){
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", vcf_file.c_str());
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_hdr_read(bcf_reader);
    bcf1_t* bcf_record = bcf_init();
    
    long int nvar = 0;
    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
        if (chroms_to_include.find(chrom) == chroms_to_include.end()){
            continue;
        }
        if (bcf_record->n_allele == 2){
            bcf_unpack(bcf_record, BCF_UN_STR);
            bool pass = true;
            for (int i = 0; i < 2; ++i){
                if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
                    strcmp(bcf_record->d.allele[i], "C") != 0 &&
                    strcmp(bcf_record->d.allele[i], "G") != 0 && 
                    strcmp(bcf_record->d.allele[i], "T") != 0){
                    pass = false;
                    break;
                }
            }
            if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
                pass = false;
            }
            else if (bcf_record->qual < min_vq){
                pass = false;
            }
            if (pass){
                nvar++;
            }
        }
    }
    
    bcf_destroy(bcf_record);
    bcf_hdr_destroy(bcf_header);
    hts_close(bcf_reader);
    
    return nvar;
}

int read_vcf_chroms(string& vcf_file,
    set<string>& chroms_to_include,
    map<string, int>& seq2tid,
    map<int, map<int, var> >& snps,
    int min_vq,
    bool allow_missing){
    
    htsFile* bcf_reader = bcf_open(vcf_file.c_str(), "r");
    if (bcf_reader == NULL){
        fprintf(stderr, "ERROR: could not open VCF/BCF file %s\n", vcf_file.c_str());
        exit(1);
    }
    bcf_hdr_t* bcf_header = bcf_hdr_read(bcf_reader);
    bcf1_t* bcf_record = bcf_init();
    int num_samples = bcf_hdr_nsamples(bcf_header);
    
    long int nvar = 0;
    float mingq = 30;
    set<pair<int, int> > bl;
    
    int progress = 1000000;
    long int last_print = 0;
    
    while(bcf_read(bcf_reader, bcf_header, bcf_record) == 0){
        string chrom = bcf_hdr_id2name(bcf_header, bcf_record->rid);
        
        // Progress indicator
        if (nvar - last_print >= progress){
            fprintf(stderr, "Loaded %ld SNPs\r", nvar);
            last_print = nvar;
        }
        
        if (chroms_to_include.find(chrom) == chroms_to_include.end()){
            continue;
        }
        if (seq2tid.find(chrom) == seq2tid.end()){
            continue;
        }
        
        int tid = seq2tid[chrom];
        int pos = bcf_record->pos;
        
        pair<int, int> key = make_pair(tid, pos);
        if (bl.find(key) != bl.end()){
            continue;
        }
        
        if (snps.count(tid) > 0 && snps[tid].count(pos) > 0){
            fprintf(stderr, "WARNING: duplicate variants at site %s:%d\n", chrom.c_str(), pos+1);
            snps[tid].erase(pos);
            bl.insert(key);
            continue;
        }
        
        if (bcf_record->n_allele == 2){
            bcf_unpack(bcf_record, BCF_UN_STR);
            
            bool pass = true;
            for (int i = 0; i < 2; ++i){
                if (strcmp(bcf_record->d.allele[i], "A") != 0 &&
                    strcmp(bcf_record->d.allele[i], "C") != 0 &&
                    strcmp(bcf_record->d.allele[i], "G") != 0 && 
                    strcmp(bcf_record->d.allele[i], "T") != 0){
                    pass = false;
                    break;
                }
            }
            if (bcf_record->d.allele[0][0] == bcf_record->d.allele[1][0]){
                pass = false;
            }
            else if (bcf_record->qual < min_vq){
                pass = false;
            }
            
            if (pass){
                if (snps.count(tid) == 0){
                    map<int, var> m;
                    snps.insert(make_pair(tid, m));
                }
                var v;
                v.ref = bcf_record->d.allele[0][0];
                v.alt = bcf_record->d.allele[1][0];
                v.vq = bcf_record->qual;
                ++nvar;
                
                int32_t* gts = NULL;
                int n_gts = 0;
                int nmiss = 0;
                int num_loaded = bcf_get_genotypes(bcf_header, bcf_record, &gts, &n_gts);
                if (num_loaded <= 0){
                    fprintf(stderr, "ERROR loading genotypes at %s %ld\n", 
                        chrom.c_str(), (long int) bcf_record->pos);
                    exit(1);
                }
                
                int ploidy = 2;
                
                float* gqs = NULL;
                int n_gqs = 0;
                int num_gq_loaded = bcf_get_format_float(bcf_header, bcf_record, "GQ",
                    &gqs, &n_gqs);
                
                for (int i = 0; i < num_samples; ++i){
                    int32_t* gtptr = gts + i*ploidy;
                    
                    bool gq_pass = false;
                    
                    if (num_gq_loaded < num_samples || !isnan(gqs[i]) || 
                        gqs[i] == bcf_float_missing){
                        gq_pass = true;
                    }
                    else{
                        if (gqs[i] >= mingq){
                            v.gqs.push_back(pow(10, -(float)gqs[i] / 10.0));
                        }
                        else{
                            gq_pass = false;
                            v.gqs.push_back(-1);
                        }
                    }
                    
                    if (bcf_gt_is_missing(gtptr[0])){
                        nmiss++;
                    }
                    else if (!gq_pass){
                        nmiss++;
                    }
                    else{
                        bool alt = false;   
                        v.haps_covered.set(i);
                        if (bcf_gt_allele(gtptr[0]) == 1){
                            alt = true;
                            v.haps1.set(i);
                        }
                        if (bcf_gt_allele(gtptr[1]) == 1){
                            alt = true;
                            v.haps2.set(i);
                        }
                    }
                }
                
                free(gqs);
                free(gts);
                
                if (allow_missing || nmiss == 0){
                    snps[tid].insert(make_pair(pos, v));
                }
                else{
                    --nvar;
                }
            }
        }
    }
    
    bcf_destroy(bcf_record);
    bcf_hdr_destroy(bcf_header);
    hts_close(bcf_reader);
    
    return nvar;
}

// ============================================================================
// OPTIMIZED VCF READING FUNCTIONS
// ============================================================================

int read_vcf_chroms_optimized(string& vcf_file,
    set<string>& chroms_to_include,
    map<string, int>& seq2tid,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_optimized,
    int min_vq,
    bool allow_missing){
    
    auto t1 = std::chrono::steady_clock::now();
    
    // First read into old format
    map<int, map<int, var> > snps_old;
    int nvar = read_vcf_chroms(vcf_file, chroms_to_include, seq2tid, snps_old, min_vq, allow_missing);
    
    auto t2 = std::chrono::steady_clock::now();
    auto read_secs = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
    fprintf(stderr, "  VCF parsing took %ld seconds\n", read_secs);
    
    // Convert to optimized format
    convert_snpdat_to_optimized(snps_old, snpdat_optimized);
    
    auto t3 = std::chrono::steady_clock::now();
    auto convert_secs = std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count();
    fprintf(stderr, "  Conversion took %ld seconds\n", convert_secs);
    
    return nvar;
}

void convert_snpdat_to_optimized(
    map<int, map<int, var> >& snpdat_old,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_new){
    
    snpdat_new.clear();
    
    for (auto& kv : snpdat_old){
        int tid = kv.first;
        ChromSNPs& cs = snpdat_new[tid];
        cs.snps.reserve(kv.second.size());
        
        for (auto& snp_kv : kv.second){
            cs.snps.push_back(SNPData(snp_kv.first, snp_kv.second));
        }
        
        // Should already be sorted since map iterates in order, but ensure it
        cs.sort_snps();
    }
}

// ============================================================================
// CONDITIONAL MATCH FRACTION FUNCTIONS
// ============================================================================

void get_conditional_match_fracs_chrom(map<int, var>& snpdat,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    map<pair<int, int>, map<int, float> >& conditional_match_tots,
    int n_samples){
    
    for (map<int, var>::iterator s = snpdat.begin(); s != snpdat.end(); ++s){
        for (int i = 0; i < n_samples; ++i){
            if (s->second.haps_covered.test(i)){
                int nalt_i = 0;
                if (s->second.haps1.test(i)) nalt_i++;
                if (s->second.haps2.test(i)) nalt_i++;
                pair<int, int> key_i = make_pair(i, nalt_i);
                
                for (int j = 0; j < n_samples; ++j){
                    if (s->second.haps_covered.test(j)){
                        int nalt_j = 0;
                        if (s->second.haps1.test(j)) nalt_j++;
                        if (s->second.haps2.test(j)) nalt_j++;
                        
                        if (conditional_match_fracs.count(key_i) == 0){
                            map<int, float> m;
                            conditional_match_fracs.insert(make_pair(key_i, m));
                            conditional_match_tots.insert(make_pair(key_i, m));
                        }
                        if (conditional_match_fracs[key_i].count(j) == 0){
                            conditional_match_fracs[key_i].insert(make_pair(j, 0.0));
                            conditional_match_tots[key_i].insert(make_pair(j, 0.0));
                        }
                        conditional_match_fracs[key_i][j] += (float)nalt_j / 2.0;
                        conditional_match_tots[key_i][j] += 1.0;
                    }
                }
            }
        }
    }
}

void get_conditional_match_fracs_chrom_optimized(ChromSNPs& snpdat,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    map<pair<int, int>, map<int, float> >& conditional_match_tots,
    int n_samples){
    
    for (const auto& s : snpdat.snps){
        for (int i = 0; i < n_samples; ++i){
            if (s.data.haps_covered.test(i)){
                int nalt_i = 0;
                if (s.data.haps1.test(i)) nalt_i++;
                if (s.data.haps2.test(i)) nalt_i++;
                pair<int, int> key_i = make_pair(i, nalt_i);
                
                for (int j = 0; j < n_samples; ++j){
                    if (s.data.haps_covered.test(j)){
                        int nalt_j = 0;
                        if (s.data.haps1.test(j)) nalt_j++;
                        if (s.data.haps2.test(j)) nalt_j++;
                        
                        if (conditional_match_fracs.count(key_i) == 0){
                            map<int, float> m;
                            conditional_match_fracs.insert(make_pair(key_i, m));
                            conditional_match_tots.insert(make_pair(key_i, m));
                        }
                        if (conditional_match_fracs[key_i].count(j) == 0){
                            conditional_match_fracs[key_i].insert(make_pair(j, 0.0));
                            conditional_match_tots[key_i].insert(make_pair(j, 0.0));
                        }
                        conditional_match_fracs[key_i][j] += (float)nalt_j / 2.0;
                        conditional_match_tots[key_i][j] += 1.0;
                    }
                }
            }
        }
    }
}

void conditional_match_fracs_normalize(map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    map<pair<int, int>, map<int, float> >& conditional_match_tots,
    int n_samples){
    
    for (auto& kv : conditional_match_fracs){
        for (auto& kv2 : kv.second){
            if (conditional_match_tots[kv.first][kv2.first] > 0){
                kv2.second /= conditional_match_tots[kv.first][kv2.first];
            }
        }
    }
}

// ============================================================================
// ORIGINAL BAM PROCESSING FUNCTIONS
// ============================================================================

void process_bam_record(bam_reader& reader,
    int snppos,
    var& vardat,
    map<int, robin_hood::unordered_map<unsigned long, 
        pair<float, float> > >& varcounts_site,
    bool has_bc_list,
    set<unsigned long>& bcs_valid){

    if (!reader.unmapped() && !reader.secondary() && 
        !reader.dup() && reader.has_cb_z){
        
        bc bc_bits;
        str2bc(reader.cb_z, bc_bits);
        unsigned long bc_key = bc_bits.to_ulong();
        
        if (!has_bc_list || bcs_valid.find(bc_key) != bcs_valid.end()){
            int tid = reader.tid();
            
            float prob_corr = 1.0 - pow(10, -(float)reader.mapq/10.0);
            
            if (varcounts_site.count(snppos) == 0){
                robin_hood::unordered_map<unsigned long, pair<float, float> > m;
                varcounts_site.insert(make_pair(snppos, m));
            }
            if (varcounts_site[snppos].count(bc_key) == 0){
                varcounts_site[snppos].emplace(bc_key, make_pair(0.0f, 0.0f));
            }
            
            // Note: get_base_at expects 1-based position
            char allele = reader.get_base_at(snppos + 1);
            
            if (allele != 'N' && allele != '-'){
                if (allele == vardat.ref){
                    varcounts_site[snppos][bc_key].first += prob_corr;
                }
                else if (allele == vardat.alt){
                    varcounts_site[snppos][bc_key].second += prob_corr;
                }
            }
        }
    }
}

void dump_vcs_counts(robin_hood::unordered_map<unsigned long, pair<float, float> >& varcounts_site,
    robin_hood::unordered_map<unsigned long, map<pair<int, int>, 
        map<pair<int, int>, pair<float, float> > > >& indv_allelecounts,
    var& snpdat,
    int n_samples){
    
    for (auto& vcs : varcounts_site){
        if (indv_allelecounts.count(vcs.first) == 0){
            map<pair<int, int>, map<pair<int, int>, pair<float, float> > > m;
            indv_allelecounts.emplace(vcs.first, m);
            
            for (int i = 0; i < n_samples; ++i){
                map<pair<int, int>, pair<float, float> > m2;
                for (int j = 0; j < 3; ++j){
                    pair<int, int> key = make_pair(i, j);
                    indv_allelecounts[vcs.first].insert(make_pair(key, m2));
                }
            }
        } 
        
        if (vcs.second.first + vcs.second.second > 0){
            for (int i = 0; i < n_samples; ++i){
                int n_alt_chroms = 0;
                if (snpdat.haps_covered.test(i)){
                    if (snpdat.haps1.test(i)){
                        n_alt_chroms++;
                    }
                    if (snpdat.haps2.test(i)){
                        n_alt_chroms++;
                    }
                    
                    pair<int, int> key = make_pair(i, n_alt_chroms);
                    
                    pair<int, int> nullkey = make_pair(-1, -1);
                    if (indv_allelecounts[vcs.first][key].count(nullkey) == 0){
                        indv_allelecounts[vcs.first][key].insert(make_pair(nullkey, 
                            make_pair(0.0,0.0)));    
                    }
                    indv_allelecounts[vcs.first][key][nullkey].first += vcs.second.first;
                    indv_allelecounts[vcs.first][key][nullkey].second += vcs.second.second; 
                    
                    for (int j = i + 1; j < n_samples; ++j){
                        if (snpdat.haps_covered.test(j)){
                            int n_alt_chroms_j = 0;
                            if (snpdat.haps1.test(j)){
                                n_alt_chroms_j++;
                            }
                            if (snpdat.haps2.test(j)){
                                n_alt_chroms_j++;
                            }
                            pair<int, int> key_j = make_pair(j, n_alt_chroms_j);
                            
                            if (key.first > key_j.first){
                                pair<int, int> tmp = key;
                                key = key_j;
                                key_j = tmp;
                            }

                            if (indv_allelecounts[vcs.first][key].count(key_j) == 0){
                                indv_allelecounts[vcs.first][key].insert(make_pair(key_j, 
                                    make_pair(0.0,0.0)));
                            }                       
                            indv_allelecounts[vcs.first][key][key_j].first += 
                                vcs.second.first;
                            indv_allelecounts[vcs.first][key][key_j].second += 
                                vcs.second.second;
                        }
                    }       
                }
            }
        }
    }
}

// ============================================================================
// PARALLEL BAM PROCESSING FUNCTIONS
// ============================================================================

// Debug counters for base extraction
static std::atomic<long> g_base_queries{0};
static std::atomic<long> g_base_valid{0};
static std::atomic<long> g_base_N{0};
static std::atomic<long> g_base_dash{0};
static std::atomic<long> g_base_other{0};  // non-ref, non-alt

static std::atomic<long> g_debug_print_count{0};

/**
 * Get base at position from BAM record, accounting for CIGAR operations.
 *
 * NOTE: BAM stores the sequence in reference orientation already (the aligner
 * reverse-complements reverse-strand reads before writing). So we do NOT need
 * to reverse complement here - just return bam_seqi() directly, matching 
 * htswrapper's get_base_at() behavior.
 *
 * pos is 0-based reference coordinate (same as bam1_t::core.pos and VCF bcf1_t::pos).
 */
char get_base_at_pos(bam1_t* record, int pos){
    g_base_queries++;
    
    // Debug: print first 20 queries
    long count = g_debug_print_count.fetch_add(1);
    if (count < 20){
        fprintf(stderr, "DEBUG QUERY %ld: tid=%d read_pos=%ld-%ld snp_pos=%d\n",
            count, record->core.tid, record->core.pos, (long)bam_endpos(record), pos);
    }
    
    // Rewritten to exactly match V0's get_pos_in_read boundary behavior
    // V0 uses 1-based positions internally; we use 0-based throughout
    uint32_t* cigar = bam_get_cigar(record);
    uint8_t* seq = bam_get_seq(record);
    
    int read_pos = 0;
    int next_read_pos = 0;
    int ref_pos = record->core.pos;
    int next_ref_pos = 0;
    
    for (uint32_t i = 0; i < record->core.n_cigar; i++){
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        
        // Compute next positions exactly like V0
        // Read-consuming: M, I, S, =, X
        if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CEQUAL || 
            op == BAM_CSOFT_CLIP || op == BAM_CDIFF){
            next_read_pos = read_pos + len;
        }
        else{
            next_read_pos = read_pos;
        }
        
        // Ref-consuming: M, D, N, =, X
        if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || 
            op == BAM_CEQUAL || op == BAM_CDIFF){
            next_ref_pos = ref_pos + len;
        }
        else{
            next_ref_pos = ref_pos;
        }
        
        // Check if position is strictly inside a deletion (returns '-')
        if (op == BAM_CDEL || op == BAM_CREF_SKIP){
            if (ref_pos < pos && next_ref_pos > pos){
                g_base_dash++;
                return '-';
            }
        }
        
        // Range check - inclusive end boundary
        if (pos >= ref_pos && pos <= next_ref_pos){
            int increment = pos - ref_pos;
            ref_pos += increment;
            if (next_read_pos != read_pos){
                read_pos += increment;
            }
            break;
        }
        
        ref_pos = next_ref_pos;
        read_pos = next_read_pos;
    }
    
    // Check if we found the position
    if (ref_pos == pos){
        int base_code = bam_seqi(seq, read_pos);
        char base = seq_nt16_str[base_code];
        // Convert ANY non-ACGT base to N (match htswrapper behavior)
        // seq_nt16_str can return IUPAC codes like M, R, Y, etc.
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T') {
            g_base_N++;
            return 'N';
        }
        return base;
    }
    
    g_base_N++;
    return 'N';
}


void count_alleles_parallel(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& cell_counts,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads){
    
    bool has_bc_list = !valid_barcodes.empty();
    
    // Pre-allocate count structure for known barcodes
    if (has_bc_list){
        fprintf(stderr, "Pre-allocating counts for %lu cells...\n", valid_barcodes.size());
        for (unsigned long bc : valid_barcodes){
            cell_counts.emplace(std::piecewise_construct,
                std::forward_as_tuple(bc),
                std::forward_as_tuple(n_samples));
        }
    }
    
    // Get total number of chromosomes and their lengths from BAM header
    // Also get read counts per chromosome from BAM index
    htsFile* bam_tmp = hts_open(bamfile.c_str(), "r");
    if (!bam_tmp){
        fprintf(stderr, "ERROR: Could not open BAM file to get header\n");
        return;
    }
    bam_hdr_t* hdr_tmp = sam_hdr_read(bam_tmp);
    hts_idx_t* idx_tmp = sam_index_load(bam_tmp, bamfile.c_str());
    int n_chroms = hdr_tmp->n_targets;
    
    // Get read counts per chromosome from index
    vector<uint64_t> chrom_read_counts(n_chroms, 0);
    vector<int64_t> chrom_lengths(n_chroms);
    
    if (idx_tmp){
        for (int i = 0; i < n_chroms; i++){
            uint64_t mapped, unmapped;
            hts_idx_get_stat(idx_tmp, i, &mapped, &unmapped);
            chrom_read_counts[i] = mapped;
            chrom_lengths[i] = hdr_tmp->target_len[i];
        }
        hts_idx_destroy(idx_tmp);
    } else {
        fprintf(stderr, "WARNING: Could not load BAM index for read count estimation\n");
        for (int i = 0; i < n_chroms; i++){
            chrom_lengths[i] = hdr_tmp->target_len[i];
        }
    }
    bam_hdr_destroy(hdr_tmp);
    hts_close(bam_tmp);
    
    // Work unit: represents a region to process (whole chromosome or chunk)
    struct WorkUnit {
        int tid;
        int start_pos;      // BAM region start (0 for whole chrom)
        int end_pos;        // BAM region end (INT_MAX for whole chrom)
        size_t snp_start;   // Index into ChromSNPs vector
        size_t snp_end;     // Index into ChromSNPs vector (exclusive)
        bool has_snps;
        uint64_t est_reads; // Estimated reads for sorting
    };
    
    // Build work units - chunk based on SNP density OR read density
    vector<WorkUnit> work_units;
    long total_snps = 0;
    int chroms_with_snps = 0;
    int chroms_chunked_by_snp = 0;
    int chroms_chunked_by_reads = 0;
    
    // Thresholds for chunking
    const size_t CHUNK_SNP_THRESHOLD = 100000;       // Chunk if >100k SNPs
    const uint64_t CHUNK_READ_THRESHOLD = 10000000;  // Chunk if >10M reads
    
    for (int tid = 0; tid < n_chroms; tid++){
        auto it = snpdat_all.find(tid);
        int64_t chrom_len = chrom_lengths[tid];
        uint64_t chrom_reads = chrom_read_counts[tid];
        
        if (it == snpdat_all.end() || it->second.empty()){
            // No SNPs - check if high read density warrants chunking
            if (chrom_reads > CHUNK_READ_THRESHOLD){
                // High read count - split into chunks by position
                size_t n_chunks = (chrom_reads + CHUNK_READ_THRESHOLD - 1) / CHUNK_READ_THRESHOLD;
                // Cap at reasonable number of chunks
                n_chunks = std::min(n_chunks, (size_t)20);
                int64_t chunk_size = (chrom_len + n_chunks - 1) / n_chunks;
                uint64_t reads_per_chunk = chrom_reads / n_chunks;
                
                for (size_t c = 0; c < n_chunks; c++){
                    int start_pos = c * chunk_size;
                    int end_pos = (c == n_chunks - 1) ? INT_MAX : (int)((c + 1) * chunk_size);
                    work_units.push_back({tid, start_pos, end_pos, 0, 0, false, reads_per_chunk});
                }
                chroms_chunked_by_reads++;
            }
            else{
                // Low read count - single unit
                work_units.push_back({tid, 0, INT_MAX, 0, 0, false, chrom_reads});
            }
        }
        else{
            ChromSNPs& chrom_snps = it->second;
            size_t n_snps = chrom_snps.snps.size();
            total_snps += n_snps;
            chroms_with_snps++;
            
            // Decide chunking strategy: SNP density or read density, whichever is higher
            size_t chunks_by_snp = (n_snps > CHUNK_SNP_THRESHOLD) ? 
                (n_snps + CHUNK_SNP_THRESHOLD - 1) / CHUNK_SNP_THRESHOLD : 1;
            size_t chunks_by_reads = (chrom_reads > CHUNK_READ_THRESHOLD) ?
                (chrom_reads + CHUNK_READ_THRESHOLD - 1) / CHUNK_READ_THRESHOLD : 1;
            
            // Cap read-based chunks
            chunks_by_reads = std::min(chunks_by_reads, (size_t)20);
            
            if (chunks_by_snp >= chunks_by_reads && chunks_by_snp > 1){
                // Chunk by SNP count
                size_t snps_per_chunk = (n_snps + chunks_by_snp - 1) / chunks_by_snp;
                uint64_t reads_per_chunk = chrom_reads / chunks_by_snp;
                
                for (size_t c = 0; c < chunks_by_snp; c++){
                    size_t snp_start = c * snps_per_chunk;
                    size_t snp_end = std::min(snp_start + snps_per_chunk, n_snps);
                    
                    int start_pos = (snp_start == 0) ? 0 : chrom_snps.snps[snp_start].pos;
                    int end_pos = (snp_end >= n_snps) ? INT_MAX : chrom_snps.snps[snp_end - 1].pos + 1000;
                    
                    work_units.push_back({tid, start_pos, end_pos, snp_start, snp_end, true, reads_per_chunk});
                }
                chroms_chunked_by_snp++;
            }
            else if (chunks_by_reads > 1){
                // Chunk by read density - split position space
                int64_t chunk_size = (chrom_len + chunks_by_reads - 1) / chunks_by_reads;
                uint64_t reads_per_chunk = chrom_reads / chunks_by_reads;
                
                size_t snp_idx = 0;
                for (size_t c = 0; c < chunks_by_reads; c++){
                    int start_pos = c * chunk_size;
                    int end_pos = (c == chunks_by_reads - 1) ? INT_MAX : (int)((c + 1) * chunk_size);
                    
                    // Find SNPs in this position range
                    size_t snp_start = snp_idx;
                    while (snp_idx < n_snps && chrom_snps.snps[snp_idx].pos < (c + 1) * chunk_size){
                        snp_idx++;
                    }
                    size_t snp_end = snp_idx;
                    
                    work_units.push_back({tid, start_pos, end_pos, snp_start, snp_end, true, reads_per_chunk});
                }
                chroms_chunked_by_reads++;
            }
            else{
                // No chunking needed - single unit
                work_units.push_back({tid, 0, INT_MAX, 0, n_snps, true, chrom_reads});
            }
        }
    }
    
    // Sort work units by estimated reads (highest first for best load balancing)
    std::sort(work_units.begin(), work_units.end(),
              [](const WorkUnit& a, const WorkUnit& b){
                  return a.est_reads > b.est_reads;
              });
    
    fprintf(stderr, "Processing %d chromosomes (%d with SNPs, %ld total SNPs) using %d threads...\n",
        n_chroms, chroms_with_snps, total_snps, n_threads);
    fprintf(stderr, "  Split into %lu work units (%d by SNP density, %d by read density)\n", 
        work_units.size(), chroms_chunked_by_snp, chroms_chunked_by_reads);
    if (work_units.size() > 0){
        const WorkUnit& largest = work_units[0];
        fprintf(stderr, "  Largest work unit: %lu SNPs, ~%luM reads\n", 
            largest.snp_end - largest.snp_start, largest.est_reads / 1000000);
    }
    
    // Progress tracking
    atomic<long> snps_processed(0);
    atomic<int> units_done(0);
    atomic<long> reads_processed(0);
    
    // Per-thread cell counts for memory-efficient accumulation
    // With fixed-point int64 arithmetic, merge order doesn't matter (integer addition is associative)
    omp_set_num_threads(n_threads);
    vector<robin_hood::unordered_map<unsigned long, CellCounts>> thread_counts(n_threads);
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        auto& local_counts = thread_counts[thread_id];
        
        // Each thread gets its own BAM reader
        htsFile* bam_fp = hts_open(bamfile.c_str(), "r");
        if (!bam_fp){
            fprintf(stderr, "ERROR: Thread %d could not open BAM file\n", omp_get_thread_num());
        }
        else{
            hts_set_threads(bam_fp, htslib_threads);
            
            bam_hdr_t* header = sam_hdr_read(bam_fp);
            hts_idx_t* idx = sam_index_load(bam_fp, bamfile.c_str());
            bam1_t* record = bam_init1();
            
            if (!idx){
                fprintf(stderr, "ERROR: Could not load BAM index\n");
            }
            else{
                // Process work units with dynamic scheduling
                #pragma omp for schedule(dynamic, 1)
                for (size_t i = 0; i < work_units.size(); i++){
                    WorkUnit& wu = work_units[i];
                    int tid = wu.tid;
                    
                    // Query BAM for this region
                    hts_itr_t* iter = sam_itr_queryi(idx, tid, wu.start_pos, wu.end_pos);
                    if (!iter) continue;
                    
                    long local_snps = 0;
                    long local_reads = 0;
                    long local_all_reads = 0;
                    
                    if (!wu.has_snps){
                        // No SNPs on this chromosome - just count reads
                        while (sam_itr_next(bam_fp, iter, record) >= 0){
                            if (record->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)){
                                continue;
                            }
                            local_all_reads++;
                        }
                        hts_itr_destroy(iter);
                        reads_processed += local_all_reads;
                        
                        int done = ++units_done;
                        if (done % 100 == 0 || done == (int)work_units.size()){
                            fprintf(stderr, "\rProgress: %d/%lu units, %ld/%ld SNPs, %ld reads",
                                done, work_units.size(), snps_processed.load(), total_snps,
                                reads_processed.load());
                        }
                        continue;
                    }
                    
                    // Has SNPs - do full processing for this chunk
                    auto snp_it = snpdat_all.find(tid);
                    ChromSNPs& chrom_snps = snp_it->second;
                    
                    // Get iterators for just our chunk of SNPs
                    auto snp_iter = chrom_snps.snps.begin() + wu.snp_start;
                    auto snp_chunk_end = chrom_snps.snps.begin() + wu.snp_end;
                    
                    // Get chromosome name for progress reporting
                    const char* chrom_name = header->target_name[tid];
                    long chunk_snp_count = wu.snp_end - wu.snp_start;
                    
                    while (sam_itr_next(bam_fp, iter, record) >= 0){
                        // Skip filtered reads - match V1 behavior exactly
                        // V1 filters: unmapped, secondary, qcfail, dup
                        if (record->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)){
                            continue;
                        }
                        
                        // Count ALL reads that pass flag filter (before CB check) - matches V1
                        local_all_reads++;
                        
                        // Progress within large chunks (every 5M reads)
                        if (local_all_reads % 5000000 == 0){
                            fprintf(stderr, "\r  [%s:%d-%d] %ldM reads, %ld/%ld SNPs...          ",
                                chrom_name, wu.start_pos, wu.end_pos, 
                                local_all_reads / 1000000, local_snps, chunk_snp_count);
                        }
                        
                        int read_start = record->core.pos;
                        int read_end = bam_endpos(record);
                        
                        // Advance SNP iterator past SNPs before this read (within our chunk)
                        while (snp_iter != snp_chunk_end && snp_iter->pos < read_start){
                            ++snp_iter;
                            ++local_snps;
                        }
                        
                        // Extract cell barcode
                        uint8_t* cb_tag = bam_aux_get(record, "CB");
                        if (!cb_tag) continue;
                        
                        const char* cb_str = bam_aux2Z(cb_tag);
                        bc cb_bits;
                        str2bc(cb_str, cb_bits);
                        unsigned long bc_key = cb_bits.to_ulong();
                        
                        // Skip if not in whitelist
                        if (has_bc_list && valid_barcodes.find(bc_key) == valid_barcodes.end()){
                            continue;
                        }
                        
                        local_reads++;
                        
                        // Get mapping quality probability and scale to fixed-point
                        float prob_correct = 1.0f - powf(10.0f, -(float)record->core.qual / 10.0f);
                        int64_t prob_scaled = (int64_t)(prob_correct * FIXED_POINT_SCALE);
                        
                        // Process all SNPs overlapping this read (within our chunk)
                        for (auto snp_check = snp_iter; 
                             snp_check != snp_chunk_end && snp_check->pos <= read_end; 
                             ++snp_check){
                            
                            char allele = get_base_at_pos(record, snp_check->pos);
                            if (allele == 'N' || allele == '-') continue;
                            
                            int64_t ref_add = 0, alt_add = 0;
                            
                            if (allele == snp_check->data.ref){
                                ref_add = prob_scaled;
                                g_base_valid++;  // Count ref matches
                            }
                            else if (allele == snp_check->data.alt){
                                alt_add = prob_scaled;
                                g_base_valid++;  // Count alt matches
                            }
                            else {
                                g_base_other++;  // non-ref, non-alt base
                            }
                            
                            if (ref_add > 0 || alt_add > 0){
                                // Accumulate to per-thread storage (no locks needed)
                                auto it = local_counts.find(bc_key);
                                if (it == local_counts.end()){
                                    local_counts.emplace(bc_key, CellCounts(n_samples));
                                    it = local_counts.find(bc_key);
                                }
                                CellCounts& cc = it->second;
                                
                                // Update for all individuals (fixed-point arithmetic)
                                for (int indv = 0; indv < n_samples; ++indv){
                                    if (!snp_check->data.haps_covered.test(indv)) continue;
                                    
                                    int nalt = 0;
                                    if (snp_check->data.haps1.test(indv)) nalt++;
                                    if (snp_check->data.haps2.test(indv)) nalt++;
                                    
                                    // Store total for this individual
                                    cc.add_total(indv, nalt, ref_add, alt_add);
                                    
                                    // Pairwise counts
                                    for (int indv2 = indv + 1; indv2 < n_samples; ++indv2){
                                        if (!snp_check->data.haps_covered.test(indv2)) continue;
                                        
                                        int nalt2 = 0;
                                        if (snp_check->data.haps1.test(indv2)) nalt2++;
                                        if (snp_check->data.haps2.test(indv2)) nalt2++;
                                        
                                        cc.add(indv, nalt, indv2, nalt2, ref_add, alt_add);
                                    }
                                }
                            }
                        }
                    }
                    
                    // Count remaining SNPs in this chunk
                    while (snp_iter != snp_chunk_end){
                        ++snp_iter;
                        ++local_snps;
                    }
                    
                    snps_processed += local_snps;
                    reads_processed += local_all_reads;
                    int done = ++units_done;
                    
                    if (done % 10 == 0 || done == (int)work_units.size()){
                        fprintf(stderr, "\rProgress: %d/%lu units, %ld/%ld SNPs, %ld reads          ",
                            done, work_units.size(), snps_processed.load(), total_snps,
                            reads_processed.load());
                    }
                    
                    hts_itr_destroy(iter);
                }
                
                hts_idx_destroy(idx);
            }
            
            bam_destroy1(record);
            bam_hdr_destroy(header);
            hts_close(bam_fp);
        }
    }
    
    // Merge per-thread counts (order doesn't matter with fixed-point int64 arithmetic)
    fprintf(stderr, "\nMerging per-thread counts (fixed-point)...\n");
    
    // Efficient merge: iterate each thread's map once and merge directly
    // This is O(total_entries) instead of O(cells * threads) lookups
    for (int t = 0; t < n_threads; t++){
        for (auto& kv : thread_counts[t]){
            unsigned long bc_key = kv.first;
            auto it = cell_counts.find(bc_key);
            if (it == cell_counts.end()){
                cell_counts.emplace(std::piecewise_construct,
                    std::forward_as_tuple(bc_key),
                    std::forward_as_tuple(n_samples));
                it = cell_counts.find(bc_key);
            }
            it->second.counts.merge(kv.second);
        }
        // Free memory incrementally
        thread_counts[t].clear();
        
        if ((t + 1) % 10 == 0 || t == n_threads - 1){
            fprintf(stderr, "  Merged thread %d/%d\n", t + 1, n_threads);
        }
    }
    
    // Final cleanup
    thread_counts.clear();
    thread_counts.shrink_to_fit();
    
    fprintf(stderr, "Completed: %d chromosomes (%lu work units), %ld SNPs, %ld reads, %lu cells\n",
        n_chroms, work_units.size(), snps_processed.load(), reads_processed.load(), 
        cell_counts.size());
    
    // Debug: print base extraction stats
    fprintf(stderr, "DEBUG base extraction stats:\n");
    fprintf(stderr, "  reads processed: %ld\n", reads_processed.load());
    fprintf(stderr, "  queries: %ld\n", g_base_queries.load());
    fprintf(stderr, "  valid (ref/alt): %ld\n", g_base_valid.load());
    fprintf(stderr, "  N:       %ld\n", g_base_N.load());
    fprintf(stderr, "  dash:    %ld\n", g_base_dash.load());
    fprintf(stderr, "  other:   %ld\n", g_base_other.load());
    fprintf(stderr, "  queries per read: %.2f\n", (double)g_base_queries.load() / reads_processed.load());
}

void count_alleles_single_threaded(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    robin_hood::unordered_map<unsigned long, CellCounts>& cell_counts,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    map<pair<int, int>, map<int, float> >& conditional_match_fracs,
    map<pair<int, int>, map<int, float> >& conditional_match_tots,
    bool compute_conditional){
    
    bool has_bc_list = !valid_barcodes.empty();
    
    string bamfile_copy = bamfile;  // bam_reader needs non-const
    bam_reader reader;
    reader.set_file(bamfile_copy);
    reader.set_cb();
    
    long total_snps = 0;
    for (auto& kv : snpdat_all){
        total_snps += kv.second.size();
    }
    
    fprintf(stderr, "Processing %ld SNPs single-threaded...\n", total_snps);
    
    int curtid = -1;
    vector<SNPData>::iterator cursnp;
    vector<SNPData>::iterator cursnp_end;
    ChromSNPs* active_chrom = nullptr;
    
    long snps_processed = 0;
    long reads_processed = 0;
    int progress = 10000;
    
    while (reader.next()){
        if (reader.unmapped() || reader.secondary() || reader.qcfail() || reader.dup()){
            continue;
        }
        
        // Count ALL reads after flag filter (like V1)
        reads_processed++;
        
        if (curtid != reader.tid()){
            // New chromosome
            curtid = reader.tid();
            
            auto it = snpdat_all.find(curtid);
            if (it != snpdat_all.end() && !it->second.empty()){
                active_chrom = &(it->second);
                cursnp = active_chrom->snps.begin();
                cursnp_end = active_chrom->snps.end();
                
                if (compute_conditional){
                    get_conditional_match_fracs_chrom_optimized(*active_chrom,
                        conditional_match_fracs, conditional_match_tots, n_samples);
                }
            }
            else{
                active_chrom = nullptr;
            }
        }
        
        if (!active_chrom) continue;
        
        // Advance SNP iterator BEFORE CB check (match V1 behavior)
        while (cursnp != cursnp_end && cursnp->pos < reader.reference_start){
            ++cursnp;
            ++snps_processed;
        }
        
        if (!reader.has_cb_z) continue;
        
        bc cb_bits;
        str2bc(reader.cb_z, cb_bits);
        unsigned long bc_key = cb_bits.to_ulong();
        
        if (has_bc_list && valid_barcodes.find(bc_key) == valid_barcodes.end()){
            continue;
        }
        
        // Get mapping quality probability and scale to fixed-point
        float prob_correct = 1.0f - powf(10.0f, -(float)reader.mapq / 10.0f);
        int64_t prob_scaled = (int64_t)(prob_correct * FIXED_POINT_SCALE);
        
        // Process overlapping SNPs
        for (auto snp_check = cursnp;
             snp_check != cursnp_end && snp_check->pos <= reader.reference_end;
             ++snp_check){
            
            char allele = reader.get_base_at(snp_check->pos + 1);
            if (allele == 'N' || allele == '-') continue;
            
            int64_t ref_add = 0, alt_add = 0;
            
            if (allele == snp_check->data.ref){
                ref_add = prob_scaled;
            }
            else if (allele == snp_check->data.alt){
                alt_add = prob_scaled;
            }
            
            if (ref_add > 0 || alt_add > 0){
                // Get or create cell counts
                auto it = cell_counts.find(bc_key);
                if (it == cell_counts.end()){
                    cell_counts.emplace(bc_key, CellCounts(n_samples));
                    it = cell_counts.find(bc_key);
                }
                
                // Update counts (fixed-point arithmetic)
                for (int indv = 0; indv < n_samples; ++indv){
                    if (!snp_check->data.haps_covered.test(indv)) continue;
                    
                    int nalt = 0;
                    if (snp_check->data.haps1.test(indv)) nalt++;
                    if (snp_check->data.haps2.test(indv)) nalt++;
                    
                    it->second.add_total(indv, nalt, ref_add, alt_add);
                    
                    for (int indv2 = indv + 1; indv2 < n_samples; ++indv2){
                        if (!snp_check->data.haps_covered.test(indv2)) continue;
                        
                        int nalt2 = 0;
                        if (snp_check->data.haps1.test(indv2)) nalt2++;
                        if (snp_check->data.haps2.test(indv2)) nalt2++;
                        
                        it->second.add(indv, nalt, indv2, nalt2, ref_add, alt_add);
                    }
                }
            }
        }
        
        if (reads_processed % progress == 0){
            fprintf(stderr, "Processed %ld reads, %ld SNPs\r", reads_processed, snps_processed);
        }
    }
    
    fprintf(stderr, "\nCompleted: %ld SNPs, %ld reads, %lu cells\n",
        snps_processed, reads_processed, cell_counts.size());
}

void finalize_parallel_counts(
    robin_hood::unordered_map<unsigned long, AlignedCellCounts>& parallel_counts,
    robin_hood::unordered_map<unsigned long, CellCounts>& final_counts){
    
    final_counts.clear();
    
    for (auto& kv : parallel_counts){
        final_counts.emplace(kv.first, std::move(kv.second.counts));
    }
}

// ============================================================================
// SHARED MEMORY VCF FUNCTIONS
// ============================================================================

bool create_shared_vcf(
    const string& vcf_file,
    const string& shm_name,
    set<string>& chroms_to_include,
    map<string, int>& seq2tid,
    int min_vq){
    
    // Load VCF data
    robin_hood::unordered_map<int, ChromSNPs> snpdat;
    string vcf_file_copy = vcf_file;  // read_vcf_chroms_optimized needs non-const
    int nvar = read_vcf_chroms_optimized(vcf_file_copy, chroms_to_include, seq2tid, snpdat, min_vq);
    
    if (nvar == 0){
        fprintf(stderr, "ERROR: No variants loaded from VCF\n");
        return false;
    }
    
    // Read sample names from VCF
    vector<string> samples;
    read_vcf_samples(vcf_file_copy, samples);
    fprintf(stderr, "Found %lu samples in VCF\n", samples.size());
    
    // Calculate total size needed
    size_t total_size = sizeof(SharedVCFHeader);
    int n_chroms = 0;
    
    for (auto& kv : snpdat){
        total_size += kv.second.snps.size() * sizeof(SNPData);
        n_chroms++;
    }
    
    fprintf(stderr, "Creating shared memory segment: %s (%.2f GB, %d SNPs, %d chroms)\n",
        shm_name.c_str(), (double)total_size / (1024.0 * 1024.0 * 1024.0), nvar, n_chroms);
    
    // Create shared memory
    int shm_fd = shm_open(shm_name.c_str(), O_CREAT | O_RDWR, 0666);
    if (shm_fd < 0){
        perror("shm_open");
        return false;
    }
    
    if (ftruncate(shm_fd, total_size) < 0){
        perror("ftruncate");
        close(shm_fd);
        shm_unlink(shm_name.c_str());
        return false;
    }
    
    void* ptr = mmap(NULL, total_size, PROT_READ | PROT_WRITE, MAP_SHARED, shm_fd, 0);
    if (ptr == MAP_FAILED){
        perror("mmap");
        close(shm_fd);
        shm_unlink(shm_name.c_str());
        return false;
    }
    
    // Serialize VCF data
    SharedVCFHeader* header = (SharedVCFHeader*)ptr;
    header->total_size = total_size;
    header->n_chromosomes = n_chroms;
    header->n_snps_total = nvar;
    
    // Store sample names
    header->n_samples = samples.size();
    for (size_t i = 0; i < samples.size() && i < 512; i++){
        strncpy(header->sample_names[i], samples[i].c_str(), 63);
        header->sample_names[i][63] = '\0';
    }
    fprintf(stderr, "Stored %d sample names in shared memory\n", header->n_samples);
    
    char* data_ptr = (char*)ptr + sizeof(SharedVCFHeader);
    size_t offset = sizeof(SharedVCFHeader);
    int chrom_idx = 0;
    
    for (auto& kv : snpdat){
        header->chrom_offsets[chrom_idx] = offset;
        header->chrom_snp_counts[chrom_idx] = kv.second.snps.size();
        header->chrom_tids[chrom_idx] = kv.first;
        
        // Copy SNP data
        size_t copy_size = kv.second.snps.size() * sizeof(SNPData);
        memcpy(data_ptr, kv.second.snps.data(), copy_size);
        data_ptr += copy_size;
        offset += copy_size;
        chrom_idx++;
    }
    
    // Sync to ensure data is written
    msync(ptr, total_size, MS_SYNC);
    
    fprintf(stderr, "Shared memory created successfully\n");
    
    // Keep mapping but close fd (other processes will open separately)
    close(shm_fd);
    
    return true;
}

bool attach_shared_vcf(
    const string& shm_name,
    robin_hood::unordered_map<int, ChromSNPs>& snpdat_all,
    vector<string>& samples){
    
    int shm_fd = shm_open(shm_name.c_str(), O_RDONLY, 0);
    if (shm_fd < 0){
        perror("shm_open");
        return false;
    }
    
    struct stat sb;
    if (fstat(shm_fd, &sb) < 0){
        perror("fstat");
        close(shm_fd);
        return false;
    }
    
    void* ptr = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, shm_fd, 0);
    if (ptr == MAP_FAILED){
        perror("mmap");
        close(shm_fd);
        return false;
    }
    
    SharedVCFHeader* header = (SharedVCFHeader*)ptr;
    
    fprintf(stderr, "Attached to shared VCF: %s (%d SNPs, %d chromosomes)\n",
        shm_name.c_str(), header->n_snps_total, header->n_chromosomes);
    
    // Retrieve sample names
    samples.clear();
    for (int i = 0; i < header->n_samples && i < 512; i++){
        samples.push_back(string(header->sample_names[i]));
    }
    fprintf(stderr, "Loaded %lu sample names from shared VCF\n", samples.size());
    
    // Deserialize
    snpdat_all.clear();
    
    for (int i = 0; i < header->n_chromosomes; i++){
        int tid = header->chrom_tids[i];
        size_t n_snps = header->chrom_snp_counts[i];
        
        SNPData* snp_ptr = (SNPData*)((char*)ptr + header->chrom_offsets[i]);
        
        ChromSNPs& cs = snpdat_all[tid];
        cs.snps.reserve(n_snps);
        
        for (size_t j = 0; j < n_snps; j++){
            cs.snps.push_back(snp_ptr[j]);
        }
    }
    
    // Note: We keep the mapping active - caller should call detach when done
    close(shm_fd);
    
    return true;
}

void detach_shared_vcf(const string& shm_name){
    // The mmap is automatically unmapped when the process exits
    // This function is a placeholder for explicit cleanup if needed
}

void destroy_shared_vcf(const string& shm_name){
    if (shm_unlink(shm_name.c_str()) < 0){
        perror("shm_unlink");
    }
    else{
        fprintf(stderr, "Shared memory destroyed: %s\n", shm_name.c_str());
    }
}

// ============================================================================
// HET VCF FUNCTIONS (NEW - for ploidy detection)
// ============================================================================

int load_het_vcf(
    const string& het_vcf_file,
    const set<string>& chroms_to_include,
    map<string, int>& seq2tid,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    int min_vq){
    
    // This is essentially read_vcf_chroms_optimized but for het VCF
    // We can reuse that function since the format is the same
    string vcf_file_copy = het_vcf_file;
    set<string> chroms_copy = chroms_to_include;
    
    return read_vcf_chroms_optimized(vcf_file_copy, chroms_copy, seq2tid, het_snpdat, min_vq);
}

void count_het_alleles_parallel(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    robin_hood::unordered_map<unsigned long, CellCounts>& het_counts,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads){
    
    // This is essentially count_alleles_parallel but with a different output structure
    // We can use the same parallel counting infrastructure
    
    robin_hood::unordered_map<unsigned long, AlignedCellCounts> parallel_counts;
    
    // Use the same counting function - it works with any SNP set
    count_alleles_parallel(bamfile, het_snpdat, parallel_counts,
        valid_barcodes, n_samples, n_threads, htslib_threads);
    
    // Finalize into the output structure
    finalize_parallel_counts(parallel_counts, het_counts);
    
    fprintf(stderr, "Het allele counting complete: %lu cells\n", het_counts.size());
}

/**
 * Extended het allele counting that collects per-site data and/or Welford stats.
 * 
 * Unlike count_het_alleles_parallel which aggregates by genotype pair (losing per-site info),
 * this preserves per-site information needed for accurate het balance variance.
 */
void count_het_alleles_extended(
    const string& bamfile,
    robin_hood::unordered_map<int, ChromSNPs>& het_snpdat,
    robin_hood::unordered_map<unsigned long, CellHetData>& het_data,
    vector<pair<int, int>>& idx_to_site,
    const set<unsigned long>& valid_barcodes,
    int n_samples,
    int n_threads,
    int htslib_threads,
    HetBalanceMethod method) {
    
    bool collect_persite = (method == HetBalanceMethod::PERSITE);
    bool collect_welford = (method == HetBalanceMethod::WELFORD);
    bool has_bc_list = !valid_barcodes.empty();
    
    // Build global site index for persite method
    map<pair<int, int>, int32_t> site_to_idx;
    idx_to_site.clear();
    int32_t global_idx = 0;
    
    for (auto& kv : het_snpdat) {
        int tid = kv.first;
        for (const auto& snp : kv.second.snps) {
            site_to_idx[{tid, snp.pos}] = global_idx;
            idx_to_site.push_back({tid, snp.pos});
            global_idx++;
        }
    }
    
    long total_het_sites = global_idx;
    const char* method_name = collect_persite ? "per-site" : "Welford";
    fprintf(stderr, "Processing %ld het sites with %s method\n", total_het_sites, method_name);
    
    // Pre-allocate for known barcodes
    if (has_bc_list) {
        fprintf(stderr, "Pre-allocating het data for %lu cells...\n", valid_barcodes.size());
        for (unsigned long bc : valid_barcodes) {
            het_data.emplace(piecewise_construct,
                forward_as_tuple(bc),
                forward_as_tuple(n_samples));
        }
    }
    
    // Get chromosome info from BAM header
    htsFile* bam_tmp = hts_open(bamfile.c_str(), "r");
    if (!bam_tmp) {
        fprintf(stderr, "ERROR: Could not open BAM file\n");
        return;
    }
    bam_hdr_t* hdr_tmp = sam_hdr_read(bam_tmp);
    int n_chroms = hdr_tmp->n_targets;
    bam_hdr_destroy(hdr_tmp);
    hts_close(bam_tmp);
    
    // Build work units by chromosome
    vector<int> chrom_work;
    for (int tid = 0; tid < n_chroms; tid++) {
        if (het_snpdat.find(tid) != het_snpdat.end() && !het_snpdat[tid].empty()) {
            chrom_work.push_back(tid);
        }
    }
    
    fprintf(stderr, "Processing %lu chromosomes with %d threads...\n", chrom_work.size(), n_threads);
    
    mutex output_mutex;
    atomic<int> chroms_done(0);
    atomic<long> reads_processed(0);
    atomic<long> sites_hit(0);
    
    #pragma omp parallel num_threads(n_threads)
    {
        // Thread-local het data
        robin_hood::unordered_map<unsigned long, CellHetData> local_het_data;
        
        if (has_bc_list) {
            for (unsigned long bc : valid_barcodes) {
                local_het_data.emplace(piecewise_construct,
                    forward_as_tuple(bc),
                    forward_as_tuple(n_samples));
            }
        }
        
        // Open BAM for this thread
        htsFile* bam_fp = hts_open(bamfile.c_str(), "r");
        if (!bam_fp) {
            fprintf(stderr, "ERROR: Thread could not open BAM\n");
        }
        else {
            if (htslib_threads > 1) {
                hts_set_threads(bam_fp, htslib_threads);
            }
            
            bam_hdr_t* header = sam_hdr_read(bam_fp);
            hts_idx_t* idx = sam_index_load(bam_fp, bamfile.c_str());
            bam1_t* record = bam_init1();
            
            if (idx) {
                #pragma omp for schedule(dynamic)
                for (size_t wi = 0; wi < chrom_work.size(); wi++) {
                    int tid = chrom_work[wi];
                    ChromSNPs& chrom_snps = het_snpdat[tid];
                    
                    long local_reads = 0;
                    long local_sites = 0;
                    
                    hts_itr_t* iter = sam_itr_queryi(idx, tid, 0, INT_MAX);
                    if (!iter) continue;
                    
                    auto snp_iter = chrom_snps.snps.begin();
                    auto snp_end = chrom_snps.snps.end();
                    
                    while (sam_itr_next(bam_fp, iter, record) >= 0) {
                        local_reads++;
                        
                        if (record->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
                            continue;
                        }
                        
                        uint8_t* cb_ptr = bam_aux_get(record, "CB");
                        if (!cb_ptr) continue;
                        
                        const char* cb_str = bam_aux2Z(cb_ptr);
                        if (!cb_str) continue;
                        
                        unsigned long bc_key = bc_ul((char*)cb_str);
                        
                        if (has_bc_list && valid_barcodes.find(bc_key) == valid_barcodes.end()) {
                            continue;
                        }
                        
                        int read_start = record->core.pos;
                        int read_end = bam_endpos(record);
                        
                        // Advance SNP iterator
                        while (snp_iter != snp_end && snp_iter->pos < read_start) {
                            ++snp_iter;
                        }
                        
                        float prob_correct = 1.0f - powf(10.0f, -(float)record->core.qual / 10.0f);
                        
                        // Process overlapping SNPs
                        for (auto snp_check = snp_iter;
                             snp_check != snp_end && snp_check->pos <= read_end;
                             ++snp_check) {
                            
                            char allele = get_base_at_pos(record, snp_check->pos);
                            if (allele == 'N' || allele == '-') continue;
                            
                            float ref_add = 0.0f, alt_add = 0.0f;
                            
                            if (allele == snp_check->data.ref) {
                                ref_add = prob_correct;
                            } else if (allele == snp_check->data.alt) {
                                alt_add = prob_correct;
                            } else {
                                continue;
                            }
                            
                            local_sites++;
                            
                            // Get or create cell's het data
                            auto it = local_het_data.find(bc_key);
                            if (it == local_het_data.end()) {
                                if (has_bc_list) continue;
                                local_het_data.emplace(bc_key, CellHetData(n_samples));
                                it = local_het_data.find(bc_key);
                            }
                            CellHetData& cell_data = it->second;
                            
                            // Per-site method: store site counts
                            if (collect_persite) {
                                int32_t site_idx = site_to_idx[{tid, snp_check->pos}];
                                cell_data.persite_data.add_site(site_idx, ref_add, alt_add);
                            }
                            
                            // Welford method: update running stats for het individuals
                            if (collect_welford) {
                                float depth = ref_add + alt_add;
                                if (depth > 0) {
                                    float alt_frac = alt_add / depth;
                                    
                                    for (int indiv = 0; indiv < n_samples; ++indiv) {
                                        if (snp_check->data.is_het(indiv)) {
                                            cell_data.welford_stats.add(indiv, alt_frac, depth);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                    hts_itr_destroy(iter);
                    reads_processed += local_reads;
                    sites_hit += local_sites;
                    
                    int done = ++chroms_done;
                    if (done % 5 == 0 || done == (int)chrom_work.size()) {
                        fprintf(stderr, "\rHet counting: %d/%lu chroms, %ld reads, %ld site hits    ",
                                done, chrom_work.size(), reads_processed.load(), sites_hit.load());
                    }
                }
                
                hts_idx_destroy(idx);
            }
            
            bam_destroy1(record);
            bam_hdr_destroy(header);
            hts_close(bam_fp);
        }
        
        // Merge thread-local data
        #pragma omp critical
        {
            for (auto& kv : local_het_data) {
                auto it = het_data.find(kv.first);
                if (it != het_data.end()) {
                    it->second.merge(kv.second);
                } else {
                    het_data.emplace(kv.first, std::move(kv.second));
                }
            }
        }
    }
    
    fprintf(stderr, "\nHet counting complete: %lu cells, %ld reads, %ld site hits\n",
            het_data.size(), reads_processed.load(), sites_hit.load());
    
    // Report stats
    if (collect_persite) {
        long total_sites_stored = 0;
        long max_sites = 0;
        for (const auto& kv : het_data) {
            long n = kv.second.persite_data.size();
            total_sites_stored += n;
            if (n > max_sites) max_sites = n;
        }
        fprintf(stderr, "Per-site: %ld total entries, max %ld/cell, avg %.1f/cell\n",
                total_sites_stored, max_sites, 
                het_data.size() > 0 ? (double)total_sites_stored / het_data.size() : 0.0);
    }
    
    if (collect_welford) {
        long max_sites = 0;
        for (const auto& kv : het_data) {
            for (int i = 0; i < n_samples; ++i) {
                long n = kv.second.welford_stats.get(i).n;
                if (n > max_sites) max_sites = n;
            }
        }
        fprintf(stderr, "Max het sites per individual per cell: %ld\n", max_sites);
    }
}
