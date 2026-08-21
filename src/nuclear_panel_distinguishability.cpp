#include <htslib/vcf.h>
#include <getopt.h>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

struct PairStats {
    uint64_t joint_callable = 0;
    uint64_t discordant = 0;
    long double expected_information = 0.0L;
};

static void usage(FILE* out) {
    std::fprintf(out,
        "nuclear_panel_distinguishability --vcf PANEL.bcf --output FILE.tsv\n"
        "Compute pairwise donor distinguishability on a biallelic diploid panel.\n"
        "expected_information is sum((alt_dosage1-alt_dosage2)^2) with dosage in [0,1].\n");
}

int main(int argc, char** argv) {
    std::string vcf_path, output_path;
    static struct option opts[] = {
        {"vcf", required_argument, 0, 'v'},
        {"output", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "v:o:h", opts, nullptr)) != -1) {
        if (c == 'v') vcf_path = optarg;
        else if (c == 'o') output_path = optarg;
        else { usage(c == 'h' ? stdout : stderr); return c == 'h' ? 0 : 2; }
    }
    if (vcf_path.empty() || output_path.empty()) { usage(stderr); return 2; }

    htsFile* fp = bcf_open(vcf_path.c_str(), "r");
    if (!fp) { std::fprintf(stderr, "ERROR: cannot open %s\n", vcf_path.c_str()); return 1; }
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    if (!hdr) { std::fprintf(stderr, "ERROR: cannot read header from %s\n", vcf_path.c_str()); bcf_close(fp); return 1; }
    const int ns = bcf_hdr_nsamples(hdr);
    if (ns < 2 || ns > 32) {
        std::fprintf(stderr, "ERROR: expected 2-32 donor samples, found %d\n", ns);
        bcf_hdr_destroy(hdr); bcf_close(fp); return 1;
    }

    // Two bits per donor: 0=0/0, 1=0/1, 2=1/1, 3=missing/non-diploid/non-biallelic.
    std::unordered_map<uint64_t, uint64_t> pattern_counts;
    pattern_counts.reserve(4096);
    bcf1_t* rec = bcf_init();
    int32_t* gt = nullptr;
    int ngt = 0;
    uint64_t accepted_sites = 0;
    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR);
        if (rec->n_allele != 2) continue;
        const int nret = bcf_get_genotypes(hdr, rec, &gt, &ngt);
        if (nret <= 0 || nret % ns != 0) continue;
        const int stride = nret / ns;
        uint64_t pattern = 0;
        for (int s = 0; s < ns; ++s) {
            unsigned code = 3;
            if (stride >= 2) {
                const int32_t g0 = gt[s * stride];
                const int32_t g1 = gt[s * stride + 1];
                if (!bcf_gt_is_missing(g0) && !bcf_gt_is_missing(g1) &&
                    g0 != bcf_int32_vector_end && g1 != bcf_int32_vector_end) {
                    const int a0 = bcf_gt_allele(g0);
                    const int a1 = bcf_gt_allele(g1);
                    if ((a0 == 0 || a0 == 1) && (a1 == 0 || a1 == 1)) {
                        code = static_cast<unsigned>(a0 + a1);
                    }
                }
            }
            pattern |= (static_cast<uint64_t>(code) << (2 * s));
        }
        ++pattern_counts[pattern];
        ++accepted_sites;
    }
    if (gt) std::free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    if (accepted_sites == 0) {
        std::fprintf(stderr, "ERROR: no biallelic sites with readable GT data in %s\n", vcf_path.c_str());
        return 1;
    }

    const size_t npairs = static_cast<size_t>(ns) * (ns - 1) / 2;
    std::vector<PairStats> stats(npairs);
    auto pair_index = [ns](int i, int j) -> size_t {
        // i < j, compact upper triangle.
        return static_cast<size_t>(i) * (2 * ns - i - 1) / 2 + static_cast<size_t>(j - i - 1);
    };
    for (const auto& kv : pattern_counts) {
        const uint64_t pattern = kv.first;
        const uint64_t n = kv.second;
        unsigned code[32];
        for (int s = 0; s < ns; ++s) code[s] = static_cast<unsigned>((pattern >> (2 * s)) & 3ULL);
        for (int i = 0; i < ns; ++i) {
            if (code[i] == 3) continue;
            for (int j = i + 1; j < ns; ++j) {
                if (code[j] == 3) continue;
                PairStats& st = stats[pair_index(i, j)];
                st.joint_callable += n;
                if (code[i] != code[j]) st.discordant += n;
                const long double d = (static_cast<long double>(code[i]) - static_cast<long double>(code[j])) / 2.0L;
                st.expected_information += static_cast<long double>(n) * d * d;
            }
        }
    }

    // Reopen only the header to recover sample names after aggregation.
    fp = bcf_open(vcf_path.c_str(), "r");
    hdr = fp ? bcf_hdr_read(fp) : nullptr;
    if (!hdr) { if (fp) bcf_close(fp); std::fprintf(stderr, "ERROR: cannot reopen VCF header\n"); return 1; }
    std::ofstream out(output_path.c_str());
    if (!out) { std::fprintf(stderr, "ERROR: cannot open %s\n", output_path.c_str()); bcf_hdr_destroy(hdr); bcf_close(fp); return 1; }
    out << "donor1\tdonor2\tjoint_callable_sites\tdiscordant_genotype_sites\texpected_information\tconfusability_status\n";
    for (int i = 0; i < ns; ++i) {
        for (int j = i + 1; j < ns; ++j) {
            const PairStats& st = stats[pair_index(i, j)];
            const char* status = st.joint_callable == 0 ? "NO_JOINT_CALLABLE_SITES" :
                                 (st.discordant == 0 ? "INDISTINGUISHABLE_ON_PANEL" : "DISTINGUISHABLE_ON_PANEL");
            out << hdr->samples[i] << '\t' << hdr->samples[j] << '\t'
                << st.joint_callable << '\t' << st.discordant << '\t'
                << static_cast<double>(st.expected_information) << '\t' << status << '\n';
        }
    }
    out.close();
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    return out ? 0 : 1;
}
