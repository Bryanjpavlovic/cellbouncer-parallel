#include "mt_evidence.h"

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace mt_evidence {
namespace {

bool valid_fmt(int32_t v){ return v != bcf_int32_missing && v != bcf_int32_vector_end && v >= 0; }
int fmt_value(const int32_t* data, int total, int ns, int sample, int offset){
    if (!data || total <= 0 || ns <= 0 || total % ns != 0) return -1;
    const int stride = total / ns;
    if (offset < 0 || offset >= stride) return -1;
    int32_t v = data[sample * stride + offset];
    return valid_fmt(v) ? (int)v : -1;
}
int gt_state(const int32_t* gt, int total, int ns, int sample){
    if (!gt || total <= 0 || ns <= 0 || total % ns != 0) return -1;
    const int stride = total / ns;
    int state=-1; bool saw=false;
    for (int k=0;k<stride;++k){
        int32_t v=gt[sample*stride+k];
        if (v==bcf_int32_vector_end) break;
        if (bcf_gt_is_missing(v)) return -1;
        int a=bcf_gt_allele(v); if(a<0||a>1) return -3;
        if(!saw){state=a;saw=true;} else if(state!=a) return -2;
    }
    return saw?state:-1;
}

int base_at_ref_pos(const bam1_t* rec, int64_t target_pos, int* baseq){
    int64_t ref = rec->core.pos;
    int q = 0;
    uint32_t* cigar = bam_get_cigar(rec);
    for (uint32_t i=0;i<rec->core.n_cigar;++i){
        int op=bam_cigar_op(cigar[i]); int len=bam_cigar_oplen(cigar[i]);
        bool consumes_ref = bam_cigar_type(op)&2;
        bool consumes_q = bam_cigar_type(op)&1;
        if(consumes_ref && consumes_q){
            if(target_pos>=ref && target_pos<ref+len){
                int qi=q+(int)(target_pos-ref);
                if(qi<0||qi>=rec->core.l_qseq) return -1;
                uint8_t* seq=bam_get_seq(rec); uint8_t* qual=bam_get_qual(rec);
                if(baseq) *baseq = qual ? qual[qi] : 0;
                return bam_seqi(seq,qi);
            }
        } else if(consumes_ref && !consumes_q){
            if(target_pos>=ref && target_pos<ref+len) return -1;
        }
        if(consumes_ref) ref+=len;
        if(consumes_q) q+=len;
    }
    return -1;
}

char nt16_to_base(int code){
    static const char* lut="=ACMGRSVTWYHKDBN";
    return (code>=0&&code<16)?lut[code]:'N';
}

struct MoleculeKey {
    std::string barcode;
    int site=-1;
    std::string umi;
    bool operator==(const MoleculeKey& o) const { return site==o.site && barcode==o.barcode && umi==o.umi; }
};
struct MoleculeHash {
    size_t operator()(const MoleculeKey& k) const {
        size_t h=std::hash<std::string>()(k.barcode);
        h ^= std::hash<std::string>()(k.umi) + 0x9e3779b9 + (h<<6) + (h>>2);
        h ^= std::hash<int>()(k.site) + 0x9e3779b9 + (h<<6) + (h>>2);
        return h;
    }
};

} // namespace

std::string canonical_pair_key(const std::string& a,const std::string& b){ return a<=b ? a+"\x1f"+b : b+"\x1f"+a; }
std::vector<std::string> split_genotype(const std::string& g){
    std::vector<std::string> out; size_t start=0;
    while(start<=g.size()){ size_t p=g.find('+',start); std::string x=g.substr(start,p==std::string::npos?std::string::npos:p-start); if(!x.empty()) out.push_back(x); if(p==std::string::npos) break; start=p+1; }
    std::sort(out.begin(),out.end()); return out;
}

Panel load_panel(const std::string& path,const std::string& chrom,int min_depth,double homoplasmy_af,const std::unordered_set<int64_t>* allowed){
    htsFile* fp=bcf_open(path.c_str(),"r"); if(!fp) throw std::runtime_error("Could not open mt panel: "+path);
    bcf_hdr_t* hdr=bcf_hdr_read(fp); if(!hdr){hts_close(fp);throw std::runtime_error("Could not read mt panel header");}
    Panel panel; panel.chrom=chrom; int ns=bcf_hdr_nsamples(hdr);
    for(int i=0;i<ns;++i){ panel.samples.push_back(hdr->samples[i]); panel.sample_index[hdr->samples[i]]=i; }
    int rid=bcf_hdr_name2id(hdr,chrom.c_str()); if(rid<0){bcf_hdr_destroy(hdr);hts_close(fp);throw std::runtime_error("Mitochondrial contig missing from panel: "+chrom);}
    bcf1_t* rec=bcf_init(); int32_t *ro=nullptr,*ao=nullptr,*ad=nullptr,*gt=nullptr; int nro=0,nao=0,nad=0,ngt=0;
    while(bcf_read(fp,hdr,rec)>=0){
        if (rec->rid != rid) continue;
        if (allowed && !allowed->empty() && !allowed->count(rec->pos)) continue;
        bcf_unpack(rec,BCF_UN_STR); if(rec->n_allele!=2||strlen(rec->d.allele[0])!=1||strlen(rec->d.allele[1])!=1) continue;
        int ret_ro=bcf_get_format_int32(hdr,rec,"RO",&ro,&nro); int ret_ao=bcf_get_format_int32(hdr,rec,"AO",&ao,&nao); int ret_ad=bcf_get_format_int32(hdr,rec,"AD",&ad,&nad); int ret_gt=bcf_get_genotypes(hdr,rec,&gt,&ngt);
        Site site; site.rid=rid; site.pos=rec->pos; site.ref=rec->d.allele[0][0]; site.alt=rec->d.allele[1][0]; site.state.assign(ns,-1);
        for(int i=0;i<ns;++i){
            int rd=-1,xd=-1; if(ret_ro>0&&ret_ao>0){rd=fmt_value(ro,ret_ro,ns,i,0);xd=fmt_value(ao,ret_ao,ns,i,0);} if((rd<0||xd<0)&&ret_ad>0){rd=fmt_value(ad,ret_ad,ns,i,0);xd=fmt_value(ad,ret_ad,ns,i,1);} if(rd<0||xd<0||rd+xd<min_depth) continue;
            double af=(double)xd/(double)(rd+xd); int ds=af<=1.0-homoplasmy_af?0:(af>=homoplasmy_af?1:-1); if(ds<0) continue; int gs=gt_state(gt,ret_gt,ns,i); if(gs==-2||gs==-3||(gs>=0&&gs!=ds)) continue; site.state[i]=(int8_t)ds;
        }
        panel.sites.push_back(std::move(site));
    }
    free(ro);free(ao);free(ad);free(gt);bcf_destroy(rec);bcf_hdr_destroy(hdr);hts_close(fp);
    std::sort(panel.sites.begin(),panel.sites.end(),[](const Site&a,const Site&b){return a.pos<b.pos;});
    if (panel.sites.empty()) {
        throw std::runtime_error("No usable mt panel sites after depth/state filtering");
    }
    return panel;
}

CellCounts count_bam_once(const std::string& bam_path,const Panel& panel,const std::unordered_set<std::string>& targets,const std::string& barcode_tag,const std::string& umi_tag,int min_mapq,int min_baseq,int threads,CountStats& stats){
    samFile* fp=sam_open(bam_path.c_str(),"r"); if(!fp) throw std::runtime_error("Could not open BAM: "+bam_path); if(threads>1) hts_set_threads(fp,threads);
    bam_hdr_t* hdr=sam_hdr_read(fp); hts_idx_t* idx=sam_index_load(fp,bam_path.c_str()); if(!hdr||!idx){if(hdr)bam_hdr_destroy(hdr);if(idx)hts_idx_destroy(idx);sam_close(fp);throw std::runtime_error("BAM header/index unavailable: "+bam_path);} int tid=bam_name2id(hdr,panel.chrom.c_str()); if(tid<0){hts_idx_destroy(idx);bam_hdr_destroy(hdr);sam_close(fp);throw std::runtime_error("BAM lacks mitochondrial contig: "+panel.chrom);}
    hts_itr_t* itr=sam_itr_queryi(idx,tid,0,hdr->target_len[tid]); bam1_t* rec=bam_init1();
    std::unordered_map<MoleculeKey,uint8_t,MoleculeHash> mol;
    while(sam_itr_next(fp,itr,rec)>=0){
        ++stats.reads_seen; uint16_t f=rec->core.flag; if(f&(BAM_FUNMAP|BAM_FSECONDARY|BAM_FSUPPLEMENTARY|BAM_FQCFAIL|BAM_FDUP)) continue; if(rec->core.qual<min_mapq) continue;
        uint8_t* cbp=bam_aux_get(rec,barcode_tag.c_str()); uint8_t* ubp=bam_aux_get(rec,umi_tag.c_str()); if(!cbp||!ubp) continue; const char* cb=bam_aux2Z(cbp); const char* ub=bam_aux2Z(ubp); if(!cb||!ub||!targets.count(cb)) continue; ++stats.reads_accepted;
        int64_t rs=rec->core.pos,re=bam_endpos(rec); auto it=std::lower_bound(panel.sites.begin(),panel.sites.end(),rs,[](const Site&s,int64_t p){return s.pos<p;});
        for(;it!=panel.sites.end()&&it->pos<re;++it){ int bq=0; int code=base_at_ref_pos(rec,it->pos,&bq); if(code<0||bq<min_baseq) continue; char base=nt16_to_base(code); int allele=base==it->ref?0:(base==it->alt?1:-1); if(allele<0) continue; ++stats.observations; int si=(int)(&(*it)-panel.sites.data()); MoleculeKey key{cb,si,ub}; auto m=mol.find(key); if(m==mol.end()) mol.emplace(std::move(key),(uint8_t)allele); else if(m->second!=(uint8_t)allele) m->second=2; }
    }
    CellCounts out; for(const auto& kv:mol){ if(kv.second>1){++stats.conflicting_molecules;continue;} auto& c=out[kv.first.barcode][kv.first.site]; if(kv.second==0)++c.ref;else++c.alt; ++stats.molecules; }
    bam_destroy1(rec);hts_itr_destroy(itr);hts_idx_destroy(idx);bam_hdr_destroy(hdr);sam_close(fp);return out;
}

double binomial_log_pmf(uint64_t alt, uint64_t ref, double q) {
    const double n = static_cast<double>(alt + ref);
    const double k = static_cast<double>(alt);
    q = std::max(1e-12, std::min(1.0 - 1e-12, q));
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) -
           std::lgamma(n - k + 1.0) + k * std::log(q) +
           (n - k) * std::log1p(-q);
}

double beta_binomial_log_pmf(uint64_t alt, uint64_t ref, double q, double rho) {
    if (rho <= 1e-9 || alt + ref <= 1) return binomial_log_pmf(alt, ref, q);
    q = std::max(1e-12, std::min(1.0 - 1e-12, q));
    rho = std::max(1e-9, std::min(0.999999, rho));
    const double concentration = 1.0 / rho - 1.0;
    const double alpha = std::max(1e-12, q * concentration);
    const double beta = std::max(1e-12, (1.0 - q) * concentration);
    const double n = static_cast<double>(alt + ref);
    const double k = static_cast<double>(alt);
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) -
           std::lgamma(n - k + 1.0) +
           std::lgamma(k + alpha) + std::lgamma(n - k + beta) -
           std::lgamma(n + alpha + beta) -
           std::lgamma(alpha) - std::lgamma(beta) +
           std::lgamma(alpha + beta);
}

Score score_genotype(const Panel& panel,
                     const std::unordered_map<int,AlleleCount>& counts,
                     const std::string& genotype,
                     double err,
                     double rho) {
    Score out;
    auto comps = split_genotype(genotype);
    if (comps.empty() || comps.size() > 2) return out;
    std::vector<int> idx;
    for (const auto& c : comps) {
        auto it = panel.sample_index.find(c);
        if (it == panel.sample_index.end()) return out;
        idx.push_back(it->second);
    }
    if (idx.size() == 1) idx.push_back(idx[0]);

    std::unordered_set<int> used;
    for (const auto& kv : counts) {
        const int si = kv.first;
        if (si < 0 || si >= static_cast<int>(panel.sites.size())) continue;
        const int a = panel.sites[si].state[idx[0]];
        const int b = panel.sites[si].state[idx[1]];
        if (a < 0 || b < 0) continue;
        double p = (a + b) / 2.0;
        p = err + (1.0 - 2.0 * err) * p;
        out.log_likelihood += beta_binomial_log_pmf(kv.second.alt, kv.second.ref, p, rho);
        out.molecules += kv.second.total();
        used.insert(si);
    }
    out.sites = static_cast<int>(used.size());
    out.scoreable = out.molecules > 0;
    return out;
}

double fit_overdispersion_rho(const Panel& panel,
                              const std::unordered_map<int,AlleleCount>& counts,
                              const std::string& genotype,
                              double error_rate,
                              double rho_max) {
    rho_max = std::max(0.0, std::min(0.95, rho_max));
    if (rho_max <= 1e-9) return 0.0;
    double best_rho = 0.0;
    double best_ll = score_genotype(panel, counts, genotype, error_rate, 0.0).log_likelihood;
    // Deterministic coarse-to-fine one-dimensional profile fit.
    for (int pass = 0; pass < 3; ++pass) {
        const double half = pass == 0 ? rho_max : rho_max / std::pow(10.0, pass);
        const double lo = std::max(0.0, best_rho - half);
        const double hi = std::min(rho_max, best_rho + half);
        for (int i = 0; i <= 40; ++i) {
            const double rho = lo + (hi - lo) * static_cast<double>(i) / 40.0;
            const Score sc = score_genotype(panel, counts, genotype, error_rate, rho);
            if (sc.scoreable && sc.log_likelihood > best_ll) {
                best_ll = sc.log_likelihood;
                best_rho = rho;
            }
        }
    }
    return best_rho;
}

double fusion_parent2_fraction(const Panel& panel,const std::unordered_map<int,AlleleCount>& counts,const std::string& p1,const std::string& p2,uint64_t* nmol,int* nsites){
    auto a=panel.sample_index.find(p1),b=panel.sample_index.find(p2); if(a==panel.sample_index.end()||b==panel.sample_index.end()) return std::numeric_limits<double>::quiet_NaN(); uint64_t s1=0,s2=0;int sites=0;
    for(const auto& kv:counts){int si=kv.first;if(si<0||si>=(int)panel.sites.size())continue;int x=panel.sites[si].state[a->second],y=panel.sites[si].state[b->second];if(x<0||y<0||x==y)continue;++sites;if(x==0&&y==1){s1+=kv.second.ref;s2+=kv.second.alt;}else{s1+=kv.second.alt;s2+=kv.second.ref;}}
    if (nmol) *nmol = s1 + s2;
    if (nsites) *nsites = sites;
    if (s1 + s2 == 0) return std::numeric_limits<double>::quiet_NaN();
    return static_cast<double>(s2) / static_cast<double>(s1 + s2);
}

} // namespace mt_evidence
