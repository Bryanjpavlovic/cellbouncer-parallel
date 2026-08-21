#include "mt_evidence.h"

#include <zlib.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {
struct Options {
    std::string bam, vcf, candidate_manifest, library_id, site_manifest;
    std::string haplotype_groups, haplotype_pairwise, output_prefix;
    std::string mito_chrom="chrM", barcode_tag="CB", umi_tag="UB";
    int min_mapq=30,min_baseq=20,min_depth=10,min_molecules=10,min_sites=2,threads=4;
    double homoplasmy_af=0.95,error_rate=0.005;
};
struct Candidate { std::string library,barcode,hypothesis_id,state_notation,donor,current; };

std::string trim(std::string s){size_t a=s.find_first_not_of(" \t\r\n");if(a==std::string::npos)return"";size_t b=s.find_last_not_of(" \t\r\n");return s.substr(a,b-a+1);}
std::vector<std::string> split(const std::string&s,char d){std::vector<std::string>o;std::stringstream ss(s);std::string x;while(std::getline(ss,x,d))o.push_back(x);return o;}
std::string fmt(double x){if(!std::isfinite(x))return"NA";std::ostringstream o;o<<std::setprecision(12)<<x;return o.str();}
std::string canonical_library_key(const std::string& raw){
    std::string s=trim(raw);
    if(s.size()>3 && (s[0]=='l'||s[0]=='L') && (s[1]=='i'||s[1]=='I') && (s[2]=='b'||s[2]=='B')){
        bool numeric=true;for(size_t i=3;i<s.size();++i)if(s[i]<'0'||s[i]>'9'){numeric=false;break;}
        if(numeric)s=s.substr(3);
    }
    return s;
}
bool same_library_key(const std::string& a,const std::string& b){
    return a==b || canonical_library_key(a)==canonical_library_key(b);
}

void usage(FILE* f,int code){
    fprintf(f,"mt_identity_score --bam RNA.bam --vcf MT.bcf --candidate_manifest FILE --library_id libN --site_manifest FILE --haplotype_pairwise FILE --output_prefix PREFIX [options]\n"
              "Targeted post-hoc mitochondrial donor verification for proposed singlet shifts and orthogonal component support of disjoint real biological lines. No ambient-RNA inputs are accepted.\n"
              "If a library has no fusion-site manifest rows but has donor pairwise metadata, scoring automatically falls back to the retained MT panel in donor-only mode.\n"
              "  --haplotype_groups FILE   optional provenance sidecar\n"
              "  --mito_chrom NAME         [chrM]\n  --min_depth N [10]\n  --homoplasmy_af X [0.95]\n"
              "  --min_molecules N [10]\n  --min_sites N [2]\n  --error_rate X [0.005]\n"
              "  --min_mapq N [30]\n  --min_baseq N [20]\n  --threads N [4]\n");
    exit(code);
}
Options parse(int argc,char**argv){
    Options o; static option opts[]={
        {"bam",1,nullptr,'b'},{"vcf",1,nullptr,'v'},{"candidate_manifest",1,nullptr,'c'},{"library_id",1,nullptr,'l'},
        {"site_manifest",1,nullptr,1000},{"haplotype_groups",1,nullptr,1001},{"haplotype_pairwise",1,nullptr,1002},{"output_prefix",1,nullptr,'o'},
        {"mito_chrom",1,nullptr,1003},{"barcode_tag",1,nullptr,1004},{"umi_tag",1,nullptr,1005},{"min_mapq",1,nullptr,1006},{"min_baseq",1,nullptr,1007},
        {"min_depth",1,nullptr,1008},{"homoplasmy_af",1,nullptr,1009},{"min_molecules",1,nullptr,1010},{"min_sites",1,nullptr,1011},{"error_rate",1,nullptr,1012},{"threads",1,nullptr,'t'},{"help",0,nullptr,'h'},{nullptr,0,nullptr,0}};
    int ch;while((ch=getopt_long(argc,argv,"b:v:c:l:o:t:h",opts,nullptr))!=-1){switch(ch){
        case'b':o.bam=optarg;break;case'v':o.vcf=optarg;break;case'c':o.candidate_manifest=optarg;break;case'l':o.library_id=optarg;break;case'o':o.output_prefix=optarg;break;case't':o.threads=atoi(optarg);break;case'h':usage(stdout,0);break;
        case 1000:o.site_manifest=optarg;break;case 1001:o.haplotype_groups=optarg;break;case 1002:o.haplotype_pairwise=optarg;break;case 1003:o.mito_chrom=optarg;break;case 1004:o.barcode_tag=optarg;break;case 1005:o.umi_tag=optarg;break;case 1006:o.min_mapq=atoi(optarg);break;case 1007:o.min_baseq=atoi(optarg);break;case 1008:o.min_depth=atoi(optarg);break;case 1009:o.homoplasmy_af=atof(optarg);break;case 1010:o.min_molecules=atoi(optarg);break;case 1011:o.min_sites=atoi(optarg);break;case 1012:o.error_rate=atof(optarg);break;default:usage(stderr,2);}}
    if(o.bam.empty()||o.vcf.empty()||o.candidate_manifest.empty()||o.library_id.empty()||o.site_manifest.empty()||o.haplotype_pairwise.empty()||o.output_prefix.empty())usage(stderr,2);
    if (o.barcode_tag.size()!=2 || o.umi_tag.size()!=2 || o.min_depth<1 ||
        o.min_sites<1 || o.min_molecules<0 || o.threads<1 ||
        o.homoplasmy_af<=0.5 || o.homoplasmy_af>1 ||
        o.error_rate<=0 || o.error_rate>=0.5) {
        throw std::runtime_error("invalid option value");
    }
    return o;
}

std::vector<Candidate> load_candidates(const std::string& path,const std::string& lib,std::unordered_set<std::string>& barcodes){
    gzFile g=gzopen(path.c_str(),"rb");if(!g)throw std::runtime_error("Could not open candidate manifest: "+path);char buf[1<<20];if(!gzgets(g,buf,sizeof(buf)))throw std::runtime_error("Empty candidate manifest");auto h=split(trim(buf),'\t');std::map<std::string,int>ix;for(int i=0;i<(int)h.size();++i)ix[h[i]]=i;for(auto n:{"library","barcode","hypothesis_id","state_notation","donor_genotype","current_donor_genotype"})if(!ix.count(n))throw std::runtime_error(std::string("Candidate manifest missing column: ")+n);
    std::vector<Candidate> all;
    std::unordered_map<std::string,bool> has_identity_shift;
    std::unordered_map<std::string,std::set<std::string>> donors_by_barcode;
    std::unordered_map<std::string,std::string> current_by_barcode;
    std::unordered_map<std::string,std::string> library_by_barcode;
    while(gzgets(g,buf,sizeof(buf))){
        auto f=split(trim(buf),'\t');if(f.size()!=h.size())continue;if(!same_library_key(f[ix["library"]],lib))continue;
        Candidate c{f[ix["library"]],f[ix["barcode"]],f[ix["hypothesis_id"]],f[ix["state_notation"]],f[ix["donor_genotype"]],f[ix["current_donor_genotype"]]};
        all.push_back(c);
        current_by_barcode[c.barcode]=c.current;
        library_by_barcode[c.barcode]=c.library;
        if(c.donor!=c.current)has_identity_shift[c.barcode]=true;
        for(const auto& d:mt_evidence::split_genotype(c.current))if(!d.empty())donors_by_barcode[c.barcode].insert(d);
        for(const auto& d:mt_evidence::split_genotype(c.donor))if(!d.empty())donors_by_barcode[c.barcode].insert(d);
    }
    gzclose(g);

    // Preserve the original singlet-hypothesis rows because reconciliation uses
    // their hypothesis IDs for direct singlet verification.  In addition, score
    // every donor component appearing in the current or candidate biological
    // identities.  These donor-level probes let reconciliation ask whether an
    // independently proposed real A+B line carries mitochondria from A or B
    // without ever scoring or constructing a compound mitochondrial genotype.
    std::vector<Candidate> out;
    std::unordered_map<std::string,std::set<std::string>> represented_donors;
    for(const auto& c:all){
        if(!has_identity_shift[c.barcode])continue;
        auto dc=mt_evidence::split_genotype(c.donor);
        if(dc.size()!=1)continue;
        out.push_back(c);
        represented_donors[c.barcode].insert(c.donor);
        barcodes.insert(c.barcode);
    }
    for(const auto& kv:donors_by_barcode){
        const std::string& bc=kv.first;
        if(!has_identity_shift[bc])continue;
        barcodes.insert(bc);
        for(const auto& donor:kv.second){
            if(represented_donors[bc].count(donor))continue;
            Candidate c{library_by_barcode[bc],bc,"MT_COMPONENT:"+donor,"D["+donor+"]",donor,current_by_barcode[bc]};
            out.push_back(c);
            represented_donors[bc].insert(donor);
        }
    }
    return out;
}

std::unordered_set<int64_t> manifest_positions(const std::string& path,const std::string& lib){
    std::ifstream in(path);if(!in)throw std::runtime_error("Could not open site manifest: "+path);std::string line;if(!std::getline(in,line))throw std::runtime_error("Empty site manifest");auto h=split(trim(line),'\t');std::map<std::string,int>ix;for(int i=0;i<(int)h.size();++i)ix[h[i]]=i;if(!ix.count("library_id")||!ix.count("pos"))throw std::runtime_error("Site manifest missing library_id/pos");std::unordered_set<int64_t> p;while(std::getline(in,line)){if(line.empty())continue;auto f=split(trim(line),'\t');if(f.size()!=h.size()||!same_library_key(f[ix["library_id"]],lib))continue;long long pos=atoll(f[ix["pos"]].c_str());if(pos>0)p.insert(pos-1);}return p;
}

struct PairResolution { int distinguishing=-1; std::string status; };
std::unordered_map<std::string,PairResolution> load_pairwise(const std::string&path,const std::string&lib){
    std::ifstream in(path);if(!in)throw std::runtime_error("Could not open haplotype pairwise: "+path);std::string line;if(!std::getline(in,line))throw std::runtime_error("Empty haplotype pairwise");auto h=split(trim(line),'\t');std::map<std::string,int>ix;for(int i=0;i<(int)h.size();++i)ix[h[i]]=i;for(auto n:{"library_id","donor1","donor2","distinguishing_sites","status"})if(!ix.count(n))throw std::runtime_error(std::string("Haplotype pairwise missing: ")+n);std::unordered_map<std::string,PairResolution> out;while(std::getline(in,line)){auto f=split(trim(line),'\t');if(f.size()!=h.size()||!same_library_key(f[ix["library_id"]],lib))continue;out[mt_evidence::canonical_pair_key(f[ix["donor1"]],f[ix["donor2"]])]={atoi(f[ix["distinguishing_sites"]].c_str()),f[ix["status"]]};}return out;
}

bool unresolved_transition(const std::string& a,const std::string& b,const std::unordered_map<std::string,PairResolution>&pr){
    auto ca=mt_evidence::split_genotype(a),cb=mt_evidence::split_genotype(b);std::set<std::string> sa(ca.begin(),ca.end()),sb(cb.begin(),cb.end()),diff;std::set_symmetric_difference(sa.begin(),sa.end(),sb.begin(),sb.end(),std::inserter(diff,diff.begin()));
    if(diff.size()==2){auto it=diff.begin();std::string x=*it++;std::string y=*it;auto r=pr.find(mt_evidence::canonical_pair_key(x,y));if(r!=pr.end()&&(r->second.distinguishing<=0||r->second.status.find("UNRESOL")!=std::string::npos||r->second.status.find("IDENTICAL")!=std::string::npos))return true;}
    auto paircheck=[&](const std::vector<std::string>&c){if(c.size()!=2||c[0]==c[1])return false;auto r=pr.find(mt_evidence::canonical_pair_key(c[0],c[1]));return r!=pr.end()&&(r->second.distinguishing<=0||r->second.status.find("UNRESOL")!=std::string::npos||r->second.status.find("IDENTICAL")!=std::string::npos);};return paircheck(ca)||paircheck(cb);
}
}

int main(int argc,char**argv){
    try{
        Options o=parse(argc,argv);std::unordered_set<std::string> target_barcodes;auto candidates=load_candidates(o.candidate_manifest,o.library_id,target_barcodes);
        if(candidates.empty()){
            std::string score_path=o.output_prefix+".mt_identity_scores.tsv.gz";gzFile out=gzopen(score_path.c_str(),"wb");if(!out)throw std::runtime_error("Could not write "+score_path);
            gzprintf(out,"library\tbarcode\thypothesis_id\tstate_notation\tdonor_genotype\tmt_log_likelihood\tmt_delta_hypothesis_vs_current\tmt_delta_vs_best_other_singlet\tmt_rank_within_singlet_candidates\tmt_sites_used\tmt_molecules_used\tmt_best_identity\tmt_second_identity\tmt_haplotype_resolution\tmt_fit_status\tschema_version\n");gzclose(out);
            std::string fusion_path=o.output_prefix+".mt_fusion_scores.tsv.gz";gzFile fo=gzopen(fusion_path.c_str(),"wb");if(!fo)throw std::runtime_error("Could not write "+fusion_path);gzprintf(fo,"library\tbarcode\thypothesis_id\tparent1\tparent2\tparent2_fraction\tinformative_sites\tinformative_molecules\tfit_status\tschema_version\n");gzclose(fo);
            std::string cache=o.output_prefix+".mt_candidate_counts.tsv.gz";gzFile co=gzopen(cache.c_str(),"wb");if(!co)throw std::runtime_error("Could not write "+cache);gzprintf(co,"library\tbarcode\tchrom\tpos\tref\talt\tref_molecules\talt_molecules\n");gzclose(co);
            std::ofstream qc(o.output_prefix+".mt_qc.tsv");qc<<"metric\tvalue\n"<<"library_id\t"<<o.library_id<<"\n"<<"candidate_cells\t0\n"<<"panel_sites\t0\n"<<"reads_seen\t0\n"<<"reads_accepted\t0\n"<<"molecules\t0\n"<<"conflicting_molecules\t0\n"<<"ambient_rna_used\tno\n";
            return 0;
        }
        auto allowed=manifest_positions(o.site_manifest,o.library_id);
        auto pairwise=load_pairwise(o.haplotype_pairwise,o.library_id);
        const bool donor_only_mode=allowed.empty();
        if(donor_only_mode&&pairwise.empty())
            throw std::runtime_error("No fusion-site manifest rows or donor pairwise metadata for "+o.library_id);
        if(donor_only_mode)
            std::cerr<<"INFO: "<<o.library_id<<" has no fusion-site manifest rows; using donor-only MT identity mode over the retained mitochondrial panel\n";
        const std::unordered_set<int64_t>* allowed_positions=donor_only_mode?nullptr:&allowed;
        auto panel=mt_evidence::load_panel(o.vcf,o.mito_chrom,o.min_depth,o.homoplasmy_af,allowed_positions);
        mt_evidence::CountStats stats;auto counts=mt_evidence::count_bam_once(o.bam,panel,target_barcodes,o.barcode_tag,o.umi_tag,o.min_mapq,o.min_baseq,o.threads,stats);
        std::vector<mt_evidence::Score> scores(candidates.size());std::unordered_map<std::string,mt_evidence::Score> current_scores;std::unordered_map<std::string,std::vector<size_t>> by_bc;
        for(size_t i=0;i<candidates.size();++i){auto hit=counts.find(candidates[i].barcode);static const std::unordered_map<int,mt_evidence::AlleleCount> empty;const auto&cc=hit==counts.end()?empty:hit->second;scores[i]=mt_evidence::score_genotype(panel,cc,candidates[i].donor,o.error_rate);by_bc[candidates[i].barcode].push_back(i);auto curc=mt_evidence::split_genotype(candidates[i].current);if(curc.size()==1&&candidates[i].donor==candidates[i].current)current_scores[candidates[i].barcode]=scores[i];}
        std::vector<int> ranks(candidates.size(),0);std::vector<double> delta_other(candidates.size(),std::numeric_limits<double>::quiet_NaN());std::vector<std::string> competitor(candidates.size());std::unordered_map<std::string,std::pair<std::string,std::string>> best2;
        for(auto&kv:by_bc){
            auto ids=kv.second;
            std::sort(ids.begin(),ids.end(),[&](size_t a,size_t b){return scores[a].log_likelihood>scores[b].log_likelihood;});
            int rank=0;std::vector<size_t> usable;
            for(size_t i:ids)if(scores[i].scoreable&&scores[i].molecules>=(uint64_t)o.min_molecules&&scores[i].sites>=o.min_sites){ranks[i]=++rank;usable.push_back(i);}
            if(!usable.empty())best2[kv.first].first=candidates[usable[0]].donor;
            if(usable.size()>1)best2[kv.first].second=candidates[usable[1]].donor;
            for(size_t i:usable){
                size_t other=(size_t)-1;
                for(size_t j:usable){if(j!=i){other=j;break;}}
                if(other!=(size_t)-1){delta_other[i]=scores[i].log_likelihood-scores[other].log_likelihood;competitor[i]=candidates[other].donor;}
            }
        }
        std::string score_path=o.output_prefix+".mt_identity_scores.tsv.gz";gzFile out=gzopen(score_path.c_str(),"wb");if(!out)throw std::runtime_error("Could not write "+score_path);gzprintf(out,"library\tbarcode\thypothesis_id\tstate_notation\tdonor_genotype\tmt_log_likelihood\tmt_delta_hypothesis_vs_current\tmt_delta_vs_best_other_singlet\tmt_rank_within_singlet_candidates\tmt_sites_used\tmt_molecules_used\tmt_best_identity\tmt_second_identity\tmt_haplotype_resolution\tmt_fit_status\tschema_version\n");
        for(size_t i=0;i<candidates.size();++i){
            auto cur=current_scores.find(candidates[i].barcode);bool enough=scores[i].scoreable&&scores[i].molecules>=(uint64_t)o.min_molecules&&scores[i].sites>=o.min_sites;
            double delta_current=(enough&&cur!=current_scores.end()&&cur->second.scoreable)?scores[i].log_likelihood-cur->second.log_likelihood:std::numeric_limits<double>::quiet_NaN();
            bool has_comp=!competitor[i].empty();bool unr=has_comp&&unresolved_transition(candidates[i].donor,competitor[i],pairwise);
            std::string status=!scores[i].scoreable?"MT_UNAVAILABLE":(!enough?"MT_INSUFFICIENT":(!has_comp?"MT_NO_SINGLET_COMPARATOR":"PASS"));
            std::string res=!has_comp?"MT_HAPLOTYPE_UNRESOLVED":(unr?"MT_HAPLOTYPE_UNRESOLVED":"MT_HAPLOTYPE_RESOLVED");auto b=best2[candidates[i].barcode];
            gzprintf(out,"%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%d\t%llu\t%s\t%s\t%s\t%s\tmt_identity_scores_v3_donor_probe\n",candidates[i].library.c_str(),candidates[i].barcode.c_str(),candidates[i].hypothesis_id.c_str(),candidates[i].state_notation.c_str(),candidates[i].donor.c_str(),fmt(enough?scores[i].log_likelihood:std::numeric_limits<double>::quiet_NaN()).c_str(),fmt(unr?std::numeric_limits<double>::quiet_NaN():delta_current).c_str(),fmt((unr||!has_comp)?std::numeric_limits<double>::quiet_NaN():delta_other[i]).c_str(),ranks[i],scores[i].sites,(unsigned long long)scores[i].molecules,b.first.c_str(),b.second.c_str(),res.c_str(),status.c_str());
        }
        gzclose(out);
        // Fusion mitochondrial composition is deliberately outside identity reconciliation.
        // The separate mt_fusion_ratio executable owns that analysis after nuclear identity is fixed.
        std::string fusion_path=o.output_prefix+".mt_fusion_scores.tsv.gz";gzFile fo=gzopen(fusion_path.c_str(),"wb");if(!fo)throw std::runtime_error("Could not write "+fusion_path);gzprintf(fo,"library\tbarcode\thypothesis_id\tparent1\tparent2\tparent2_fraction\tinformative_sites\tinformative_molecules\tfit_status\tschema_version\n");gzclose(fo);
        std::string cache=o.output_prefix+".mt_candidate_counts.tsv.gz";gzFile co=gzopen(cache.c_str(),"wb");gzprintf(co,"library\tbarcode\tchrom\tpos\tref\talt\tref_molecules\talt_molecules\n");for(const auto&cell:counts)for(const auto&kv:cell.second){const auto&s=panel.sites[kv.first];gzprintf(co,"%s\t%s\t%s\t%lld\t%c\t%c\t%llu\t%llu\n",o.library_id.c_str(),cell.first.c_str(),panel.chrom.c_str(),(long long)(s.pos+1),s.ref,s.alt,(unsigned long long)kv.second.ref,(unsigned long long)kv.second.alt);}gzclose(co);
        std::ofstream qc(o.output_prefix+".mt_qc.tsv");qc<<"metric\tvalue\n"<<"library_id\t"<<o.library_id<<"\n"<<"candidate_cells\t"<<target_barcodes.size()<<"\n"<<"site_selection_mode\t"<<(donor_only_mode?"DONOR_ONLY_RETAINED_PANEL":"FUSION_SITE_MANIFEST")<<"\n"<<"fusion_manifest_positions\t"<<allowed.size()<<"\n"<<"pairwise_donor_comparisons\t"<<pairwise.size()<<"\n"<<"panel_sites\t"<<panel.sites.size()<<"\n"<<"reads_seen\t"<<stats.reads_seen<<"\n"<<"reads_accepted\t"<<stats.reads_accepted<<"\n"<<"molecules\t"<<stats.molecules<<"\n"<<"conflicting_molecules\t"<<stats.conflicting_molecules<<"\n"<<"ambient_rna_used\tno\n";
        return 0;
    }catch(const std::exception&e){fprintf(stderr,"ERROR: %s\n",e.what());return 1;}
}
