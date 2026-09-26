// Bounded native containers for the targeted cache consumer only.
// Included after JointLinkedUnit is defined. Spill files live beside the task
// output on the configured BeeGFS/RAM filesystem and are removed on close.
#ifndef CELLBOUNCER_JOINT_DOUBLET_BOUNDED_H
#define CELLBOUNCER_JOINT_DOUBLET_BOUNDED_H
#include <memory>
#include <fcntl.h>
#include <mutex>

static string joint_spill_root;
static size_t joint_spill_limit = 16ULL*1024*1024;
static atomic<unsigned long> joint_spill_serial(0);

static void joint_pread_exact(int fd,void* destination,size_t bytes,uint64_t offset){
    char* out=static_cast<char*>(destination);
    while (bytes){
        const ssize_t got=pread(fd,out,bytes,static_cast<off_t>(offset));
        if (got<0 && errno==EINTR) continue;
        if (got<=0) throw runtime_error("truncated/failed bounded evidence read");
        out+=got; bytes-=got; offset+=got;
    }
}
static void joint_write_exact(int fd,const void* source,size_t bytes){
    const char* in=static_cast<const char*>(source);
    while (bytes){
        const ssize_t wrote=write(fd,in,bytes);
        if (wrote<0 && errno==EINTR) continue;
        if (wrote<=0) throw runtime_error("failed bounded evidence spill write");
        in+=wrote; bytes-=wrote;
    }
}
// Immutable raw slices never load an entire high-coverage cell. Each iterator
// owns at most 4096 records; copies are independent and safe across workers.
template<class Record> class JointRecordView {
    struct File { int fd=-1; ~File(){if(fd>=0) close(fd);} };
    shared_ptr<File> file;
    uint64_t first=0,n=0;
public:
    JointRecordView(){}
    JointRecordView(const string& path,uint64_t offset,uint64_t count):file(new File),first(offset),n(count){
        file->fd=open(path.c_str(),O_RDONLY);
        if(file->fd<0) throw runtime_error("could not open indexed evidence: "+path);
        struct stat st;
        if(fstat(file->fd,&st)!=0 || offset>uint64_t(st.st_size)/sizeof(Record) ||
           count>uint64_t(st.st_size)/sizeof(Record)-offset)
            throw runtime_error("indexed evidence slice exceeds file: "+path);
    }
    uint64_t size() const{return n;} bool empty() const{return n==0;}
    class iterator {
        const JointRecordView* view; uint64_t i;
        mutable vector<Record> block; mutable uint64_t start=UINT64_MAX;
    public:
        iterator(const JointRecordView* v,uint64_t p):view(v),i(p){}
        bool operator!=(const iterator& other)const{return i!=other.i;}
        iterator& operator++(){++i;return *this;}
        const Record& operator*()const{
            if(start==UINT64_MAX || i<start || i>=start+block.size()){
                start=i; block.resize(static_cast<size_t>(min<uint64_t>(4096,view->n-i)));
                joint_pread_exact(view->file->fd,block.data(),block.size()*sizeof(Record),(view->first+i)*sizeof(Record));
            }
            return block[static_cast<size_t>(i-start)];
        }
    };
    iterator begin()const{return iterator(this,0);} iterator end()const{return iterator(this,n);}
};

// Observed fields are interned once per active cell/control-parent workspace.
// Spilled candidate terms contain compact references and candidate coefficients,
// not copies of counts, genomic labels or background/locked probabilities.
struct JointObservedPool {
    typedef tuple<int32_t,int32_t,string,string,string,long double,long double,
                  long double,long double,long double> Key;
    map<Key,uint64_t> ids;
    vector<JointSiteUnit> observed;
    mutex lock;
    uint64_t intern(const JointSiteUnit& v){
        lock_guard<mutex> guard(lock);
        const Key key=make_tuple(v.tid,v.pos,v.contig,v.ref_allele,v.alt_allele,
            v.ref,v.alt,v.q_locked,v.q_ambient,isfinite(v.probability_intercept)?v.probability_intercept:-1.0L);
        auto found=ids.find(key);if(found!=ids.end())return found->second;
        const uint64_t id=observed.size();JointSiteUnit raw=v;
        raw.q_second=0;raw.probability_slope=0;raw.compiled_site_id=UINT32_MAX;
        observed.push_back(move(raw));ids.emplace(key,id);return id;
    }
    JointSiteUnit get(uint64_t id){
        lock_guard<mutex> guard(lock);
        if(id>=observed.size())throw runtime_error("invalid observed-evidence reference");
        return observed[id];
    }
};
static shared_ptr<JointObservedPool> joint_current_observed_pool;

// Candidate workspaces spill whole linked units, never individual SNP draws.
// No candidate observations accumulate in the library-wide result queue.
class JointStoredUnits {
    struct State {
        shared_ptr<JointObservedPool> pool;
        State():pool(joint_current_observed_pool?joint_current_observed_pool:shared_ptr<JointObservedPool>(new JointObservedPool)){}
        vector<JointLinkedUnit> memory;
        vector<uint64_t> offsets;
        size_t bytes=0;
        int fd=-1;
        uint64_t end=0;
        ~State(){if(fd>=0) close(fd);}
    };
    shared_ptr<State> state;
    template<class T> static void put(string& data,const T& v){data.append(reinterpret_cast<const char*>(&v),sizeof(v));}
    static void put_string(string& data,const string& v){uint64_t n=v.size();put(data,n);data+=v;}
    template<class T> static void get(const string& data,size_t& p,T& v){
        if(p>data.size() || sizeof(v)>data.size()-p)throw runtime_error("invalid unit spill record");
        memcpy(&v,data.data()+p,sizeof(v));p+=sizeof(v);
    }
    static string get_string(const string& data,size_t& p){
        uint64_t n;get(data,p,n);
        if(p>data.size() || n>data.size()-p)throw runtime_error("invalid unit spill string");
        string v=data.substr(p,n);p+=n;return v;
    }
    string pack(const JointLinkedUnit& u){
        string data;put(data,u.molecule);put_string(data,u.molecule_id);put(data,u.basis);
        put_string(data,u.parent_origin);put(data,u.compiled_unit_id);put(data,u.genotype_distinguishable);put(data,u.fold);
        uint64_t n=u.sites.size();put(data,n);
        for(const JointSiteUnit& v:u.sites){
            const uint64_t reference=state->pool->intern(v);put(data,reference);
            put(data,v.compiled_site_id);put(data,v.q_second);put(data,v.probability_slope);
        }return data;
    }
    JointLinkedUnit unpack(const string& data)const{
        size_t p=0;JointLinkedUnit u;get(data,p,u.molecule);u.molecule_id=get_string(data,p);get(data,p,u.basis);
        u.parent_origin=get_string(data,p);get(data,p,u.compiled_unit_id);get(data,p,u.genotype_distinguishable);get(data,p,u.fold);
        uint64_t n;get(data,p,n);
        if(n>data.size()/32)throw runtime_error("invalid spilled unit site count");
        u.sites.resize(n);
        for(JointSiteUnit& v:u.sites){
            uint64_t reference;get(data,p,reference);v=state->pool->get(reference);
            get(data,p,v.compiled_site_id);get(data,p,v.q_second);get(data,p,v.probability_slope);
        }
        if(p!=data.size())throw runtime_error("trailing unit spill bytes");return u;
    }
    void append_disk(const JointLinkedUnit& u){
        string data=pack(u);uint64_t bytes=data.size();state->offsets.push_back(state->end);
        joint_write_exact(state->fd,&bytes,sizeof(bytes));joint_write_exact(state->fd,data.data(),bytes);
        state->end+=sizeof(bytes)+bytes;
    }
    void spill(){
        if(state->fd>=0)return;
        if(joint_spill_root.empty())throw runtime_error("bounded unit spill root is unset");
        const string path=joint_spill_root+".units."+to_string((long long)getpid())+"."+to_string(joint_spill_serial++);
        state->fd=open(path.c_str(),O_CREAT|O_EXCL|O_RDWR,0600);
        if(state->fd<0)throw runtime_error("cannot create bounded unit spill: "+path);
        if(unlink(path.c_str())!=0)throw runtime_error("cannot unlink task-owned unit spill");
        for(const JointLinkedUnit& u:state->memory)append_disk(u);
        vector<JointLinkedUnit>().swap(state->memory);
    }
public:
    JointStoredUnits():state(new State){}
    JointStoredUnits(const vector<JointLinkedUnit>& units):state(new State){for(const auto& u:units)push_back(u);}
    void reserve(size_t){} // Capacity is governed by bytes, never the cell count.
    size_t size()const{return state->fd<0?state->memory.size():state->offsets.size();}
    bool empty()const{return size()==0;}
    void push_back(const JointLinkedUnit& u){
        if(!state.unique())throw runtime_error("attempt to mutate shared immutable unit store");
        size_t bytes=sizeof(u)+u.sites.size()*sizeof(JointSiteUnit)+u.molecule_id.size()+u.parent_origin.size();
        for(const auto& v:u.sites)bytes+=v.contig.size()+v.ref_allele.size()+v.alt_allele.size();
        if(state->fd<0 && state->bytes+bytes>joint_spill_limit)spill();
        if(state->fd<0)state->memory.push_back(u);else append_disk(u);
        state->bytes+=bytes;
    }
    JointLinkedUnit operator[](size_t i)const{
        if(i>=size())throw runtime_error("bounded unit index out of range");
        if(state->fd<0)return state->memory[i];
        uint64_t bytes;const uint64_t p=state->offsets[i];joint_pread_exact(state->fd,&bytes,sizeof(bytes),p);
        if(bytes>state->end-p-sizeof(bytes))throw runtime_error("invalid spilled unit length");
        string data(bytes,'\0');joint_pread_exact(state->fd,&data[0],bytes,p+sizeof(bytes));return unpack(data);
    }
    class iterator {
        const JointStoredUnits* view;size_t i;mutable JointLinkedUnit unit;
    public:
        iterator(const JointStoredUnits* v,size_t p):view(v),i(p){}
        bool operator!=(const iterator& other)const{return i!=other.i;}
        iterator& operator++(){++i;return *this;}
        const JointLinkedUnit& operator*()const{
            if(view->state->fd<0)return view->state->memory[i];
            unit=(*view)[i];return unit;
        }
    };
    iterator begin()const{return iterator(this,0);}iterator end()const{return iterator(this,size());}
};
#endif
