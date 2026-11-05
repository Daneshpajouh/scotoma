// Exact ball convolution, with a radial constant block and ternary coordinates.
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using Big = boost::multiprecision::cpp_int;
using Wide = __int128_t;
using Clock = std::chrono::steady_clock;

uint64_t choose(int n, int k) {
    if (k < 0 || k > n) return 0;
    uint64_t a = 1;
    for (int j=1; j<=k; ++j) a=a*(n-j+1)/j;
    return a;
}
uint64_t power(uint64_t a, int n) { uint64_t b=1; while(n--) b*=a; return b; }
Big natural(const std::string& s) {
    if(s.empty() || s.size()>100000) throw std::runtime_error("bad integer length");
    Big x=0;
    for(char c:s) {
        if(c<'0' || c>'9') throw std::runtime_error("expected nonnegative decimal integer");
        x=x*10+(c-'0');
    }
    return x;
}
std::string decimal(const Big& x) { return x.str(); }
std::string decimal(Wide x) {
    if(!x) return "0";
    bool neg=x<0;
    __uint128_t u=neg ? __uint128_t(-(x+1))+1 : __uint128_t(x);
    std::string s;
    while(u) { s.push_back(char('0'+u%10)); u/=10; }
    if(neg) s.push_back('-');
    std::reverse(s.begin(),s.end());
    return s;
}
Wide to_wide(const Big& b) {
    // Bound is checked in arbitrary precision before any conversion.
    return b.convert_to<Wide>();
}
template<class I> I convert(const Big& b);
template<> Big convert<Big>(const Big& b) { return b; }
template<> Wide convert<Wide>(const Big& b) { return to_wide(b); }

template<class F> void parallel(unsigned n, F f) {
    if(n==1) { f(0); return; }
    std::vector<std::thread> pool;
    for(unsigned t=1;t<n;++t) pool.emplace_back([&,t]{f(t);});
    f(0);
    for(auto& t:pool) t.join();
}
struct Term { uint64_t id; int w, nx; Big q; };
struct Input {
    int N, r, m, L, C;
    uint64_t P, count;
    Big D, M, bound;
    std::vector<Term> terms;
};
Input read_input(const std::string& path) {
    std::ifstream f(path);
    std::string magic, den;
    int nterms;
    Input a;
    if(!(f>>magic>>a.N>>a.r>>a.m>>den>>nterms) || magic!="ZX2_DENSITY_V1" ||
       a.N<1 || a.N>22 || a.r<1 || a.r>a.N+1 || a.m<0 || a.m>a.N ||
       nterms<1 || nterms>200000) throw std::runtime_error("bad density header");
    a.D=natural(den);
    if(a.D<=0) throw std::runtime_error("zero denominator");
    a.L=a.r-1; a.C=a.N-a.L;
    a.P=power(3,a.L); a.count=a.P*(a.C+1);
    a.M=0;
    for(int k=0;k<nterms;++k) {
        std::string s, q;
        int w;
        if(!(f>>s>>w>>q) || (int(s.size())!=a.L && !(a.L==0 && s=="-")) || w<0 || w>a.C)
            throw std::runtime_error("bad density term");
        // A dash represents the empty prefix, for r=1.
        if(a.L==0 && s=="-") s.clear();
        uint64_t id=0, stride=1;
        int nx=0, best=a.N+1;
        for(char c:s) {
            int d=c=='X'?2:c-'0';
            if(d<0 || d>2) throw std::runtime_error("bad prefix symbol");
            id+=d*stride; stride*=3; nx+=d==2;
        }
        for(int j=0;j<a.r;++j) {
            int d=w;
            for(int i=0;i<a.L;++i) d+=s[i]!=(i<j?'1':'0');
            best=std::min(best,d);
        }
        if(best>a.m) throw std::runtime_error("positive support outside U");
        Big weight=natural(q);
        if(weight<=0) throw std::runtime_error("nonpositive density");
        a.M+=weight*power(2,nx)*choose(a.C,w)*power(3,w);
        a.terms.push_back({id,w,nx,weight});
    }
    std::string tail;
    if(f>>tail) throw std::runtime_error("trailing density input");
    uint64_t volume=0;
    for(int d=0;d<=a.m;++d) volume+=choose(a.N,d)*power(3,d);
    a.bound=a.M*volume*std::max(uint64_t(4),power(4,a.C));
    return a;
}

// K[s][w] is the Fourier coefficient of a weight-w sphere at frequency weight s.
std::vector<std::vector<int64_t>> radial_matrix(int C) {
    std::vector<std::vector<int64_t>> K(C+1,std::vector<int64_t>(C+1));
    for(int s=0;s<=C;++s) for(int w=0;w<=C;++w)
        for(int j=0;j<=w;++j) {
            int64_t x=int64_t(choose(s,j)*choose(C-s,w-j)*power(3,w-j));
            K[s][w]+=(j%2?-x:x);
        }
    // A cheap exact check of the complete radial inverse, every run.
    for(int i=0;i<=C;++i) for(int j=0;j<=C;++j) {
        Wide z=0;
        for(int k=0;k<=C;++k) z+=Wide(K[i][k])*K[k][j];
        if(z!=(i==j?Wide(power(4,C)):Wide(0)))
            throw std::runtime_error("radial inverse identity failed");
    }
    return K;
}

template<class I> void prefix_transform(std::vector<I>& a, int L,
                                       unsigned nt, bool inverse) {
    const uint64_t triples=a.size()/3;
    uint64_t stride=1;
    for(int axis=0;axis<L;++axis,stride*=3) {
        std::atomic<bool> bad{false};
        parallel(nt,[&](unsigned t){
            uint64_t pos=triples*t/nt, end=triples*(t+1)/nt;
            while(pos<end) {
                uint64_t group=pos/stride, off=pos%stride;
                uint64_t len=std::min(end-pos,stride-off), base=group*3*stride+off;
                for(uint64_t j=0;j<len;++j) {
                    uint64_t i=base+j;
                    const I x=a[i], y=a[i+stride], z=a[i+2*stride];
                    if(!inverse) {
                        a[i]=x+y+2*z; a[i+stride]=x-y; a[i+2*stride]=x+y-2*z;
                    } else {
                        const I p=x+2*y+z, q=x-2*y+z, v=x-z;
                        if(p%4!=0 || q%4!=0 || v%4!=0) bad=true;
                        a[i]=p/4; a[i+stride]=q/4; a[i+2*stride]=v/4;
                    }
                }
                pos+=len;
            }
        });
        if(bad) throw std::runtime_error("nonintegral prefix inverse");
        std::cerr<<(inverse?"inverse":"forward")<<" prefix axis "<<axis+1<<"/"<<L<<std::endl;
    }
}

template<class I> void suffix_transform(std::vector<I>& a, const Input& in,
                                        unsigned nt, bool inverse,
                                        const std::vector<std::vector<int64_t>>& K) {
    std::atomic<bool> bad{false};
    parallel(nt,[&](unsigned t){
        std::array<I,23> row, result;
        uint64_t lo=in.P*t/nt, hi=in.P*(t+1)/nt;
        auto tick=Clock::now();
        for(uint64_t id=lo;id<hi;++id) {
            for(int w=0;w<=in.C;++w) row[w]=a[w*in.P+id];
            for(int s=0;s<=in.C;++s) {
                result[s]=0;
                for(int w=0;w<=in.C;++w) result[s]+=row[w]*K[s][w];
                if(inverse) {
                    if(result[s]%power(4,in.C)!=0) bad=true;
                    result[s]/=power(4,in.C);
                }
            }
            for(int s=0;s<=in.C;++s) a[s*in.P+id]=result[s];
            if(t==0 && (id-lo)%1048576==0 && Clock::now()-tick>std::chrono::seconds(10)) {
                std::cerr<<(inverse?"inverse":"forward")<<" suffix thread0 "<<id-lo<<"/"<<hi-lo<<std::endl;
                tick=Clock::now();
            }
        }
    });
    if(bad) throw std::runtime_error("nonintegral radial inverse");
    std::cerr<<(inverse?"inverse":"forward")<<" constant-block transform complete"<<std::endl;
}

// Enumerate ternary digits with an odometer, avoiding division at every entry.
struct Digits {
    std::array<int,22> d{};
    int L, nonzero=0, nx=0;
    Digits(uint64_t id,int l):L(l) {
        for(int i=0;i<L;++i,id/=3) { d[i]=id%3; nonzero+=d[i]!=0; nx+=d[i]==2; }
    }
    void next() {
        for(int i=0;i<L;++i) {
            if(d[i]==2) { d[i]=0; --nonzero; --nx; }
            else { if(d[i]==0) ++nonzero; else ++nx; ++d[i]; break; }
        }
    }
};

template<class I> void execute(const Input& in, const std::string& outpath,
                              unsigned nt, uint64_t maxbytes, const std::string& backend) {
    // Include per-object allocation and conservative allocator overhead for cpp_int.
    uint64_t element=sizeof(I);
    if(backend=="cpp_int") element+=8*((boost::multiprecision::msb(in.bound)+64)/64)+32;
    if(in.count>maxbytes/element) throw std::runtime_error("memory budget refused before allocation");
    std::cerr<<"allocate classes="<<in.count<<" principal_estimate_bytes="<<in.count*element
             <<" backend="<<backend<<" threads="<<nt<<std::endl;
    auto start=Clock::now();
    std::vector<I> a(in.count);
    std::array<uint64_t,23> two{};
    std::array<uint64_t,23> sphere{};
    for(int x=0;x<=in.L;++x) two[x]=power(2,x);
    for(int w=0;w<=in.C;++w) sphere[w]=choose(in.C,w)*power(3,w);
    for(size_t i=0;i<in.terms.size();++i) {
        const auto& t=in.terms[i];
        a[t.w*in.P+t.id]+=convert<I>(t.q);
        if((i+1)%128==0 || i+1==in.terms.size())
            std::cerr<<"populated "<<i+1<<"/"<<in.terms.size()<<std::endl;
    }
    std::vector<I> totals(nt,0);
    parallel(nt,[&](unsigned t){
        uint64_t lo=in.P*t/nt, hi=in.P*(t+1)/nt;
        Digits d(lo,in.L);
        for(uint64_t id=lo;id<hi;++id,d.next())
            for(int w=0;w<=in.C;++w)
                totals[t]+=a[w*in.P+id]*(two[d.nx]*sphere[w]);
    });
    I dense=0; for(const I& x:totals) dense+=x;
    if(dense!=convert<I>(in.M)) throw std::runtime_error("dense measure mass differs from input mass");
    std::cerr<<"independent dense mass recount PASS"<<std::endl;
    auto K=radial_matrix(in.C);
    prefix_transform(a,in.L,nt,false);
    suffix_transform(a,in,nt,false,K);
    std::vector<int64_t> lambda(in.N+1);
    for(int k=0;k<=in.N;++k) for(int v=0;v<=in.m;++v) for(int j=0;j<=v;++j) {
        int64_t z=int64_t(choose(k,j)*choose(in.N-k,v-j)*power(3,v-j));
        lambda[k]+=(j%2?-z:z);
    }
    parallel(nt,[&](unsigned t){
        uint64_t lo=in.P*t/nt, hi=in.P*(t+1)/nt;
        Digits d(lo,in.L);
        for(uint64_t id=lo;id<hi;++id,d.next())
            for(int s=0;s<=in.C;++s) a[s*in.P+id]*=lambda[d.nonzero+s];
    });
    std::cerr<<"ball eigenvalues applied"<<std::endl;
    suffix_transform(a,in,nt,true,K);
    prefix_transform(a,in.L,nt,true);
    std::vector<I> maxima(nt,-1), sums(nt,0);
    std::vector<uint64_t> args(nt,0), counts(nt,0), ambient(nt,0);
    std::atomic<bool> negative{false};
    parallel(nt,[&](unsigned t){
        uint64_t lo=in.P*t/nt, hi=in.P*(t+1)/nt;
        Digits d(lo,in.L);
        auto tick=Clock::now();
        for(uint64_t id=lo;id<hi;++id,d.next()) {
            for(int w=0;w<=in.C;++w) {
                uint64_t index=w*in.P+id;
                const I& v=a[index];
                if(v<0) negative=true;
                if(v>maxima[t] || (v==maxima[t] && index<args[t])) {maxima[t]=v;args[t]=index;}
                uint64_t multiplicity=two[d.nx]*sphere[w];
                sums[t]+=v*multiplicity;
                ambient[t]+=multiplicity;
                ++counts[t];
            }
            if(t==0 && (id-lo)%1048576==0 && Clock::now()-tick>std::chrono::seconds(10)) {
                std::cerr<<"maximum scan thread0 "<<id-lo<<"/"<<hi-lo<<std::endl;
                tick=Clock::now();
            }
        }
    });
    if(negative) throw std::runtime_error("negative final ball load");
    I maximum=-1, check=0;
    uint64_t arg=0, count=0, words=0;
    for(unsigned t=0;t<nt;++t) {
        check+=sums[t]; count+=counts[t]; words+=ambient[t];
        if(maxima[t]>maximum || (maxima[t]==maximum && args[t]<arg)) {maximum=maxima[t];arg=args[t];}
    }
    uint64_t volume=0;
    for(int d=0;d<=in.m;++d) volume+=choose(in.N,d)*power(3,d);
    if(count!=in.count || words!=power(4,in.N) || check!=convert<I>(in.M*volume))
        throw std::runtime_error("complete scan / total ball-load identity failed");
    std::string prefix;
    for(uint64_t z=arg%in.P;prefix.size()<size_t(in.L);z/=3)
        prefix.push_back(z%3==2?'X':char('0'+z%3));
    double seconds=std::chrono::duration<double>(Clock::now()-start).count();
    std::ofstream out(outpath);
    if(!out) throw std::runtime_error("cannot open output");
    out<<"{\n  \"algorithm\": \"ZX2_exact_radial_tensor_v1\",\n"
       <<"  \"backend\": \""<<backend<<"\",\n"
       <<"  \"N\": "<<in.N<<", \"r\": "<<in.r<<", \"m\": "<<in.m<<",\n"
       <<"  \"classes_checked\": "<<count<<",\n"
       <<"  \"denominator\": \""<<decimal(in.D)<<"\",\n"
       <<"  \"mass_numerator\": \""<<decimal(in.M)<<"\",\n"
       <<"  \"maximum_numerator\": \""<<decimal(maximum)<<"\",\n"
       <<"  \"argmax\": {\"prefix\": \""<<prefix<<"\", \"suffix_weight\": "<<arg/in.P<<"},\n"
       <<"  \"weighted_load_sum\": \""<<decimal(check)<<"\",\n"
       <<"  \"ambient_words\": "<<words<<",\n"
       <<"  \"principal_estimate_bytes\": "<<in.count*element<<",\n"
       <<"  \"seconds\": "<<seconds<<",\n"
       <<"  \"threads\": "<<nt;
    if(in.count<=10000) {
        out<<",\n  \"loads\": [";
        for(uint64_t i=0;i<in.count;++i) {if(i) out<<",";out<<"\""<<decimal(a[i])<<"\"";}
        out<<"]";
    }
    out<<"\n}\n";
    if(!out) throw std::runtime_error("failed to write complete output");
    std::cerr<<"COMPLETE classes="<<count<<" maximum_numerator="<<decimal(maximum)
             <<" denominator="<<decimal(in.D)<<" argmax="<<prefix<<","<<arg/in.P
             <<" seconds="<<seconds<<std::endl;
}

int main(int argc,char** argv) {
    try {
        if(argc!=6) throw std::runtime_error("usage: exact_tensor INPUT OUTPUT THREADS MAX_BYTES auto|big|wide");
        unsigned nt=std::stoul(argv[3]);
        if(nt<1 || nt>256) throw std::runtime_error("bad thread count");
        uint64_t bytes=std::stoull(argv[4]);
        Input in=read_input(argv[1]);
        bool slurm=std::getenv("SLURM_JOB_ID")!=nullptr;
#ifdef __APPLE__
        slurm=false;
#endif
        if(!slurm && (nt!=1 || bytes>1073741824ULL || in.count>2000000ULL))
            throw std::runtime_error("local guard: one native thread, <=1 GiB, <=2 million classes");
        std::string mode=argv[5];
        if(mode!="auto" && mode!="big" && mode!="wide") throw std::runtime_error("bad arithmetic mode");
        // The spare sign bits cover expression temporaries and all reductions.
        bool wide=in.bound<(Big(1)<<120);
        if(mode=="wide" && !wide) throw std::runtime_error("signed-128 bound refused");
        if(mode=="big") wide=false;
        if(wide) execute<Wide>(in,argv[2],nt,bytes,"signed128");
        else execute<Big>(in,argv[2],nt,bytes,"cpp_int");
        return 0;
    } catch(const std::exception& e) {std::cerr<<"FAIL: "<<e.what()<<std::endl;return 1;}
}
