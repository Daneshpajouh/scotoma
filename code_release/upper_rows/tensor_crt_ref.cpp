#define main unused_tensor_main
#include "tensor_exact_ref.cpp"
#undef main
#include <cstdio>
using V=uint32_t;
struct Info{int d0=0,end=0,mn=0,cnt=1,arg=0,nx=0;};
Info info(U id,int len){Info q;std::vector<int>d(len);for(int k=len-1;k>=0;k--){d[k]=id%3;id/=3;}for(int k=0;k<len;k++){int a=d[k];q.d0+=a!=1;q.nx+=a==2;q.end+=(a==1)-(a==0);if(q.end<q.mn){q.mn=q.end;q.cnt=1;q.arg=k+1;}else if(q.end==q.mn)q.cnt++;}return q;}
int priv(const Info&a,const Info&b){int t=a.end+b.mn,mn=std::min(a.mn,t),cnt=(a.mn==mn?a.cnt:0)+(t==mn?b.cnt:0)-(a.end==mn);if(a.d0+b.d0+mn!=9||cnt!=1)return -1;return a.mn<t?a.arg:t<a.mn?9+b.arg:9;}

template<U C> inline V red(U x){constexpr U M=(U(1)<<31)-1,P=(U(1)<<31)-C;x=(x&M)+(x>>31)*C;x=(x&M)+(x>>31)*C;if(x>=P)x-=P;if(x>=P)x-=P;return V(x);}
template<U C> V power(V a,U b){V v=1;while(b){if(b&1)v=red<C>(U(v)*a);a=red<C>(U(a)*a);b>>=1;}return v;}
template<U C> std::vector<V> pass(int n,int l,int m,int nt,U size,const std::vector<Word>&words){
 constexpr U P=(U(1)<<31)-C;int c=n-l,h=c+1;U T=size/h;std::vector<V>a(size,0);for(auto w:words)a[w.id*h+w.w]=red<C>(a[w.id*h+w.w]+U(w.a%P));
 std::vector<V>K(h*h);for(int t=0;t<h;t++)for(int w=0;w<h;w++){I v=kp(c,t,w)%I(P);if(v<0)v+=P;K[t*h+w]=V(v);}
 auto tail=[&]{if(!c)return;par(T,nt,[&](U b,U e,int){for(U z=b;z<e;z++){V v[21]={};for(int t=0;t<h;t++)for(int w=0;w<h;w++)v[t]=red<C>(v[t]+U(red<C>(U(K[t*h+w])*a[z*h+w])));for(int t=0;t<h;t++)a[z*h+t]=v[t];}});};
 auto tr=[&]{for(U stride=h;stride<size;stride*=3){par(size/3,nt,[&](U b,U e,int){for(U i=b;i<e;i++){U g=i/stride,k=i%stride,base=g*stride*3;U x=a[base+k],y=a[base+stride+k],v=a[base+2*stride+k];a[base+k]=red<C>(x+y+2*v);a[base+stride+k]=red<C>(x+y+2*P-2*v);a[base+2*stride+k]=red<C>(x+P-y);}});}};
 tail();tr();std::vector<V>lambda(n+1);for(int s=0;s<=n;s++){I z=0;for(int i=0;i<=m;i++)z+=kp(n,s,i);z%=P;if(z<0)z+=P;lambda[s]=V(z);}
 par(T,nt,[&](U b,U e,int){for(U z=b;z<e;z++){U v=z;int s=0;for(int k=0;k<l;k++){s+=v%3!=0;v/=3;}for(int t=0;t<h;t++)a[z*h+t]=red<C>(U(a[z*h+t])*lambda[s+t]);}});
 tr();tail();V inv=power<C>(power<C>(4,n),P-2);par(size,nt,[&](U b,U e,int){for(U z=b;z<e;z++)a[z]=red<C>(U(a[z])*inv);});return a;
}
#ifndef TENSOR_CRT_NO_MAIN
int main(int argc,char**argv){
 if(argc<3)return 2;std::ifstream in(argv[1]);std::ofstream out(argv[2]);int n,r,m,nt;std::string ds,qs;U len;in>>n>>r>>m>>ds>>qs>>len>>nt;I D=parse(ds),Q=parse(qs);int l=r-1,c=n-l,h=c+1;U T=1;for(int k=0;k<l;k++)T*=3;U size=T*h;std::vector<Word>words;for(U k=0;k<len;k++){U id;int w;std::string s;in>>id>>w>>s;words.push_back({id,w,parse(s)});}if(!in)return 3;
 constexpr U P1=(U(1)<<31)-1,P2=(U(1)<<31)-19,MOD=P1*P2;I denI=D*Q,B=r*denI;if(argc>3&&std::string(argv[3])=="price"){B=0;for(auto w:words){U id=w.id;int nx=0;for(int k=0;k<l;k++){nx+=id%3==2;id/=3;}B+=w.a*(I(1)<<nx);}}if(B>=I(MOD)){std::cerr<<"CRT modulus too small\n";return 4;}U den=U(denI),bound=U(B);
 auto start=std::chrono::steady_clock::now();std::string tmp=std::string(argv[2])+".residue";
 {auto a=pass<1>(n,l,m,nt,size,words);std::ofstream f(tmp,std::ios::binary);f.write((char*)a.data(),size*sizeof(V));if(!f)return 5;}
 std::cerr<<"first prime seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 auto a=pass<19>(n,l,m,nt,size,words);std::cerr<<"second prime seconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 std::ifstream f(tmp,std::ios::binary);std::ifstream validation;bool pricing=argc>3&&std::string(argv[3])=="price";bool validate=argc>3&&!pricing;if(validate)validation.open(argv[3],std::ios::binary);std::vector<U>pats;for(int j=0;j<=l;j++){U id=0;for(int k=0;k<l;k++)id=id*3+(k>=j);pats.push_back(id*h);}std::sort(pats.begin(),pats.end());
 std::vector<Info>left,right;if(pricing){if(n!=19||r!=20||m!=9)return 10;left.resize(19683);right.resize(59049);for(U z=0;z<left.size();z++)left[z]=info(z,9);for(U z=0;z<right.size();z++)right[z]=info(z,10);}
 using Item=std::pair<U,U>;std::vector<std::vector<Item>>allbad(nt);std::vector<std::priority_queue<Item>>lowest(nt*20);std::vector<U>pcnt(nt*20,0),card(nt*20,0);
 std::vector<std::priority_queue<Item,std::vector<Item>,std::greater<Item>>>top(nt);std::vector<U>best(nt,0),arg(nt,0),bad(nt,0);std::atomic<U>errors{0},checks{0};U chunk=U(1)<<23;std::vector<V>buf(chunk);std::vector<I>old;if(validate)old.resize(chunk);V inv=power<19>(18,P2-2);
 for(U off=0;off<size;off+=chunk){U count=std::min(chunk,size-off);f.read((char*)buf.data(),count*sizeof(V));if(!f)return 6;if(validate){validation.read((char*)old.data(),count*sizeof(I));if(!validation)return 7;}
  par(count,nt,[&](U b,U e,int t){auto pp=std::lower_bound(pats.begin(),pats.end(),off+b);U next=pp==pats.end()?size:*pp;for(U j=b;j<e;j++){U z=off+j;int region=-1;if(pricing){region=priv(left[z/59049],right[z%59049]);if(region<0)continue;}int64_t d=int64_t(a[z])-buf[j];if(d<0)d+=P2;if(d<0)d+=P2;if(U(d)>=P2)d-=P2;U k=red<19>(U(d)*inv),v=P1*k+buf[j];if(v>bound)errors++;if(validate){if(I(v)!=old[j])errors++;checks++;}if(pricing){int ix=t*20+region;pcnt[ix]++;card[ix]+=U(1)<<(left[z/59049].nx+right[z%59049].nx);auto&hp=lowest[ix];if(hp.size()<64)hp.emplace(v,z);else if(Item{v,z}<hp.top()){hp.pop();hp.emplace(v,z);}continue;}if(z==next){if(v!=den)errors++;++pp;next=pp==pats.end()?size:*pp;continue;}if(v>best[t]){best[t]=v;arg[t]=z;}if(v>den){bad[t]++;if(allbad[t].size()<200000)allbad[t].emplace_back(v,z);auto&hp=top[t];if(hp.size()<128)hp.emplace(v,z);else if(Item{v,z}>hp.top()){hp.pop();hp.emplace(v,z);}}}});
 }
 U v=0,z=0,badtotal=0;for(int t=0;t<nt;t++){badtotal+=bad[t];if(best[t]>v||(best[t]==v&&arg[t]<z)){v=best[t];z=arg[t];}}
 out<<"classes "<<size<<"\nnonpattern_max "<<frac(v,den)<<"\narg "<<z/h<<" "<<z%h<<"\nviolations "<<badtotal<<"\nerrors "<<errors<<"\ncrt_checks "<<checks<<"\nseconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 std::vector<Item>viol;for(auto&vv:allbad)viol.insert(viol.end(),vv.begin(),vv.end());for(auto&hp:top)while(!hp.empty()){viol.push_back(hp.top());hp.pop();}std::sort(viol.rbegin(),viol.rend());viol.erase(std::unique(viol.begin(),viol.end()),viol.end());if(viol.size()>100000)viol.resize(100000);for(auto[w,zz]:viol)out<<"badcenter "<<zz/h<<" "<<zz%h<<" "<<frac(w,den)<<"\n";
 if(pricing){U total=0;for(int j=0;j<20;j++){std::vector<Item>vv;U cnt=0,sz=0;for(int t=0;t<nt;t++){int ix=t*20+j;cnt+=pcnt[ix];sz+=card[ix];while(!lowest[ix].empty()){vv.push_back(lowest[ix].top());lowest[ix].pop();}}std::sort(vv.begin(),vv.end());if(vv.empty())return 11;total+=cnt;out<<"private_region "<<j<<" "<<cnt<<" "<<sz<<"\nprivate_min "<<j<<" "<<frac(vv[0].first,den)<<"\n";if(vv.size()>100)vv.resize(100);for(auto[w,zz]:vv)out<<"private_column "<<j<<" "<<zz<<" "<<frac(w,den)<<"\n";}out<<"private_orbits "<<total<<"\n";}
 std::remove(tmp.c_str());return errors?8:0;
}

#endif
