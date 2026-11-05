#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <tuple>
#include <string>
#include <thread>
#include <vector>
using I=__int128_t; using U=uint64_t;
std::string str(I x){if(!x)return "0";bool neg=x<0;if(neg)x=-x;std::string s;while(x){s.push_back('0'+x%10);x/=10;}if(neg)s.push_back('-');std::reverse(s.begin(),s.end());return s;}
I parse(std::string s){I v=0;bool neg=s[0]=='-';for(size_t k=neg;k<s.size();++k)v=v*10+s[k]-'0';return neg?-v:v;}
I gcd(I a,I b){while(b){I r=a%b;a=b;b=r;}return a;}
std::string frac(I a,I b){I d=gcd(a,b);return str(a/d)+"/"+str(b/d);}
long long choose(int n,int k){if(k<0||k>n)return 0;long long z=1;for(int j=1;j<=k;j++)z=z*(n-j+1)/j;return z;}
long long kp(int n,int t,int w){long long z=0;for(int a=0;a<=w;a++){long long p=1;for(int k=0;k<w-a;k++)p*=3;z+=(a%2?-1:1)*choose(t,a)*choose(n-t,w-a)*p;}return z;}
template<class F> void par(U count,int nt,F f){std::atomic<U> next{0};std::vector<std::thread> ts;for(int t=0;t<nt;t++)ts.emplace_back([&,t]{for(;;){U a=next.fetch_add(256);if(a>=count)break;f(a,std::min(a+256,count),t);}});for(auto &t:ts)t.join();}
struct Word{U id;int w;I a;};
int main(int argc,char**argv){
 if(argc<3)return 2;std::ifstream in(argv[1]);std::ofstream out(argv[2]);int n,r,m,nt;std::string ds,qs;U len;in>>n>>r>>m>>ds>>qs>>len>>nt;I D=parse(ds),Q=parse(qs);int l=r-1,c=n-l,h=c+1;U T=1;for(int k=0;k<l;k++)T*=3;U size=T*h;std::vector<I>a(size,0);std::vector<Word>words;
 for(U k=0;k<len;k++){U id;int w;std::string s;in>>id>>w>>s;I v=parse(s);a[id*h+w]+=v;words.push_back({id,w,v});}if(!in){std::cerr<<"bad input\n";return 3;}
 std::vector<long long>K(h*h);for(int t=0;t<h;t++)for(int w=0;w<h;w++)K[t*h+w]=kp(c,t,w);
 auto start=std::chrono::steady_clock::now();
 auto tail=[&]{if(!c)return;par(T,nt,[&](U b,U e,int){for(U z=b;z<e;z++){I v[21]={};for(int t=0;t<h;t++)for(int w=0;w<h;w++)v[t]+=K[t*h+w]*a[z*h+w];for(int t=0;t<h;t++)a[z*h+t]=v[t];}});};
 auto tr=[&]{for(U stride=h;stride<size;stride*=3){par(size/3,nt,[&](U b,U e,int){for(U i=b;i<e;i++){U g=i/stride,k=i%stride,base=g*stride*3;I x=a[base+k],y=a[base+stride+k],v=a[base+2*stride+k];a[base+k]=x+y+2*v;a[base+stride+k]=x+y-2*v;a[base+2*stride+k]=x-y;}});}};
 tail();tr();std::cerr<<"forward "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" seconds\n";
 std::vector<long long>lambda(n+1);for(int s=0;s<=n;s++)for(int i=0;i<=m;i++)lambda[s]+=kp(n,s,i);
 par(T,nt,[&](U b,U e,int){for(U z=b;z<e;z++){U v=z;int s=0;for(int k=0;k<l;k++){s+=v%3!=0;v/=3;}for(int t=0;t<h;t++)a[z*h+t]*=lambda[s+t];}});
 tr();tail();std::cerr<<"inverse "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" seconds\n";
 I div=I(1)<<(2*n),den=D*Q;std::vector<I>best(nt,-1);std::vector<U>arg(nt,0),bad(nt,0);using Item=std::pair<I,U>;std::vector<std::priority_queue<Item,std::vector<Item>,std::greater<Item>>>top(nt);std::atomic<int>err{0};std::vector<U>pats;for(int j=0;j<=l;j++){U id=0;for(int k=0;k<l;k++)id=id*3+(k>=j);pats.push_back(id*h);}
 par(size,nt,[&](U b,U e,int t){for(U z=b;z<e;z++){I v=a[z];if(v%div){err++;continue;}v/=div;a[z]=v;if(std::find(pats.begin(),pats.end(),z)!=pats.end()){if(v!=den)err++;continue;}if(v>den){bad[t]++;auto &hp=top[t];if(hp.size()<128)hp.emplace(v,z);else if(Item{v,z}>hp.top()){hp.pop();hp.emplace(v,z);}}if(v>best[t]){best[t]=v;arg[t]=z;}}});
 if(argc>3&&std::string(argv[3])=="direct"){
  std::atomic<U>checks{0};par(size,nt,[&](U b,U e,int){for(U z=b;z<e;z++){U co=z/h;int p=z%h;std::vector<int>ct(l);for(int k=l-1;k>=0;k--){ct[k]=co%3;co/=3;}I val=0;for(auto w:words){U v=w.id;int aa=0,bb=0;for(int k=l-1;k>=0;k--){int d=v%3;v/=3;if(d==2&&ct[k]==2)aa++;else if(d!=ct[k])bb++;}I num=0;for(int k=0;k<=aa;k++)for(int ov=0;ov<=p&&ov<=w.w;ov++)for(int eq=0;eq<=ov;eq++){if(bb+k+p+w.w-ov-eq>m)continue;long long q=choose(aa,k)*choose(p,ov)*choose(c-p,w.w-ov)*choose(ov,eq);for(int t=0;t<ov-eq;t++)q*=2;for(int t=0;t<w.w-ov;t++)q*=3;num+=q;}int nx=0;v=w.id;for(int k=0;k<l;k++){nx+=v%3==2;v/=3;}val+=w.a*num*(I(1)<<(nx-aa));}if(val!=a[z]){if(err.fetch_add(1)<3)std::cerr<<"DIRECT MISMATCH "<<z<<" "<<str(val)<<" "<<str(a[z])<<"\n";}checks++;}});out<<"direct_checks "<<checks<<"\n";
 }
 I v=-1;U z=0,badtotal=0;for(int t=0;t<nt;t++){badtotal+=bad[t];if(best[t]>v||(best[t]==v&&arg[t]<z)){v=best[t];z=arg[t];}}
 out<<"classes "<<size<<"\nnonpattern_max "<<frac(v,den)<<"\narg "<<z/h<<" "<<z%h<<"\nviolations "<<badtotal<<"\nerrors "<<err<<"\nseconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 std::vector<Item>viol;for(auto &hp:top)while(!hp.empty()){viol.push_back(hp.top());hp.pop();}std::sort(viol.rbegin(),viol.rend());if(viol.size()>1024)viol.resize(1024);for(auto [w,zz]:viol)out<<"badcenter "<<zz/h<<" "<<zz%h<<" "<<frac(w,den)<<"\n";
 if(argc>3&&std::string(argv[3])=="dump"){std::ofstream f(std::string(argv[2])+".values",std::ios::binary);f.write((char*)a.data(),size*sizeof(I));}
 if(err)return 4;return 0;
}
