#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <fstream>
#include <iostream>
#include <queue>
#include <thread>
#include <vector>
using Z=boost::multiprecision::cpp_int;
using U=uint64_t;
long long choose(int n,int k){if(k<0||k>n)return 0;long long z=1;for(int j=1;j<=k;j++)z=z*(n-j+1)/j;return z;}
long long kp(int n,int t,int w){long long z=0;for(int a=0;a<=w;a++){long long p=1;for(int k=0;k<w-a;k++)p*=3;z+=(a%2?-1:1)*choose(t,a)*choose(n-t,w-a)*p;}return z;}
template<class F> void par(U count,int nt,F f){std::atomic<U> next{0};std::vector<std::thread> ts;for(int t=0;t<nt;t++)ts.emplace_back([&,t]{for(;;){U b=next.fetch_add(1024);if(b>=count)break;f(b,std::min(b+1024,count),t);}});for(auto&t:ts)t.join();}
int main(int argc,char**argv){
 if(argc<3)return 2;std::ifstream in(argv[1]);std::ofstream out(argv[2]);int n,r,m,nt,packing,keep;Z scale;U len;in>>n>>r>>m>>nt>>packing>>keep>>scale>>len;int l=r-1,c=n-l,h=c+1;U T=1;for(int k=0;k<l;k++)T*=3;U size=T*h;std::vector<Z>a(size);
 for(U k=0;k<len;k++){U id;int w;Z v;in>>id>>w>>v;if(id>=T||w<0||w>c||v<0)return 3;a[id*h+w]+=v;}if(!in)return 4;
 std::vector<long long>K(h*h);for(int t=0;t<h;t++)for(int w=0;w<h;w++)K[t*h+w]=kp(c,t,w);
 auto start=std::chrono::steady_clock::now();
 auto tail=[&]{if(!c)return;par(T,nt,[&](U b,U e,int){for(U z=b;z<e;z++){Z v[21];for(int t=0;t<h;t++)for(int w=0;w<h;w++)v[t]+=K[t*h+w]*a[z*h+w];for(int t=0;t<h;t++)a[z*h+t]=std::move(v[t]);}});};
 auto tr=[&]{for(U stride=h;stride<size;stride*=3){par(size/3,nt,[&](U b,U e,int){for(U i=b;i<e;i++){U g=i/stride,k=i%stride,base=g*stride*3;Z x=a[base+k],y=a[base+stride+k],v=a[base+2*stride+k];a[base+k]=x+y+2*v;a[base+stride+k]=x+y-2*v;a[base+2*stride+k]=x-y;}});}};
 tail();tr();std::vector<long long>lambda(n+1);for(int s=0;s<=n;s++)for(int i=0;i<=m;i++)lambda[s]+=kp(n,s,i);
 par(T,nt,[&](U b,U e,int){for(U z=b;z<e;z++){U v=z;int s=0;for(int k=0;k<l;k++){s+=v%3!=0;v/=3;}for(int t=0;t<h;t++)a[z*h+t]*=lambda[s+t];}});
 tr();tail();Z div=Z(1)<<(2*n);std::vector<Z>best(nt);std::vector<U>arg(nt),bad(nt),checks(nt);std::vector<uint8_t>set(nt,0);using Item=std::pair<Z,U>;
 std::vector<std::priority_queue<Item>> low(nt);std::vector<std::priority_queue<Item,std::vector<Item>,std::greater<Item>>> high(nt);std::atomic<U>errors{0};
 par(T,nt,[&](U b,U e,int t){for(U id=b;id<e;id++){U tmp=id;int digit[20];for(int k=l-1;k>=0;k--){digit[k]=tmp%3;tmp/=3;}int dist=0;for(int k=0;k<l;k++)dist+=digit[k]!=1;int dmin=dist;for(int k=0;k<l;k++){dist+=(digit[k]!=0)-(digit[k]!=1);dmin=std::min(dmin,dist);}for(int w=0;w<h;w++){U z=id*h+w;if(!packing&&dmin+w>m)continue;if(a[z]%div!=0){errors++;continue;}Z v=a[z]/div;checks[t]++;if(!set[t]||(packing?v>best[t]:v<best[t])){best[t]=v;arg[t]=z;set[t]=true;}bool violation=packing?v>scale:v<scale;if(violation)bad[t]++;if(keep){if(packing){auto&hp=high[t];if(hp.size()<U(keep))hp.emplace(v,z);else if(Item(v,z)>hp.top()){hp.pop();hp.emplace(v,z);}}else{auto&hp=low[t];if(hp.size()<U(keep))hp.emplace(v,z);else if(Item(v,z)<hp.top()){hp.pop();hp.emplace(v,z);}}}}}});
 bool first=true;Z ext;U at=0,nb=0,nc=0;std::vector<Item>items;for(int t=0;t<nt;t++){nb+=bad[t];nc+=checks[t];if(set[t]&&(first||(packing?best[t]>ext:best[t]<ext))){ext=best[t];at=arg[t];first=false;}while(!high[t].empty()){items.push_back(high[t].top());high[t].pop();}while(!low[t].empty()){items.push_back(low[t].top());low[t].pop();}}
 std::sort(items.begin(),items.end());if(packing)std::reverse(items.begin(),items.end());if(items.size()>U(keep))items.resize(keep);
 out<<"scale "<<scale<<"\nextremum "<<ext<<"\narg "<<at/h<<" "<<at%h<<"\nclasses "<<size<<"\nchecked "<<nc<<"\nviolations "<<nb<<"\nerrors "<<errors<<"\nseconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";if(argc>3){bool eligible=std::string(argv[3])=="eligible";Z threshold=argc>4?Z(argv[4]):Z(0);for(U id=0;id<T;id++)for(int w=0;w<h;w++){Z v=a[id*h+w]/div;if(!eligible||v>=threshold)out<<(eligible?"eligible ":"value ")<<id<<" "<<w<<" "<<v<<"\n";}}for(auto&it:items)out<<"candidate "<<it.second/h<<" "<<it.second%h<<" "<<it.first<<"\n";return errors?5:0;
}
