#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <atomic>
#include <fstream>
#include <iostream>
#include <thread>
#include <queue>
#include <vector>
using Z=boost::multiprecision::cpp_int;using U=uint64_t;
struct Term{U zero,one,x;int w;Z mass;std::vector<Z> table;};
int main(int argc,char**argv){std::ifstream in(argv[1]);std::ofstream out(argv[2]);int N,r,m,packing,nt;Z S;size_t ns,nq;in>>N>>r>>m>>packing>>nt>>S>>ns>>nq;int L=r-1,C=N-L,H=L+1,hs=H*H*(C+1);std::vector<Term>terms(ns);for(auto&t:terms)in>>t.zero>>t.one>>t.x>>t.w>>t.mass;std::vector<U>K((C+1)*hs);for(auto&v:K)in>>v;std::vector<std::pair<U,int>>queries(nq);for(auto&q:queries)in>>q.first>>q.second;if(!in)return 3;
for(auto&t:terms){t.table.resize(hs);for(int p=0;p<=C;p++)for(int d=0;d<=L;d++)for(int both=0;both<=L;both++){U q=K[t.w*hs+p*H*H+d*H+both];t.table[p*H*H+d*H+both]=t.mass*q;}}
std::atomic<U>next{0};std::vector<Z>ext(nt);std::vector<U>arg(nt),cnt(nt),bad(nt);std::vector<uint8_t>set(nt,0);using Item=std::pair<Z,U>;std::vector<std::priority_queue<Item>>lo(nt);std::vector<std::priority_queue<Item,std::vector<Item>,std::greater<Item>>>hi(nt);std::vector<std::thread>ts;
for(int ti=0;ti<nt;ti++)ts.emplace_back([&,ti]{for(;;){U start=next.fetch_add(128);if(start>=nq)break;for(U qi=start;qi<std::min(start+128,U(nq));qi++){auto[id,w]=queries[qi];U v=id,zero=0,one=0,x=0;for(int j=L-1;j>=0;j--){int d=v%3;v/=3;U bit=U(1)<<j;if(d==0)zero|=bit;else if(d==1)one|=bit;else x|=bit;}Z val;for(auto&t:terms){int forced=L-__builtin_popcountll((zero&t.zero)|(one&t.one)|(x&t.x));int both=__builtin_popcountll(x&t.x);val+=t.table[w*H*H+forced*H+both];}cnt[ti]++;U z=id*(C+1)+w;if(!set[ti]||(packing?val>ext[ti]:val<ext[ti])){set[ti]=1;ext[ti]=val;arg[ti]=z;}if(packing?val>S:val<S)bad[ti]++;
if(packing){auto&hp=hi[ti];if(hp.size()<2048)hp.emplace(val,z);else if(Item(val,z)>hp.top()){hp.pop();hp.emplace(val,z);}}else{auto&hp=lo[ti];if(hp.size()<2048)hp.emplace(val,z);else if(Item(val,z)<hp.top()){hp.pop();hp.emplace(val,z);}}
}}});for(auto&t:ts)t.join();bool first=true;Z best;U at=0,count=0,nb=0;std::vector<Item>items;for(int ti=0;ti<nt;ti++){count+=cnt[ti];nb+=bad[ti];if(set[ti]&&(first||(packing?ext[ti]>best:ext[ti]<best))){first=false;best=ext[ti];at=arg[ti];}while(!lo[ti].empty()){items.push_back(lo[ti].top());lo[ti].pop();}while(!hi[ti].empty()){items.push_back(hi[ti].top());hi[ti].pop();}}
std::sort(items.begin(),items.end());if(packing)std::reverse(items.begin(),items.end());if(items.size()>2048)items.resize(2048);out<<"scale "<<S<<"\nextremum "<<best<<"\narg "<<at/(C+1)<<" "<<at%(C+1)<<"\nchecked "<<count<<"\nviolations "<<nb<<"\n";for(auto&it:items)out<<"candidate "<<it.second/(C+1)<<" "<<it.second%(C+1)<<" "<<it.first<<"\n";}
