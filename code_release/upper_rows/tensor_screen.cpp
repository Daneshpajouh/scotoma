#define TENSOR_CRT_NO_MAIN
#include "tensor_crt_ref.cpp"
#include <mutex>
int main(int argc,char**argv){
 std::ifstream in(argv[1]);std::ofstream out(argv[2]);int n,r,m,nt,packing;U len,scale,bound,threshold,limit;in>>n>>r>>m>>nt>>packing>>scale>>bound>>threshold>>limit>>len;int l=r-1,c=n-l,h=c+1;U T=1;for(int k=0;k<l;k++)T*=3;U size=T*h;std::vector<Word>words;for(U k=0;k<len;k++){U id;int w;std::string s;in>>id>>w>>s;words.push_back({id,w,parse(s)});}if(!in)return 3;
 constexpr U P1=(U(1)<<31)-1,P2=(U(1)<<31)-19,MOD=P1*P2;if(bound>=MOD)return 4;
 auto start=std::chrono::steady_clock::now();std::string tmp=std::string(argv[2])+".residue";
 {auto a=pass<1>(n,l,m,nt,size,words);std::ofstream f(tmp,std::ios::binary);f.write((char*)a.data(),size*sizeof(V));if(!f)return 5;}
 auto a=pass<19>(n,l,m,nt,size,words);std::ifstream f(tmp,std::ios::binary);U chunk=U(1)<<22;std::vector<V>buf(chunk);V inv=power<19>(18,P2-2);
 using Item=std::pair<U,U>;std::vector<std::vector<U>>flagged(nt);std::vector<U>cnt(nt),flags(nt);std::vector<std::priority_queue<Item>>lo(nt);std::vector<std::priority_queue<Item,std::vector<Item>,std::greater<Item>>>hi(nt);std::atomic<U>errors{0};
 for(U off=0;off<size;off+=chunk){U count=std::min(chunk,size-off);f.read((char*)buf.data(),count*sizeof(V));if(!f)return 6;
 par(count,nt,[&](U b,U e,int t){for(U j=b;j<e;j++){U z=off+j;U id=z/h;int w=z%h;if(!packing){U tmpid=id;int ds[20];for(int k=l-1;k>=0;k--){ds[k]=tmpid%3;tmpid/=3;}int d=0;for(int k=0;k<l;k++)d+=ds[k]!=1;int mn=d;for(int k=0;k<l;k++){d+=(ds[k]!=0)-(ds[k]!=1);mn=std::min(mn,d);}if(mn+w>m)continue;}
 int64_t d=int64_t(a[z])-buf[j];if(d<0)d+=P2;if(d<0)d+=P2;if(U(d)>=P2)d-=P2;U k=red<19>(U(d)*inv),v=P1*k+buf[j];if(v>bound){errors++;continue;}cnt[t]++;
 bool flag=packing?v>=threshold:v<=threshold;if(!flag)continue;flags[t]++;if(flagged[t].size()<limit)flagged[t].push_back(z);
 if(packing){auto&hp=hi[t];if(hp.size()<2048)hp.emplace(v,z);else if(Item{v,z}>hp.top()){hp.pop();hp.emplace(v,z);}}else{auto&hp=lo[t];if(hp.size()<2048)hp.emplace(v,z);else if(Item{v,z}<hp.top()){hp.pop();hp.emplace(v,z);}}
 }});
 }
 U nc=0,nf=0;std::vector<U>all;std::vector<Item>top;bool complete=true;for(int t=0;t<nt;t++){nc+=cnt[t];nf+=flags[t];complete&=flags[t]<=limit;all.insert(all.end(),flagged[t].begin(),flagged[t].end());while(!hi[t].empty()){top.push_back(hi[t].top());hi[t].pop();}while(!lo[t].empty()){top.push_back(lo[t].top());lo[t].pop();}}
 std::sort(all.begin(),all.end());std::sort(top.begin(),top.end());if(packing)std::reverse(top.begin(),top.end());if(top.size()>8192)top.resize(8192);
 out<<"classes "<<size<<"\nchecked "<<nc<<"\nflagged "<<nf<<"\ncomplete "<<complete<<"\nerrors "<<errors<<"\nseconds "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"\n";
 for(auto&it:top)out<<"priority "<<it.second/h<<" "<<it.second%h<<" "<<it.first<<"\n";
 if(complete)for(auto z:all)out<<"flag "<<z/h<<" "<<z%h<<"\n";
 std::remove(tmp.c_str());return errors?8:0;
}
