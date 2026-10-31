#include "kernel_common.hpp"
#include <fstream>
#include <iostream>
int main(int argc,char**argv) {
    if(argc!=6)return 2;
    int N=std::stoi(argv[1]),r=std::stoi(argv[2]),m=std::stoi(argv[3]);
    std::ifstream in(argv[4]);std::vector<int>ids;int id;
    while(in>>id)ids.push_back(id);
    auto T=build_kernel(N,r,m);
    std::ofstream out(argv[5],std::ios::binary);
    auto write=[&](const void*p,size_t n){out.write((const char*)p,n);};
    int hdr[5]={N,r,m,(int)ids.size(),T.n_pts};write(hdr,sizeof hdr);
    write(ids.data(),ids.size()*sizeof(int));
    write(T.point_pref.data(),T.n_pts*sizeof(int));write(T.point_w.data(),T.n_pts*sizeof(int));
    write(T.rhs.data(),T.n_pts*sizeof(int64_t));
    for(auto i:ids)write(T.scaled.data()+size_t(i)*T.n_pts,T.n_pts*sizeof(int64_t));
    std::cout<<"matrix "<<ids.size()<<' '<<T.n_pts<<'\n';
    return out?0:3;
}
