#pragma once
#include "lib/bridge_setup.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include "parallel_io.hpp"

using Matrix=std::array<std::complex<double>,9>;
static int volume=0,half_volume=0,global_volume=0,max_iterations=0;
static double target=0.0;
static std::vector<int> dimensions;
static std::vector<Matrix> gauge_store;
constexpr double light_mass=-0.279,coeff=1.160678675188,rho=0.125;

inline void initialize(int& argc,char**& argv,const char* thread_setting=nullptr) {
  if(argc!=5 && argc!=9) throw std::runtime_error("four lattice extents and an optional four-dimensional MPI grid required");
  std::vector<int> grid(4,1); int expected=1;
  int sites=1;
  for(int mu=0;mu<4;++mu) {
    const int n=std::stoi(argv[mu+1]);
    if(n<2 || n%2 || sites>std::numeric_limits<int>::max()/72/n)
      throw std::runtime_error("unsupported lattice extent or site count");
    if(argc==9) grid[mu]=std::stoi(argv[mu+5]);
    if(grid[mu]<1 || n%grid[mu] || (n/grid[mu])%2)
      throw std::runtime_error("HMC grid must give integral even local extents");
    expected*=grid[mu]; dimensions.push_back(n);sites*=n;
    hmc_mpi::global[mu]=n;
  }
  global_volume=sites;
  bridge_initialize(&argc,&argv);
  int threads=1;
  if(thread_setting) {
    if(const char* value=std::getenv(thread_setting)) {
      std::size_t consumed=0;
      threads=std::stoi(value,&consumed);
      if(threads<1 || value[consumed]!='\0')
        throw std::runtime_error(std::string(thread_setting)+" requires a positive integer");
    }
  }
  hmc_mpi::parallel_size();
  if(hmc_mpi::ranks!=expected) throw std::runtime_error("MPI size differs from the requested HMC grid");
  bridge_setup(dimensions,grid,threads,3,"stdout","Crucial");
  volume=CommonParameters::Nvol();half_volume=volume/2;gauge_store.resize(4*volume);
  for(int mu=0;mu<4;++mu) {
    dimensions[mu]=hmc_mpi::global[mu]/grid[mu];
    hmc_mpi::local[mu]=dimensions[mu];
    hmc_mpi::offset[mu]=Communicator::ipe(mu)*dimensions[mu];
  }
}
inline int parity(int site) {
  int sum=0;
  for(const int n:dimensions) {sum+=site%n; site/=n;}
  return sum%2;
}
inline Matrix eye() {Matrix m{};m[0]=m[4]=m[8]=1.;return m;}
inline void read_binary_field(std::vector<Matrix>& field,const char* path) {
  if(hmc_mpi::ranks>1) {
    hmc_mpi::transfer(path,reinterpret_cast<double*>(field.data()),72,false);
    return;
  }
  std::ifstream in(path,std::ios::binary);
  in.read(reinterpret_cast<char*>(field.data()),field.size()*sizeof(Matrix));
  if(!in || in.peek()!=EOF) throw std::runtime_error(std::string("invalid lattice field ")+path);
}
inline Matrix link_matrix(int site,int mu,double) {return gauge_store[4*site+mu];}
inline void write_binary_field(const std::vector<Matrix>& field,const char* file="moved.bin") {
  if(hmc_mpi::ranks>1) {
    hmc_mpi::transfer(file,const_cast<double*>(reinterpret_cast<const double*>(field.data())),72,true);
    return;
  }
  std::ofstream out(file,std::ios::binary);
  out.write(reinterpret_cast<const char*>(field.data()),field.size()*sizeof(Matrix));
  if(!out) throw std::runtime_error(std::string("cannot write ")+file);
}
inline void save_force(const Field_G& f,const std::string& name) {
  std::vector<Matrix> links(4*volume);
  for(int i=0;i<volume;++i) for(int mu=0;mu<4;++mu) for(int k=0;k<9;++k)
    links[4*i+mu][k]={f.cmp_r(k,i,mu),f.cmp_i(k,i,mu)};
  write_binary_field(links,name.c_str());
}
inline void done() {
  hmc_mpi::finish_request();
  if(hmc_mpi::rank==0) {std::fflush(stdout);std::cout<<"REQUEST_DONE"<<std::endl;}
}
