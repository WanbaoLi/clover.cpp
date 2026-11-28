#pragma once
#include "runtime/bridge_setup.h"
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
#if defined(USE_MPI)
#include <mpi.h>
#endif
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>

namespace hmc_mpi {
static int rank=0, ranks=1;
static std::array<int,4> global{}, local{}, offset{};

#if defined(USE_MPI)
inline void check(int code,const char* what) {
  if (code!=MPI_SUCCESS) {
    char message[MPI_MAX_ERROR_STRING]; int len=0;
    MPI_Error_string(code,message,&len);
    throw std::runtime_error(std::string(what)+": "+std::string(message,len));
  }
}

// The controller and native ranks share a regular file. Each rank owns a
// disjoint hyperrectangle in its [t,z,y,x,component] layout. Join adjacent
// rows into a single I/O operation wherever the decomposition permits it.
inline void transfer(const char* path,double* buffer,int components,bool writing) {
  static_assert(sizeof(off_t)>=8,"HMC requires large-file offsets");
  std::uint64_t global_sites=1;
  std::size_t local_sites=1,contiguous=1;
  bool adjoining=true;
  for(int d=0;d<4;++d) {
    global_sites*=global[d];
    local_sites*=local[d];
    if(adjoining) contiguous*=local[d];
    adjoining=adjoining && local[d]==global[d];
  }
  const off_t size=static_cast<off_t>(global_sites*components*sizeof(double));
  int error=0;
  if(writing) {
    if(rank==0) {
      const int fd=::open(path,O_WRONLY|O_CREAT|O_TRUNC,0600);
      if(fd<0) error=errno;
      else {
        if(::ftruncate(fd,size)!=0) error=errno;
        if(::close(fd)!=0 && !error) error=errno;
      }
    }
    check(MPI_Bcast(&error,1,MPI_INT,0,MPI_COMM_WORLD),"prepare output field");
    if(error) throw std::runtime_error(std::string(path)+": "+std::strerror(error));
  }
  const int fd=::open(path,writing ? O_WRONLY : O_RDONLY);
  if(fd<0) error=errno;
  else {
    struct stat st;
    if(::fstat(fd,&st)!=0) error=errno;
    else if(st.st_size!=size) error=EINVAL;
    for(std::size_t first=0;first<local_sites && !error;first+=contiguous) {
      std::size_t site=first;
      std::uint64_t global_site=0,stride=1;
      for(int d=0;d<4;++d) {
        global_site+=(site%local[d]+offset[d])*stride;
        site/=local[d];
        stride*=global[d];
      }
      off_t position=static_cast<off_t>(global_site*components*sizeof(double));
      char* data=reinterpret_cast<char*>(buffer+first*components);
      std::size_t remaining=contiguous*components*sizeof(double);
      while(remaining) {
        const ssize_t n=writing ? ::pwrite(fd,data,remaining,position) : ::pread(fd,data,remaining,position);
        if(n<0 && errno==EINTR) continue;
        if(n<=0) {error=n<0 ? errno : EIO;break;}
        position+=n;data+=n;remaining-=n;
      }
    }
    if(::close(fd)!=0 && !error) error=errno;
  }
  check(MPI_Allreduce(MPI_IN_PLACE,&error,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD),"complete field access");
  if(error) throw std::runtime_error(std::string(path)+": "+std::strerror(error));
}
#else
inline void transfer(const char*,double*,int,bool) {
  throw std::runtime_error("distributed field I/O requires an MPI build");
}
#endif

inline void parallel_size() {
#if defined(USE_MPI)
  check(MPI_Comm_rank(MPI_COMM_WORLD,&rank),"worker rank");
  check(MPI_Comm_size(MPI_COMM_WORLD,&ranks),"worker size");
#endif
}

inline void finish_request() {
#if defined(USE_MPI)
  if(ranks>1) check(MPI_Barrier(MPI_COMM_WORLD),"worker completion");
#endif
}

inline long double sum_determinant(long double local_value) {
#if defined(USE_MPI)
  if(ranks>1) check(MPI_Allreduce(MPI_IN_PLACE,&local_value,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD),"sum local determinant");
#endif
  return local_value;
}

inline void report(const char* format,...) {
  if (rank!=0) return;
  va_list args;
  va_start(args,format);
  std::vprintf(format,args);
  va_end(args);
}

inline bool request(int& mode,double& num,double& den,double& tolerance,int& limit) {
  if (ranks==1) return bool(std::cin>>mode>>num>>den>>tolerance>>limit) && mode>=0;
#if defined(USE_MPI)
  std::array<int,3> header{};
  std::array<double,3> params{};
  if (rank==0 && (std::cin>>mode>>num>>den>>tolerance>>limit) && mode>=0) {
    header={{1,mode,limit}}; params={{num,den,tolerance}};
  }
  // Other persistent workers occupy the same allocation and should sleep while idle.
  MPI_Request req;
  check(MPI_Ibcast(header.data(),3,MPI_INT,0,MPI_COMM_WORLD,&req),"worker command");
  int ready=0;
  while (!ready) {
    check(MPI_Test(&req,&ready,MPI_STATUS_IGNORE),"worker command wait");
    if (!ready) std::this_thread::sleep_for(std::chrono::milliseconds(5));
  }
  if (!header[0]) return false;
  check(MPI_Bcast(params.data(),3,MPI_DOUBLE,0,MPI_COMM_WORLD),"worker parameters");
  mode=header[1]; limit=header[2]; num=params[0]; den=params[1]; tolerance=params[2];
  return true;
#else
  return false;
#endif
}
}


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
