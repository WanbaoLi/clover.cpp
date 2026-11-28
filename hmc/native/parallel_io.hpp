#pragma once
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
