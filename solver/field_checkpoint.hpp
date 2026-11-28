#pragma once
#include <array>
#include <cstdint>
#include <cstring>
#include <string>

// FIPS 180-4 SHA-256 for an already resident native rank-local field.
namespace clover_checkpoint {
inline uint32_t rotate_right(uint32_t x, unsigned n) { return (x >> n) | (x << (32-n)); }

inline void hash_block(uint32_t h[8], const unsigned char block[64]) {
  static const uint32_t k[64]={
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
  uint32_t w[64];
  for (unsigned i=0;i<16;++i) {
    const unsigned char* p=block+4*i;
    w[i]=(uint32_t(p[0])<<24)|(uint32_t(p[1])<<16)|(uint32_t(p[2])<<8)|p[3];
  }
  for (unsigned i=16;i<64;++i) {
    uint32_t u=w[i-15],v=w[i-2];
    w[i]=w[i-16]+(rotate_right(u,7)^rotate_right(u,18)^(u>>3))+w[i-7]+(rotate_right(v,17)^rotate_right(v,19)^(v>>10));
  }
  uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],z=h[7];
  for (unsigned i=0;i<64;++i) {
    uint32_t t1=z+(rotate_right(e,6)^rotate_right(e,11)^rotate_right(e,25))+((e&f)^(~e&g))+k[i]+w[i];
    uint32_t t2=(rotate_right(a,2)^rotate_right(a,13)^rotate_right(a,22))+((a&b)^(a&c)^(b&c));
    z=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
  }
  h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=z;
}

inline std::string sha256(const void* data, std::size_t size) {
  uint32_t h[8]={0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
  const unsigned char* bytes=static_cast<const unsigned char*>(data);
  std::size_t full=size-size%64;
  for (std::size_t n=0;n<full;n+=64) hash_block(h,bytes+n);
  unsigned char tail[128]={};const std::size_t rest=size-full, used=rest<56?64:128;
  if(rest) std::memcpy(tail,bytes+full,rest);
  tail[rest]=0x80;
  for(unsigned n=0;n<8;++n) tail[used-1-n]=static_cast<unsigned char>((uint64_t(size)*8)>>(8*n));
  hash_block(h,tail);if(used==128) hash_block(h,tail+64);
  const char hex[]="0123456789abcdef";std::string result;
  for(uint32_t word:h) for(int shift=28;shift>=0;shift-=4) result.push_back(hex[(word>>shift)&15]);
  return result;
}
}

#include "solver_cli.hpp"
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <locale>
#include <mpi.h>

namespace clover_checkpoint {
// Every rank has agreed to this failure through an MPI reduction.
class CollectiveError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

inline std::string signature(const SolveOptions& o,int rhs) {
  std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(17);
  s<<"clover-bicgstab-checkpoint-v1\n"
   <<"operator\tWilson-Clover,Dirac,Stout1-4D,bc=1,1,1,-1,unscaled-native-FP64\n"
   <<"lattice\t"<<solve_join(o.lattice)<<"\ngrid\t"<<solve_join(o.grid)<<"\nsource\t"<<o.source
   <<"\nposition\t"<<solve_join(o.position)<<"\nrhs\t"<<rhs<<"\nmass\t"<<o.mass<<"\ncsw\t"<<o.csw<<"\nrho\t"<<o.rho<<'\n';
  return s.str();
}

inline std::string prefix(const std::string& dir,int rhs) { return dir+"/checkpoint-rhs-"+std::to_string(rhs); }

inline bool all_ok(bool local) {
  int value=local?1:0,ok=0;
  MPI_Allreduce(&value,&ok,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  return ok!=0;
}

inline void write(const SolveOptions& o,int rhs,int rank,int ranks,const std::string& coordinates,
                  const std::string& gauge_hash,const double* data,std::size_t count,double residual) {
  const std::string base=prefix(o.output,rhs),shard=base+"-"+coordinates+".f64";
  std::ofstream f(shard,std::ios::binary);const std::size_t bytes=count*sizeof(double);
  f.write(reinterpret_cast<const char*>(data),bytes);f.close();
  if(!all_ok(bool(f))) throw CollectiveError("cannot write all checkpoint shards");
  const std::string field_hash=sha256(data,bytes);
  std::array<char,256> record{};
  std::string row=std::to_string(rank)+"\t"+coordinates+"\t"+gauge_hash+"\t"+field_hash;
  if(row.size()>=record.size()) throw std::runtime_error("checkpoint rank record too long");
  std::memcpy(record.data(),row.data(),row.size());
  std::vector<char> records(rank==0?std::size_t(ranks)*record.size():0);
  MPI_Gather(record.data(),record.size(),MPI_CHAR,records.data(),record.size(),MPI_CHAR,0,MPI_COMM_WORLD);
  bool ready=true;
  if(rank==0) {
    std::ofstream manifest(base+".tmp");manifest.imbue(std::locale::classic());
    manifest<<signature(o,rhs)<<"residual\t"<<std::setprecision(17)<<residual<<'\n';
    for(int r=0;r<ranks;++r) manifest<<&records[std::size_t(r)*record.size()]<<'\n';
    manifest.close();ready=bool(manifest) && std::rename((base+".tmp").c_str(),(base+".ready").c_str())==0;
  }
  if(!all_ok(ready)) throw CollectiveError("cannot commit checkpoint manifest");
}

inline double read(const SolveOptions& o,int rhs,int rank,int ranks,const std::string& coordinates,
                   const std::string& gauge_hash,double* data,std::size_t count) {
  const std::string base=prefix(o.resume_from,rhs),expected=signature(o,rhs);
  std::ifstream manifest(base+".ready",std::ios::binary);manifest.imbue(std::locale::classic());
  std::string header(expected.size(),'\0');manifest.read(&header[0],header.size());
  bool ok=bool(manifest) && header==expected;
  std::string key,line,field_hash;double residual=-1;
  if(!std::getline(manifest,line)) ok=false;
  std::istringstream metric(line);metric.imbue(std::locale::classic());
  if(!(metric>>key>>residual) || key!="residual" || !std::isfinite(residual) || residual<0 || (metric>>key)) ok=false;
  for(int r=0;r<ranks;++r) {
    if(!std::getline(manifest,line)) {ok=false;break;}
    std::istringstream row(line);int index=-1;std::string where,gauge,field,extra;
    if(!(row>>index>>where>>gauge>>field) || index!=r || field.size()!=64 || gauge.size()!=64 || (row>>extra)) {ok=false;break;}
    if(r==rank) {ok=ok && where==coordinates && gauge==gauge_hash;field_hash=field;}
  }
  if(std::getline(manifest,line)) ok=false;
  if(!all_ok(ok)) throw CollectiveError("checkpoint manifest, source, operator, gauge or MPI layout mismatch");
  const std::size_t bytes=count*sizeof(double);
  std::ifstream shard(base+"-"+coordinates+".f64",std::ios::binary);
  shard.read(reinterpret_cast<char*>(data),bytes);
  ok=bool(shard) && shard.peek()==std::char_traits<char>::eof() && sha256(data,bytes)==field_hash;
  if(!all_ok(ok)) throw CollectiveError("checkpoint shard length or SHA-256 mismatch");
  return residual;
}
}
