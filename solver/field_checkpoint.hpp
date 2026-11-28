#pragma once
#include "checkpoint_sha256.hpp"
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
