#pragma once
#include <complex>
#include <fstream>
using FieldComplex=std::complex<double>;
void load_spinor(const std::string& name,Field_F& phi) {
  Field_F lex(volume,1); lex.set(0.);
  if(hmc_mpi::ranks>1) {
    hmc_mpi::transfer(name.c_str(),lex.ptr(0),24,false);
    Index_eo ix; ix.convertField(phi,lex,1);
    return;
  }
  std::ifstream in(name,std::ios::binary);
  for(int i=0;i<volume;++i) for(int s=0;s<4;++s) for(int c=0;c<3;++c) {
    FieldComplex z; in.read(reinterpret_cast<char*>(&z),sizeof(z)); lex.set_ri(c,s,i,0,z.real(),z.imag());
  }
  if(!in || in.peek()!=EOF) throw std::runtime_error("spinor import");
  Index_eo ix; ix.convertField(phi,lex,1);
}
void write_spinor(const std::string& name,const Field_F& phi) {
  Field_F lex(volume,1); lex.set(0.);
  Index_eo ix; ix.reverseField(lex,phi,1);
  if(hmc_mpi::ranks>1) {
    hmc_mpi::transfer(name.c_str(),lex.ptr(0),24,true);
    return;
  }
  std::ofstream out(name,std::ios::binary);
  for(int i=0;i<volume;++i) for(int s=0;s<4;++s) for(int c=0;c<3;++c) {
    FieldComplex z(lex.cmp_r(c,s,i),lex.cmp_i(c,s,i)); out.write(reinterpret_cast<const char*>(&z),sizeof(z));
  }
  if(!out) throw std::runtime_error("spinor export");
}
