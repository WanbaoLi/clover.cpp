#pragma once
#include "Force/Fermion/force_F_Clover_Nf2.h"
#include <stdexcept>
#include <array>
#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

using Z = std::complex<double>;
using Mat12 = std::array<Z,144>;
inline double factor(Mat12& block)
{
  double logdet=0;
  for(int i=0;i<12;++i) for(int j=0;j<=i;++j) {
    if(std::abs(block[12*i+j]-std::conj(block[12*j+i]))>1e-12)
      throw std::runtime_error("Clover block is not Hermitian");
    Z s=block[12*i+j];
    for(int k=0;k<j;++k) s-=block[12*i+k]*std::conj(block[12*j+k]);
    if(i==j) {
      if(s.real()<=0 || std::abs(s.imag())>1e-12) throw std::runtime_error("Nonpositive pivot");
      block[12*i+j]=std::sqrt(s.real());
      logdet+=std::log(s.real());
    } else {
      block[12*i+j]=s/block[12*j+j].real();
    }
  }
  return logdet;
}
inline std::array<Z,12> inverse_column(const Mat12& l,int col)
{
  std::array<Z,12> v{};
  for(int i=0;i<12;++i) {
    Z s=(i==col) ? 1. : 0.;
    for(int j=0;j<i;++j) s-=l[12*i+j]*v[j];
    v[i]=s/l[12*i+i].real();
  }
  for(int i=11;i>=0;--i) {
    Z s=v[i];
    for(int j=i+1;j<12;++j) s-=std::conj(l[12*j+i])*v[j];
    v[i]=s/l[12*i+i].real();
  }
  return v;
}

// One HMC force request supplies an immutable smeared gauge field. Bridge's
// smearing wrapper revisits set_config before each bilinear; the Clover term
// and its staples can be prepared once for this field and this mass.
class FixedGaugeCloverForce final : public Force_F_Clover_Nf2 {
  Field* const links_;
public:
  FixedGaugeCloverForce(Field_G* links, double kappa, double csw,
                        const std::vector<int>& boundary)
    : Force_F_Clover_Nf2("Dirac"), links_(links) {
    set_parameters(kappa,csw,boundary);
    Force_F_Clover_Nf2::set_config(links);
  }
  void set_config(Field* links) override {
    if(links!=links_) throw std::logic_error("HMC force used with a different gauge field");
  }
};
