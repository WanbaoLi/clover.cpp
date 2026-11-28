// On the unit gauge, a plane wave at temporal momentum pi/T has an analytic D norm.
#include "runtime/bridge_setup.h"
#include "operator/fopr_Clover.h"
#include "field/field_F.h"
#include "field/field_G.h"
#include <cmath>
#include <cstdio>
int main(int argc,char** argv) {
  bridge_initialize(&argc,&argv);
  const int extent=4;
  bridge_setup({4,4,4,extent},{1,1,1,1},1,3,"stdout","Crucial");
  double error=1.;
  {
    const double k=.1,p=std::acos(-1.)/extent;
    Field_G links(256,4);links.set_unit();
    Fopr_Clover op("Dirac");op.set_parameters(k,1.160678675188,{1,1,1,-1});
    op.set_config(&links);op.set_mode("D");
    Field_F in(256,1),out(256,1);
    in.set(0.);
    for(int site=0;site<256;++site)
      in.set_ri(0,0,site,0,std::cos(p*(site/64)),std::sin(p*(site/64)));
    op.mult(out,in);
    const double expected=std::pow(1.-2.*k*(3.+std::cos(p)),2)+std::pow(2.*k*std::sin(p),2);
    error=std::abs(out.norm2()/in.norm2()-expected);
    std::printf("plane_wave_norm2_absolute_error=%.17g\n",error);
  }
  bridge_finalize();
  return std::isfinite(error) && error<1e-13 ? 0 : 1;
}
