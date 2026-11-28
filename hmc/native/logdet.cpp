#include "runtime/bridge_setup.h"
#include "operator/fopr_Clover.h"
#include "operator/fopr_Smeared.h"
#include "force/force_F_Smeared.h"
#include "gauge/forceSmear_APE.h"
#include "gauge/projection_Stout_SU3.h"
#include "gauge/smear_APE.h"
#include "runtime/threadManager.h"
#include "field.hpp"
#include "local_block.hpp"
#include <fstream>
#include <iomanip>
#include <vector>

void evaluate(double mass,int nf,bool unit,double h)
{
  const double kap=0.5/(4.+mass);
  Parameters sp; sp.set_double("rho_uniform",rho);
  Projection_Stout_SU3 proj; Smear_APE smear(&proj,sp); ForceSmear_APE ds(&proj,sp);
  Director_Smear director(&smear); director.set_parameters(1);
  Fopr_Clover kernel("Dirac"); kernel.set_parameters(kap,coeff,{1,1,1,-1});
  Fopr_Smeared op(&kernel,&director);
  Field_G u(volume,4),f(volume,4),tmp(volume,4);
  Field_F e(volume,1),v(volume,1),g5v(volume,1);
  for(int i=0;i<volume;++i) for(int mu=0;mu<4;++mu) {
    const Matrix m=unit ? eye() : link_matrix(i,mu,h);
    for(int n=0;n<9;++n) u.set_ri(n,i,mu,m[n].real(),m[n].imag());
  }
  op.set_config(&u); op.set_mode("D");
  FixedGaugeCloverForce derivative(director.get_config(),kap,coeff,{1,1,1,-1});
  Force_F_Smeared force(&derivative,&ds,&director);
  for(int cb=0;cb<2;++cb) {
    // Even local x extents give exactly one site of this parity per x pair.
    std::vector<Mat12> blocks(half_volume);
    for(int j=0;j<12;++j) {
      e.set(0.);
      for(int i=0;i<volume;++i) if(parity(i)==cb) e.set_ri(j%3,j/3,i,0,1.,0.);
      op.mult(v,e);
      for(int i=0;i<volume;++i) if(parity(i)==cb) for(int k=0;k<12;++k)
        blocks[i/2][12*k+j] = Z(v.cmp_r(k%3,k/3,i),v.cmp_i(k%3,k/3,i));
    }
    long double total=0;
    for(auto& block:blocks) total-=nf*factor(block);
    total=hmc_mpi::sum_determinant(total);
    hmc_mpi::report("LOGDET %d %d %d %.17g %.17g\n",nf,cb,int(unit),h,double(total));
    if(!unit && h==0.) {
      f.set(0.); force.set_config(&u);
      for(int j=0;j<12;++j) {
        e.set(0.); v.set(0.);
        for(int i=0;i<volume;++i) if(parity(i)==cb) {
          e.set_ri(j%3,j/3,i,0,1.,0.);
          auto col=inverse_column(blocks[i/2],j);
          for(int k=0;k<12;++k) v.set_ri(k%3,k/3,i,0,col[k].real(),col[k].imag());
        }
#pragma omp parallel if(ThreadManager::get_num_threads_available() > 1)
        {
          kernel.mult_gm5(g5v,v);
          force.force_core1(tmp,g5v,e);
          axpy(f,double(nf),tmp);
        }
      }
      hmc_mpi::report("DERIV %d %d %.17g\n",nf,cb,f.cmp_r(3,0,0)-f.cmp_r(1,0,0));
      save_force(f,"nf"+std::to_string(nf)+"_cb"+std::to_string(cb)+".bin");
    }
  }
}
int main(int argc,char** argv) {
initialize(argc,argv,"CLOVER_HMC_LOGDET_THREADS");
int mode; double num,den;
  while(hmc_mpi::request(mode,num,den,target,max_iterations)) {
    read_binary_field(gauge_store,"u.bin");
evaluate(light_mass,2,false,0.); evaluate(-.24,1,false,0.);
done(); }
bridge_finalize();
}
