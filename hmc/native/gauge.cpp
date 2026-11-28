#include "lib/bridge_setup.h"
#include "Action/Gauge/action_G_Rectangle.h"
#include "field.hpp"
#include "ResourceManager/threadManager.h"
#include <cstdio>
#include <fstream>
#include <iomanip>

void do_drift(double dt) {
  std::vector<Matrix> mom((4*volume)), result((4*volume));
  read_binary_field(mom,"p.bin");
  Field_G p(volume,4),u(volume,4),out(volume,4);
#pragma omp parallel for schedule(static) if(ThreadManager::get_num_threads_available() > 1)
  for(int x=0;x<volume;++x) for(int mu=0;mu<4;++mu) for(int c=0;c<9;++c) {
    auto a=mom[x*4+mu][c], b=gauge_store[x*4+mu][c];
    p.set_ri(c,x,mu,a.real(),a.imag()); u.set_ri(c,x,mu,b.real(),b.imag());
  }
  // Bridge++ partitions sites among the threads calling this routine.
#pragma omp parallel if(ThreadManager::get_num_threads_available() > 1)
  {
    mult_exp_Field_G(out,dt,p,u,20);
  }
#pragma omp parallel for schedule(static) if(ThreadManager::get_num_threads_available() > 1)
  for(int x=0;x<volume;++x) for(int mu=0;mu<4;++mu) for(int c=0;c<9;++c)
    result[4*x+mu][c]={out.cmp_r(c,x,mu),out.cmp_i(c,x,mu)};
  write_binary_field(result);
}

int main(int argc, char **argv)
{
  initialize(argc,argv,"CLOVER_HMC_GAUGE_THREADS");
  Parameters pars;
  pars.set_double("beta", 6.2);
  pars.set_double("c_plaq", 1.0);
  pars.set_double("c_rect", -1.0/(20*0.855520*0.855520));
  pars.set_string("verbose_level", "Crucial");
  Action_G_Rectangle action(pars);
  Field_G u(volume, 4), force(volume, 4);
  action.set_config(&u);
  u.set_unit();
  const double baseline = action.calcH();
  hmc_mpi::report("BASELINE %.17g\n", baseline);
int mode; double num,den;
  while(hmc_mpi::request(mode,num,den,target,max_iterations)) {
    read_binary_field(gauge_store,"u.bin");
    if(mode==3) {
      Field_G p(volume,4); std::vector<Matrix> out((4*volume));
      const double scale=std::sqrt(0.5);
#pragma omp parallel for schedule(static) if(ThreadManager::get_num_threads_available() > 1)
      for(int x=0;x<volume;++x) for(int mu=0;mu<4;++mu) {
        for(int k=0;k<9;++k) {
          const auto a=gauge_store[4*x+mu][k]*scale;
          p.set_ri(k,x,mu,a.real(),a.imag());
        }
        auto a=p.mat(x,mu); a.at(); p.set_mat(x,mu,a);
        for(int k=0;k<9;++k) out[4*x+mu][k]={p.cmp_r(k,x,mu),p.cmp_i(k,x,mu)};
      }
      write_binary_field(out); done();continue;
    }
    if(mode==1) { do_drift(num); done(); continue; }
  for (double h : {0.0}) {
#pragma omp parallel for schedule(static) if(ThreadManager::get_num_threads_available() > 1)
    for (int site = 0; site < volume; ++site) {
      for(int mu=0;mu<4;++mu) {
        const auto mat=link_matrix(site,mu,0.);
        for(int n=0;n<9;++n) u.set_ri(n,site,mu,mat[n].real(),mat[n].imag());
      }
    }
    hmc_mpi::report("ACTION %.17g %.17g\n", h, action.calcH()-baseline);
    if (h == 0) {
      action.force(force);
      const double deriv = force.cmp_i(4,0,0) - force.cmp_i(0,0,0);
      hmc_mpi::report("DERIVATIVE %.17g\n", deriv);
      save_force(force,"force.bin");
    }
  }
    done();
  }
  bridge_finalize();
}
