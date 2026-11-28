#pragma once
#include "lib/bridge_setup.h"
#include "Field/index_eo.h"
#include "Fopr/fopr_Clover.h"
#include "Tools/gammaMatrixSet_Dirac.h"
#include "Fopr/fopr_CloverTerm_eo.h"
#include "Fopr/fopr_Wilson_eo.h"
#include "local_block.hpp"
#include "Force/Fermion/force_F_Smeared.h"
#include "Smear/forceSmear_APE.h"
#include "Smear/projection_Stout_SU3.h"
#include "Smear/smear_APE.h"
#include "field.hpp"
#include "ResourceManager/threadManager.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <memory>

struct OddSchur {
  Fopr_Wilson_eo hop;
  Fopr_CloverTerm_eo block;
  GammaMatrix gamma5;
  std::unique_ptr<Force_F_Clover_Nf2> kernel_force;
  std::unique_ptr<Force_F_Smeared> chain;
  Index_eo idx;
  OddSchur(double mass,Field_G& thin,Director_Smear& dr,ForceSmear_APE* ds) :
    hop("Dirac"),block("Dirac"),gamma5(GammaMatrixSet_Dirac().get_GM(GammaMatrixSet::GAMMA5)),
    scratch_a(half_volume,1),scratch_b(half_volume,1) {
    const std::vector<int> bc={1,1,1,-1};
    const double k=.5/(4.+mass);
    hop.set_parameters(k,bc); block.set_parameters(k,coeff,bc);
    dr.set_config(&thin);
    auto fat=dr.get_config();
    hop.set_config(fat); block.set_config(fat);
    // Heatbath inversions only apply the Schur operator. A force workspace
    // includes its own full Clover operator and smeared-link derivatives.
    if(ds) {
      kernel_force.reset(new FixedGaugeCloverForce(fat,k,coeff,bc));
      chain.reset(new Force_F_Smeared(kernel_force.get(),ds,&dr));
      chain->set_config(&thin);
    }
  }
  void apply(Field_F& out,const Field_F& in,bool dagger) {
#pragma omp parallel if(ThreadManager::get_num_threads_available() > 1)
    { apply_team(out,in,dagger); }
  }
  // Keep both Schur applications in one OpenMP team and retain the image
  // for callers that also need ||M*in|| or the RHMC force.
  void apply_normal(Field_F& out,Field_F& image,const Field_F& in) {
#pragma omp parallel if(ThreadManager::get_num_threads_available() > 1)
    {
      apply_team(image,in,false);
#pragma omp barrier
      apply_team(out,image,true);
    }
  }
  // Accumulate +2 Re[left^dag (delta D) right] on the original thin field.
  void bilinear(Field_G& out,const Field_F& left,int cl,const Field_F& right,int cr,double weight) {
    Field_F l(volume,1),r(volume,1),gl(volume,1);
    Field_G tmp(volume,4);
    // Every member enters the cooperative Clover and Stout force kernels
    // with the same fields. Finish reading tmp before the next force writes it.
#pragma omp parallel if(ThreadManager::get_num_threads_available() > 1)
    {
      l.set(0.); r.set(0.);
      idx.reverseField(l,left,cl); idx.reverseField(r,right,cr);
      mult_GM(gl,gamma5,l);
      chain->force_core1(tmp,gl,r); axpy(out,-weight,tmp);
#pragma omp barrier
      chain->force_core1(tmp,r,gl); axpy(out,-weight,tmp);
    }
  }
  void derivative(Field_G& out,const Field_F& left,const Field_F& right,double w) {
    Field_F t(half_volume,1),u(half_volume,1),v(half_volume,1),lo(half_volume,1),le(half_volume,1);
#pragma omp parallel if(ThreadManager::get_num_threads_available() > 1)
    {
      hop.Meo(t,right,0); block.mult_csw_inv(u,t,0);
      hop.Meo(t,u,1); block.mult_csw_inv(v,t,1);
      block.mult_csw_inv(lo,left,1); hop.Mdageo(t,lo,0); block.mult_csw_inv(le,t,0);
    }
    bilinear(out,lo,1,v,1,w);
    bilinear(out,lo,1,u,0,-w);
    bilinear(out,le,0,u,0,w);
    bilinear(out,le,0,right,1,-w);
  }

 private:
  void apply_team(Field_F& out,const Field_F& in,bool dagger) {
    auto& a=scratch_a;
    auto& b=scratch_b;
    if(!dagger) {
      hop.Meo(a,in,0); block.mult_csw_inv(b,a,0);
      hop.Meo(a,b,1); block.mult_csw_inv(out,a,1);
    } else {
      block.mult_csw_inv(a,in,1); hop.Mdageo(b,a,0);
      block.mult_csw_inv(a,b,0); hop.Mdageo(out,a,1);
    }
    aypx(-1.,out,in);
  }
  Field_F scratch_a,scratch_b;
};
