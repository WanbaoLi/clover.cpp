#include "schur.hpp"
#include "spinor.hpp"
#include "chrono.hpp"
#include "mg_light.hpp"

static int refresh_mode=0;
void evaluate(double h,double num,double den,bool ratio) {
  Parameters sp; sp.set_double("rho_uniform",rho);
  Projection_Stout_SU3 proj; Smear_APE smear(&proj,sp);
  std::unique_ptr<ForceSmear_APE> ds;
  if(refresh_mode!=2) ds.reset(new ForceSmear_APE(&proj,sp));
  Director_Smear dr(&smear); dr.set_parameters(1);
  Field_G u(volume,4);
  Field_F phi(half_volume,1);
  for(int i=0;i<volume;++i) {
    for(int mu=0;mu<4;++mu) {
      const auto m=link_matrix(i,mu,h);
      for(int n=0;n<9;++n) u.set_ri(n,i,mu,m[n].real(),m[n].imag());
    }
  }
  load_spinor("phi.bin",phi);
  const auto mass_pair=std::make_pair(num,den);
  if(refresh_mode==2) {
    past_solutions.erase(mass_pair);
    Field_F rhs(half_volume,1);
    {
      OddSchur M(num,u,dr,nullptr);
      M.apply(rhs,phi,true);
    }
    if(!ratio) {
      write_spinor("refreshed.bin",rhs);
      hmc_mpi::report("REFRESH 0 0\n");return;
    }
    // M^dag eta is the sole numerator input to the denominator inversion.
    // Releasing the numerator first avoids storing two sets of local inverses.
    Field_F out(half_volume,1);
    int iterations=0; double residual=0.;
    {
      OddSchur denominator(den,u,dr,nullptr);
      Field_F z(half_volume,1),p(half_volume,1),r(half_volume,1),y(half_volume,1),q(half_volume,1);
      copy(r,rhs);z.set(0.);const double initial=r.norm2();
      if(hmc_light::use_multigrid() && initial>0.) {
        const int calls=hmc_light::correct_normal(z,r,denominator,*dr.get_config(),den,target);
        denominator.apply_normal(q,y,z);copy(r,rhs);axpy(r,-1.,q);
        hmc_mpi::report("MG_NORMAL refresh %d %.17g\n",calls,r.norm()/std::sqrt(initial));
      }
      copy(p,r);double rr=r.norm2();
      while(rr>initial*target*target && iterations<max_iterations) {
        denominator.apply_normal(q,y,p);double alpha=rr/y.norm2();
        axpy(z,alpha,p);axpy(r,-alpha,q);double next=r.norm2();aypx(next/rr,p,r);rr=next;++iterations;
      }
      denominator.apply_normal(q,out,z);axpy(q,-1.,rhs);residual=q.norm()/std::sqrt(initial);
    }
    write_spinor("refreshed.bin",out);
    hmc_mpi::report("REFRESH %.17g %d\n",residual,iterations);return;
  }
  OddSchur M(num,u,dr,ds.get());
  std::unique_ptr<OddSchur> P;
  if(ratio) P.reset(new OddSchur(den,u,dr,ds.get()));
  Field_F rhs(half_volume,1),x(half_volume,1),r(half_volume,1),p(half_volume,1),y(half_volume,1),q(half_volume,1);
  if(ratio) P->apply(rhs,phi,true); else copy(rhs,phi);
  copy(r,rhs); copy(p,r); x.set(0.);
  const double initial=r.norm2();
  const int depth=refresh_mode>=11 && refresh_mode<=17 ? refresh_mode-10 : 0;
  int used=0;
  if(depth) {
    auto& history=past_solutions[mass_pair];
    if(!history.empty()) used=chronological_guess(M,history,rhs,x,r);
    copy(p,r);
  }
  double rr=r.norm2();const double starting_residual=std::sqrt(rr/initial);int it=0;
  if(hmc_light::use_multigrid() && rr>initial*target*target) {
    const int calls=hmc_light::correct_normal(x,r,M,*dr.get_config(),num,target);
    M.apply_normal(q,y,x);copy(r,rhs);axpy(r,-1.,q);
    rr=r.norm2();copy(p,r);
    hmc_mpi::report("MG_NORMAL action %d %.17g\n",calls,std::sqrt(rr/initial));
  }
  while(rr>initial*target*target && it<max_iterations) {
    M.apply_normal(q,y,p);
    const double alpha=rr/y.norm2();
    axpy(x,alpha,p); axpy(r,-alpha,q);
    const double next=r.norm2();
    aypx(next/rr,p,r); rr=next; ++it;
  }
  M.apply_normal(q,y,x); axpy(q,-1.,rhs);
  const double residual=q.norm()/rhs.norm();
  if(residual>1.5*target) {
    std::fprintf(stderr,"Schur true residual %.17g after %d iterations, target %.17g\n",residual,it,target);
    throw std::runtime_error("Schur true residual");
  }
  if(depth) {
    auto& history=past_solutions[mass_pair];
    history.push_back(x);
    while(static_cast<int>(history.size())>depth) history.pop_front();
    hmc_mpi::report("CHRONO %d %d %.17g %.17g %d\n",depth,used,starting_residual,residual,it);
  }
  // Subtract ||phi||^2 before the reduction by forming a per-component difference.
  Field_F px(half_volume,1);
  if(ratio) P->apply(px,x,false); else copy(px,x);
  axpy(px,-1.,phi);
  hmc_mpi::report("RATIO %.17g %.17g %.17g %.17g %d\n",h,dot(rhs,x),dot(phi,px),residual,it);
  if(h==0.) {
    Field_G f(volume,4); f.set(0.);
    if(ratio) P->derivative(f,phi,x,1.);
    M.derivative(f,y,x,-1.);
    hmc_mpi::report("DERIV %.17g\n",f.cmp_r(3,0,0)-f.cmp_r(1,0,0));
    save_force(f,"force.bin");
  }
}
int main(int argc,char** argv) {
  if(argc==6 && std::string(argv[1])=="--check-mg") {
    try {
      std::vector<int> dims;
      for(int i=2;i<6;++i) dims.push_back(std::stoi(argv[i]));
      hmc_light::sap_blocks(dims);
      return 0;
    } catch(const std::exception& ex) {
      std::fprintf(stderr,"%s\n",ex.what());return 2;
    }
  }
initialize(argc,argv,"CLOVER_HMC_FERMION_THREADS");
int mode; double num,den;
  while(hmc_mpi::request(mode,num,den,target,max_iterations)) {
    read_binary_field(gauge_store,"u.bin"); refresh_mode=mode;
evaluate(0.,num,den,den!=9.);
done(); }
bridge_finalize();
}
