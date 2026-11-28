#include "schur.hpp"
#include "spinor.hpp"
#include "schur.hpp"
#include <deque>
#include <map>
#include <utility>

// Least-residual chronological predictor in the real inner product on complex fields.
// Recompute each image with the current gauge. Distinct mass pairs keep separate histories.
using SolutionHistory=std::deque<Field_F>;
static std::map<std::pair<double,double>,SolutionHistory> past_solutions;

inline int chronological_guess(OddSchur& M,const SolutionHistory& stored,const Field_F& rhs,
                               Field_F& x,Field_F& residual) {
  std::vector<Field_F> directions,images;
  Field_F y(half_volume,1);
  const double rhs_norm=rhs.norm();
  for(auto prev=stored.rbegin();prev!=stored.rend();++prev) {
    Field_F v(*prev),a(half_volume,1);
    M.apply(y,v,false);M.apply(a,y,true);
    const double before=a.norm();
    if(!(before>0.)) continue;
    // Modified Gram-Schmidt twice for nearly collinear neighboring MD solutions.
    for(int pass=0;pass<2;++pass) for(size_t j=0;j<images.size();++j) {
      const double coefficient=dot(images[j],a);
      axpy(a,-coefficient,images[j]);axpy(v,-coefficient,directions[j]);
    }
    const double length=a.norm();
    if(length<1.e-10*before) continue;
    scal(a,1./length);scal(v,1./length);
    const double coefficient=dot(a,rhs);
    axpy(x,coefficient,v);
    images.push_back(a);directions.push_back(v);
  }
  // Evaluate b-MdagM*x explicitly before CG; a weak or degenerate prediction uses zero.
  M.apply(y,x,false);M.apply(residual,y,true);aypx(-1.,residual,rhs);
  const double predicted=residual.norm();
  if(!std::isfinite(predicted) || predicted>rhs_norm) {
    x.set(0.);copy(residual,rhs);return 0;
  }
  return static_cast<int>(images.size());
}

#include "schur.hpp"
#include "../../solver/clover_sap_geometry.hpp"
#include <cstdlib>
#if defined(USE_ALT_CODE) && defined(USE_ALT_QXS)
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Fopr/afopr_Clover.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"
#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt/Solver/asolver_MG.h"
#endif

namespace hmc_light {
inline bool use_multigrid() {
  const char* value=std::getenv("CLOVER_HMC_LIGHT_SOLVER");
  return value && std::string(value)=="mg";
}
inline std::vector<int> sap_blocks(const std::vector<int>& dims) {
#if defined(USE_ALT_CODE) && defined(USE_ALT_QXS)
  return clover_qxs::sap_blocks(dims, VLENXS, VLENYS);
#else
  throw std::runtime_error("light MG requires a Bridge++ QXS build");
#endif
}

// M=A_o^{-1} S for S the odd Schur complement.  The odd block of the
// inverse of the full D is S^{-1}. D^dag=gamma5 D gamma5 on the fat links.
inline int correct_normal(Field_F& x, const Field_F& residual, OddSchur& op,
                          Field& fat, double mass, double tolerance) {
#if defined(USE_ALT_CODE) && defined(USE_ALT_QXS)
  const auto block = sap_blocks(dimensions);
  using DoubleField=AField<double,QXS>;
  using FloatField=AField<float,QXS>;
  Parameters fermion,level,settings;
  fermion.set_string("fermion_type","Clover");
  fermion.set_string("gamma_matrix_type","Dirac");
  fermion.set_double("hopping_parameter",.5/(4.+mass));
  fermion.set_double("clover_coefficient",coeff);
  fermion.set_int_vector("boundary_condition",{1,1,1,-1});
  fermion.set_string("verbose_level","Crucial");
  fermion.set_int_vector("sap_block",block);
  fermion.set_int_vector("block_size",block);
  level.set_int_vector("sap_block",block);
  level.set_int("setup_number_of_vectors",24);
  level.set_int("setup_number_of_step",1);
  level.set_int("maximum_number_of_iteration",200);
  level.set_int("maximum_number_of_restart",1);
  level.set_double("convergence_criterion_squared",2e-2);
  level.set_int("smoother_number_of_iteration",4);
  level.set_double("smoother_convergence_criterion_squared",1e-3);
  level.set_string("verbose_level","Crucial");
  settings.set_Parameters("MultiGrid_Level1",level);
  settings.set_int("maximum_number_of_iteration",500);
  settings.set_int("maximum_number_of_restart",10);
  settings.set_double("convergence_criterion_squared",.01*tolerance*tolerance);
  settings.set_string("verbose_level","Crucial");

  AFopr_Clover<DoubleField> dirac(fermion);
  AFopr_Clover_dd<FloatField> packed(fermion);
  ASolver_MG<DoubleField> solver;
  dirac.set_config(&fat);dirac.set_mode("D");
  packed.set_config(&fat);
  solver.set_parameters(settings);solver.set_foprD(&dirac);solver.set_foprF(&packed);
  solver.init_solver();solver.run_setup();

  Field_F source(volume,1),full(volume,1),odd(half_volume,1),weighted(half_volume,1);
  DoubleField b(source.nin(),volume,1),solution(source.nin(),volume,1);
  int calls=0;
  auto invert=[&](const Field_F& in,Field_F& out) {
    int count=0;double diff=0;
#pragma omp parallel
    {
      dirac.convert(b,in);solution.set(0.0);
      solver.solve(solution,b,count,diff);
      dirac.reverse(out,solution);
    }
    calls+=count;
  };
  source.set(0.0);op.idx.reverseField(source,residual,1);
  mult_GM(full,op.gamma5,source);
  invert(full,source);
  mult_GM(full,op.gamma5,source);
  op.idx.convertField(odd,full,1);
  op.block.D(weighted,odd,1);      // (M^dag)^{-1}r = A_o S^{-dag} r
  op.block.D(odd,weighted,1);      // RHS for M^{-1}: A_o (M^dag)^{-1}r
  source.set(0.0);op.idx.reverseField(source,odd,1);
  invert(source,full);
  op.idx.convertField(odd,full,1);
  axpy(x,1.0,odd);
  return calls;
#else
  throw std::runtime_error("light MG requires a Bridge++ QXS build");
#endif
}
}


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
