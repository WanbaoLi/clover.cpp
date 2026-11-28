#include "schur.hpp"
#include "spinor.hpp"
#include "coefficients.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>

static int refresh_mode=0;
namespace {
struct ShiftState {
  explicit ShiftState(const Field_F& rhs) : x(half_volume,1), direction(half_volume,1) {
    x.set(0.0);
    copy(direction,rhs);
  }
  Field_F x,direction;
  double z=1.0,previous_z=1.0;
  int updates=0;
  bool active=true;
};
using States=std::vector<std::unique_ptr<ShiftState>>;

void normal(OddSchur& M,Field_F& result,Field_F& tmp,const Field_F& x,double shift) {
  M.apply(tmp,x,false);
  M.apply(result,tmp,true);
  axpy(result,shift,x);
}

int independent_cg(OddSchur& M,ShiftState& state,const Field_F& initial_residual,
                   double shift,double target_sq,int limit) {
  Field_F r(half_volume,1),tmp(half_volume,1),product(half_volume,1);
  copy(r,initial_residual);
  copy(state.direction,r);
  double rr=r.norm2();
  int steps=0;
  while (rr>target_sq && steps<limit) {
    normal(M,product,tmp,state.direction,shift);
    const double denominator=tmp.norm2()+shift*state.direction.norm2();
    if (!(denominator>0.0) || !std::isfinite(denominator)) throw std::runtime_error("invalid shifted-CG curvature");
    const double step=rr/denominator;
    axpy(state.x,step,state.direction);
    axpy(r,-step,product);
    const double next=r.norm2();
    aypx(next/rr,state.direction,r);
    rr=next;
    ++steps;
  }
  return steps;
}

// Shifted residuals stay collinear with the residual of the smallest positive
// shift. One application of M^dag M supplies one step for every active pole.
int shared_cg(OddSchur& M,States& states,const std::vector<double>& poles,
              const Field_F& rhs,double target_sq,int limit) {
  if (!std::is_sorted(poles.begin(),poles.end()) || !(poles.front()>0.0))
    throw std::runtime_error("multi-shift CG requires sorted positive poles");
  Field_F r(half_volume,1),tmp(half_volume,1),product(half_volume,1);
  copy(r,rhs);
  ShiftState& base=*states.front();
  double rr=r.norm2(),previous_step=1.0,previous_beta=0.0;
  while (rr>target_sq && base.updates<limit) {
    normal(M,product,tmp,base.direction,poles.front());
    const double curvature=tmp.norm2()+poles.front()*base.direction.norm2();
    if (!(curvature>0.0) || !std::isfinite(curvature)) throw std::runtime_error("invalid multi-shift CG curvature");
    const double step=rr/curvature;
    axpy(base.x,step,base.direction);
    axpy(r,-step,product);
    const double next=r.norm2(),beta=next/rr;
    aypx(beta,base.direction,r);
    ++base.updates;
    const double coupling=previous_beta*step/previous_step;
    for (size_t i=1;i<poles.size();++i) {
      ShiftState& s=*states[i];
      if (!s.active) continue;
      const double reciprocal=(1.0+coupling+(poles[i]-poles.front())*step)/s.z-coupling/s.previous_z;
      if (!(reciprocal>0.0) || !std::isfinite(reciprocal)) throw std::runtime_error("invalid shifted residual recurrence");
      const double z=1.0/reciprocal,ratio=z/s.z;
      axpy(s.x,step*ratio,s.direction);
      ++s.updates;
      if (z*z*next<=target_sq) s.active=false;
      else {
        scal(s.direction,beta*ratio*ratio);
        axpy(s.direction,z,r);
      }
      s.previous_z=s.z;
      s.z=z;
    }
    previous_step=step;
    previous_beta=beta;
    rr=next;
  }
  return base.updates;
}

bool shared_mode() {
  const char* env=std::getenv("CLOVER_HMC_STRANGE_SOLVER");
  const std::string choice=env?env:"shared";
  if (choice!="independent" && choice!="shared") throw std::runtime_error("CLOVER_HMC_STRANGE_SOLVER must be shared or independent");
  return choice=="shared";
}
}
void evaluate(double h,double tol) {
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
  OddSchur M(strange_mass,u,dr,ds.get());
  std::unique_ptr<Field_G> f;
  if(refresh_mode!=2) { f.reset(new Field_G(volume,4)); f->set(0.); }
  auto rational=[&](Field_F& out,const char* tag,const Rational& rt,const Field_F& rhs,bool force) {
    copy(out,rhs); scal(out,rt.c);
    const bool shared=shared_mode();
    const double norm_sq=rhs.norm2(),target_sq=norm_sq*tol*tol;
    if (!std::isfinite(norm_sq) || !(tol>0) || max_iterations<1 || rt.pole.empty() || rt.pole.size()!=rt.res.size())
      throw std::runtime_error("invalid rational solver input");
    States states;
    int steps=0,polish=0,checks=0;
    for (size_t i=0;i<rt.pole.size();++i) states.emplace_back(new ShiftState(rhs));
    if (shared) steps=shared_cg(M,states,rt.pole,rhs,target_sq,max_iterations);
    for(size_t i=0;i<rt.pole.size();++i) {
      ShiftState& s=*states[i];
      const double pole=rt.pole[i];
      if (!shared) {
        s.updates=independent_cg(M,s,rhs,pole,target_sq,max_iterations);
        steps+=s.updates;
      }
      Field_F y(half_volume,1),q(half_volume,1);
      normal(M,q,y,s.x,pole);
      ++checks;
      axpy(q,-1.,rhs);
      double relative=norm_sq? q.norm()/std::sqrt(norm_sq) : 0.0;
      if (relative>tol && s.updates<max_iterations) {
        scal(q,-1.);
        const int extra=independent_cg(M,s,q,pole,target_sq,max_iterations-s.updates);
        polish+=extra;
        s.updates+=extra;
        normal(M,q,y,s.x,pole);
        ++checks;
        axpy(q,-1.,rhs);
        relative=q.norm()/std::sqrt(norm_sq);
      }
      hmc_mpi::report("SHIFT %s %.17g %zu %.17g %.17g %d\n",tag,h,i,pole,relative,s.updates);
      if (!std::isfinite(relative)) throw std::runtime_error("nonfinite true shifted residual");
      axpy(out,rt.res[i],s.x);
      if(force) M.derivative(*f,y,s.x,-rt.res[i]);
    }
    hmc_mpi::report("KRYLOV %s %s %zu %d %d %d\n",tag,shared?"shared":"independent",rt.pole.size(),steps,polish,checks);
  };
  if(refresh_mode==2) {
    Field_F heat(half_volume,1);rational(heat,"heat",action_ipfe,phi,false);write_spinor("refreshed.bin",heat);return;
  }
  // A force-gradient midpoint needs only the MD rational function.
  std::unique_ptr<Field_F> action;
  if (refresh_mode != 3) {
    action.reset(new Field_F(half_volume,1));
    rational(*action,"action",action_pfe,phi,false);
  }
  Field_F md(half_volume,1);
  rational(md,"md",force_pfe,phi,h==0.);
  if (action) {
    const double norm = action->norm2();
    hmc_mpi::report("ACTION %.17g %.17g %.17g %.17g\n",h,norm-global_volume*6,norm,dot(phi,md));
  }
  if(h!=0.) return;
  hmc_mpi::report("DERIV %.17g\n",f->cmp_r(3,0,0)-f->cmp_r(1,0,0));
  save_force(*f,"force.bin");
}
int main(int argc,char** argv) {
initialize(argc,argv,"CLOVER_HMC_FERMION_THREADS");
int mode; double num,den;
  while(hmc_mpi::request(mode,num,den,target,max_iterations)) {
    read_binary_field(gauge_store,"u.bin"); refresh_mode=mode;
evaluate(0.,target);
done(); }
bridge_finalize();
}
