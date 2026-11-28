#include "schur.hpp"
#include "spinor.hpp"
#include <vector>
struct Rational {double c; std::vector<double> res,pole;};
const double strange_mass=-0.2400;
const Rational action_pfe={0.163770623718441,{0.000456950593274148,0.000933070019211125,0.00173260367764942,0.00325334197095044,0.00618223198007706,0.0118237931424094,0.0226893233800592,0.0436470377125287,0.0842552605967904,0.163834105667813,0.324048670587981,0.667674513073305,1.52011890693849,4.46787289574894,28.5894502783284},{5.55835590775444e-05,0.000416755864364307,0.00139575880817907,0.00378307619387897,0.00950643266149974,0.023189559723109,0.0558993621083468,0.134164184434855,0.321874947989621,0.774669290837093,1.8820254251941,4.68164410178872,12.3608907310646,38.3422607790049,204.873221861185}};
const Rational action_ipfe={6.10610118771501,{-5.90262826538435e-06,-2.63363387226834e-05,-8.62160355606352e-05,-0.000263984258286453,-0.000792810319715722,-0.00236581977385576,-0.00704746125114149,-0.0210131715847004,-0.0629242233443976,-0.190538104129215,-0.592816342814611,-1.96992441194278,-7.70705574740274,-46.55440910469,-1281.70053339288},{0.000109335909283339,0.000584211769074023,0.00181216713967916,0.00478464392272826,0.0119020708754186,0.0289155646996088,0.0695922442548162,0.166959610676697,0.400720136243831,0.965951931276981,2.35629923417205,5.92110728201649,16.0486180482883,53.7484938194392,402.99686403222}};
const Rational force_pfe={0.0303920178459613,{0.00943108618345698,0.0122499930158508,0.0187308029056777,0.0308130330025528,0.0521206555919226,0.0890870585774984,0.153090120000215,0.26493803350899,0.466760251501358,0.866223656646014,1.8819154073627,6.96033769739192},{5.23045292201785e-05,0.000569214182255549,0.00226724207135389,0.00732861083302471,0.0222608882919378,0.0662886891030569,0.196319420401789,0.582378159903323,1.74664271771668,5.42569216297222,18.850085313508,99.6213166072174}};

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
  M.apply_normal(result,tmp,x);
  axpy(result,shift,x);
}

// Both norms enter the same CG curvature. Reduce the two real scalars
// together while preserving the local site order of Field::norm2().
double shifted_curvature(const Field_F& image,const Field_F& direction,double shift) {
  const double* y=image.ptr(0);
  const double* p=direction.ptr(0);
  double local[2]={0.,0.},global[2];
  for (int k=0;k<image.size();++k) {
    local[0]+=y[k]*y[k];
    local[1]+=p[k]*p[k];
  }
  Communicator::reduce_sum(2,global,local);
  return global[0]+shift*global[1];
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
    const double denominator=shifted_curvature(tmp,state.direction,shift);
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
    const double curvature=shifted_curvature(tmp,base.direction,poles.front());
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
    Field_F y(half_volume,1),q(half_volume,1);
    for(size_t i=0;i<rt.pole.size();++i) {
      ShiftState& s=*states[i];
      const double pole=rt.pole[i];
      if (!shared) {
        s.updates=independent_cg(M,s,rhs,pole,target_sq,max_iterations);
        steps+=s.updates;
      }
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
