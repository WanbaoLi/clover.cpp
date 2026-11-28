#include "solver_cli.hpp"
#include "field_checkpoint.hpp"
#include "mpi_layout.hpp"
#include "clover_sap_geometry.hpp"
#include "gauge_layout.hpp"
#include "configure.h"
#include "bridge_setup.h"
#include "bridge_init_factory.h"
#include "lib/Communicator/communicator.h"
#include "lib/Field/field_F.h"
#include "lib/Field/field_G.h"
#include "lib/Fopr/fopr.h"
#include "lib/IO/gaugeConfig.h"
#include "lib/Measurements/Fermion/source_Local.h"
#include "lib/Measurements/Fermion/source_Wall.h"
#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/Smear/projection_Stout_SU3.h"
#include "lib/Smear/smear_APE.h"
#include "lib/Tools/timer.h"
#include "lib/Solver/solver_BiCGStab_Cmplx.h"
#include "lib/Fopr/fopr_Clover_eo.h"
#if defined(USE_ALT_CODE) && defined(USE_ALT_QXS)
#define CLOVER_HAS_QXS 1
#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/afield-inc.h"
#include "lib_alt_QXS/Fopr/afopr_Clover.h"
#include "lib_alt_QXS/Fopr/afopr_Clover_dd.h"
#include "lib_alt_QXS/inline/define_vlen.h"
#include "lib_alt/Solver/asolver_MG.h"
#endif
#include <mpi.h>
#include <sys/stat.h>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#ifdef CLOVER_HAS_QXS
using DoubleField=AField<double,QXS>;
using FloatField=AField<float,QXS>;
#endif

static void require_all(bool local,const std::string& message) {
  int bad=local ? 0 : 1, any=0;
  MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if (any) throw clover_checkpoint::CollectiveError(message);
}
static void validate_geometry(const SolveOptions& o,int ranks) {
  long long volume=1,local=1,grid_volume=1;
  for (int i=0;i<4;++i) {
    if (o.grid[i]<1 || o.lattice[i]<1 || o.lattice[i]%o.grid[i] ||
        o.position[i]<0 || o.position[i]>=o.lattice[i])
      throw std::runtime_error("invalid MPI grid, lattice dimension or source position");
    const int n=o.lattice[i]/o.grid[i];
    if ((o.method=="eo-bicgstab") && (o.lattice[i]%2 || (i==0 && n%2)))
      throw std::runtime_error("even-odd preconditioning requires even global extents and an even local X extent");
    if (o.method=="mg") {
#ifdef CLOVER_HAS_QXS
    if (o.sap_block.empty()) {
      const int pack=i==0 ? VLENXS : (i==1 ? VLENYS : 1);
      if (n%((i<2 ? 2 : 4)*pack))
        throw std::runtime_error("local size is incompatible with SAP and QXS coarse packing");
    }
#else
    throw std::runtime_error("QXS is unavailable in this build; use --method bicgstab");
#endif
    }
    volume*=o.lattice[i];local*=n;grid_volume*=o.grid[i];
    // Field::ntot() stores a signed int; an SU(3) gauge has 18 reals in four directions.
    if (volume>std::numeric_limits<int>::max() || local>std::numeric_limits<int>::max()/72)
      throw std::runtime_error("global or local lattice exceeds Bridge++ field indexing; use more MPI ranks for a large local gauge");
  }
  if (grid_volume!=ranks) throw std::runtime_error("MPI process count differs from the product of --grid");
}
static std::ofstream open_output(const std::string& path) {
  std::ofstream f(path.c_str());
  f<<std::setprecision(17);
  return f;
}
#ifdef CLOVER_HAS_QXS
template<class FineField>
static double check_operator(AFopr_Clover<FineField>& qxs,Fopr& scalar,
    Field_F& b,Field_F& x,Field_F& ref,FineField& qb,FineField& qx,const std::vector<int>& L) {
  const int nx=CommonParameters::Nx(),ny=CommonParameters::Ny(),nz=CommonParameters::Nz(),nt=CommonParameters::Nt();
  const int x0=Communicator::ipe(0)*nx,y0=Communicator::ipe(1)*ny,z0=Communicator::ipe(2)*nz,t0=Communicator::ipe(3)*nt;
  const int volume=CommonParameters::Nvol();
#pragma omp parallel for schedule(static)
  for (int site=0;site<volume;++site) {
    const int xx=x0+site%nx,y=y0+(site/nx)%ny,z=z0+(site/(nx*ny))%nz,t=t0+site/(nx*ny*nz);
    const std::uint32_t index=24U*(((std::uint32_t(t)*L[2]+z)*L[1]+y)*L[0]+xx);
    double* values=b.ptr(24*site);
    for (int c=0;c<24;++c) {
      std::uint32_t h=(index+c+1U)*2654435761U;
      h^=h>>16;
      values[c]=std::ldexp(static_cast<double>(h),-32)-0.5;
    }
  }
  scalar.mult(ref,b);
#pragma omp parallel
  { qxs.convert(qb,b);qxs.mult(qx,qb);qxs.reverse(x,qx); }
  axpy(x,-1.0,ref);
  return std::sqrt(x.norm2()/ref.norm2());
}
class QxsSystem {
public:
  virtual ~QxsSystem() = default;
  virtual double solve(Field_F&, Field_F&, int&, double&, bool, int) = 0;
  std::vector<int> block;
  double operator_error=0,setup_seconds=0;
};

template<class FineField>
class QxsSystemImpl : public QxsSystem {
  AFopr_Clover<FineField> qxs;
  FineField qb,qx;
  std::unique_ptr<AFopr_Clover_dd<FloatField>> dd;
  std::unique_ptr<AFopr_dd<FloatField>> smoother;
  ASolver_MG<FineField> mg;

public:
  QxsSystemImpl(const SolveOptions& o,const Parameters& op,Field_G& gauge,Fopr& scalar,
            Field_F& b,Field_F& x,Field_F& res)
      : qxs(op),qb(b.nin(),b.nvol(),1),qx(b.nin(),b.nvol(),1) {
    qxs.set_config(&gauge);qxs.set_mode("D");
    operator_error=check_operator(qxs,scalar,b,x,res,qb,qx,o.lattice);
    const double operator_tolerance=sizeof(typename FineField::real_t)==sizeof(float) ? 1e-6 : 1e-12;
    if (!std::isfinite(operator_error) || operator_error>=operator_tolerance) throw std::runtime_error("independent scalar/QXS operator check failed");
  Parameters mgp,level;
  std::vector<int> local(4);
  for (int axis=0; axis<4; ++axis) local[axis]=o.lattice[axis]/o.grid[axis];
  block=o.sap_block.empty() ? clover_qxs::sap_blocks(local,VLENXS,VLENYS) : o.sap_block;
  mgp.set_int("maximum_number_of_iteration",o.mg_max_iters ? o.mg_max_iters : 500);
  mgp.set_int("maximum_number_of_restart",o.mg_max_iters ? 1 : 10);
  // Reserve room for roundoff between the recursive FP32 residual and the true residual.
  const double inner_tolerance=o.precision=="single" ? 0.25*o.tolerance : o.tolerance;
  mgp.set_double("convergence_criterion_squared",inner_tolerance*inner_tolerance);
  mgp.set_string("verbose_level","General");
  level.set_int_vector("sap_block",block);level.set_int("setup_number_of_vectors",o.vectors);level.set_int("setup_number_of_step",o.setup_steps);
  level.set_int("maximum_number_of_iteration",200);level.set_int("maximum_number_of_restart",1);
  level.set_double("convergence_criterion_squared",2e-2);level.set_int("smoother_number_of_iteration",o.smooth_steps);
  level.set_double("smoother_convergence_criterion_squared",1e-3);level.set_string("verbose_level","Crucial");
  mgp.set_Parameters("MultiGrid_Level1",level);
    Parameters ddp=op;ddp.set_int_vector("sap_block",block);ddp.set_int_vector("block_size",block);
    dd.reset(new AFopr_Clover_dd<FloatField>(ddp));dd->set_config(&gauge);
    mg.set_parameters(mgp);mg.set_foprD(&qxs);mg.set_foprF(dd.get());
    if (mg.use_fopr_for_smoother()) {
      smoother.reset(mg.new_fopr_smoother(ddp));smoother->set_config(&gauge);mg.set_fopr_smoother(smoother.get());
    }
    mg.init_solver();Timer setup("mg_setup");setup.start();mg.run_setup();setup.stop();
    setup_seconds=setup.elapsed_sec();
  }
  double solve(Field_F& x,Field_F& b,int& nconv,double& diff,bool profile,int rhs) override {
    typename FineField::real_t qxs_diff=0;
#pragma omp parallel
    { qxs.convert(qb,b);qx.set(0.0); }
    if (profile) mg.reset_flop_count();
    Timer timer("source_column");timer.start();
#pragma omp parallel
    { mg.solve(qx,qb,nconv,qxs_diff); }
    diff=qxs_diff;
    timer.stop();
    if (profile) {
      vout.general("MG solve profile for RHS %d\n",rhs);
      mg.flop_count();
    }
#pragma omp parallel
    { qxs.reverse(x,qx); }
    return timer.elapsed_sec();
  }
};
#endif

class EvenOddBiCGStab {
  Fopr_Clover_eo op;
  Solver_BiCGStab_Cmplx solver;
  Field even_source, odd_source, even_solution;

 public:
  EvenOddBiCGStab(const Parameters& parameters, Field_G& gauge,
                 double tolerance, int iterations_per_cycle, int maximum_cycles)
    : op(parameters), solver(&op) {
    op.set_config(&gauge);
    const int nin=2*CommonParameters::Nc()*CommonParameters::Nd(), half=gauge.nvol()/2;
    even_source.reset(nin,half,1);
    odd_source.reset(nin,half,1);
    even_solution.reset(nin,half,1);
    op.set_mode("D");
    Parameters settings;
    settings.set_int("maximum_number_of_iteration",iterations_per_cycle);
    settings.set_int("maximum_number_of_restart",maximum_cycles);
    settings.set_double("convergence_criterion_squared",tolerance*tolerance);
    settings.set_bool("use_initial_guess",true);
    settings.set_double("Omega_tolerance",0.0);
    settings.set_string("verbose_level","Detailed");
    solver.set_parameters(settings);
  }

  void load_even(const Field& solution, const Field& source) {
#pragma omp parallel
    {
      op.extract_even(even_solution,solution);
      op.preProp(even_source,odd_source,source);
    }
  }

  void restore(Field& solution, const Field& source) {
#pragma omp parallel
    {
      op.extract_even(even_solution,solution);
      op.preProp(even_source,odd_source,source);
      // The odd component is fixed by the source and the current even field.
      op.postProp(solution,even_solution,odd_source);
    }
  }

  bool solve(Field& solution, const Field& source, int& multiplications, double& relative, bool restored) {
    bool converged=false;
    int schur_multiplications=0;
#pragma omp parallel
    {
      if (!restored) {
        even_solution.set(0.0);
        op.preProp(even_source,odd_source,source);
      }
#pragma omp barrier
      const bool ok=solver.try_solve(even_solution,even_source,schur_multiplications,relative);
#pragma omp master
      { converged=ok; }
#pragma omp barrier
      op.postProp(solution,even_solution,odd_source);
    }
    multiplications=2*schur_multiplications;
    return converged;
  }
};

static int run(SolveOptions o,int rank,int ranks) {
  const std::string requested_method=o.method;
  const std::vector<int> scalar_packing(4,1);
#ifdef CLOVER_HAS_QXS
  std::vector<int> mg_packing{2*VLENXS,2*VLENYS,4,4};
  if (!o.sap_block.empty()) {
    const int lanes[] = {VLENXS, VLENYS, 1, 1};
    for (int axis=0;axis<4;++axis) {
      const int n=o.sap_block[axis];
      if (n%lanes[axis] || n>std::numeric_limits<int>::max()/lanes[axis])
        throw std::runtime_error("--sap-block is incompatible with QXS vector packing");
      mg_packing[axis]=n*lanes[axis];
    }
  }
  if(o.method=="mg" || o.method=="auto") {
    const auto mg_grid=o.grid.empty() ? clover_layout::choose(o.lattice,ranks,mg_packing) : o.grid;
    if(!mg_grid.empty() && clover_layout::compatible(o.lattice,mg_grid,mg_packing)) {
      o.method="mg";o.grid=mg_grid;
    } else if(o.method=="mg") {
      throw std::runtime_error("no compatible MG grid for this method, rank count and QXS packing; select --method auto or bicgstab");
    }
  }
#endif
  if(o.method=="auto") o.method="bicgstab";
  if((o.profile_mg || o.mg_max_iters || !o.sap_block.empty()) && o.method!="mg")
    throw std::runtime_error("--profile-mg, --sap-block and --mg-max-iters require a multigrid solve");
  if(o.precision=="single" && o.method!="mg")
    throw std::runtime_error("--precision single requires a multigrid solve");
  if(o.checkpoint && o.method!="bicgstab" && o.method!="eo-bicgstab")
    throw std::runtime_error("--checkpoint and --resume-from require bicgstab or eo-bicgstab");
  if (o.method=="mg" && (o.bicg_iterations || o.bicg_cycles))
    throw std::runtime_error("--bicg-iters and --bicg-cycles require a BiCGStab method");
  const int iterations=o.bicg_iterations ? o.bicg_iterations : (o.method=="eo-bicgstab" ? 12000 : 1000);
  const int cycles=o.bicg_cycles ? o.bicg_cycles : (o.method=="eo-bicgstab" ? 2 : 10);
  // The EO diagnostic counts both half-lattice applications per Schur multiplication.
  const long long maximum_count=std::numeric_limits<int>::max()/(o.method=="eo-bicgstab" ? 2 : 1);
  if (2LL*iterations*cycles+1>maximum_count)
    throw std::runtime_error("BiCGStab iteration budget exceeds the solver counter range");
  if(o.grid.empty()) o.grid=clover_layout::choose(o.lattice,ranks,(o.method=="eo-bicgstab") ? std::vector<int>{2,1,1,1} : scalar_packing);
  if(o.grid.empty()) throw std::runtime_error("no MPI grid divides this lattice for the requested rank count");
  validate_geometry(o,ranks);
  std::vector<int> local(4);
  for(int d=0;d<4;++d) local[d]=o.lattice[d]/o.grid[d];
  if(!rank) std::cout<<"LAYOUT grid="<<solve_join(o.grid)<<" local="<<solve_join(local)
                    <<" ranks="<<ranks<<" threads_per_rank="<<o.threads<<" method="<<o.method<<std::endl;
  if(o.plan_layout) return 0;
  if(!o.resume_from.empty()) {
    bool available=true;
    if(!rank) for(int rhs:o.rhs)
      available=available && std::ifstream(clover_checkpoint::prefix(o.resume_from,rhs)+".ready").good();
    require_all(available,"a requested source has no committed checkpoint");
  }
  if (o.gauge!="unit") {
    std::ifstream input(o.gauge.c_str(),std::ios::binary);
    require_all(bool(input),"cannot open ILDG gauge file on every rank");
    if (!rank) validate_lime(o.gauge,o.lattice);
  }
  require_all(rank || ::mkdir(o.output.c_str(),0770)==0,"output must be a new directory with an existing writable parent");
  Parameters main;
  main.set_int_vector("lattice_size",o.lattice);main.set_int_vector("grid_size",o.grid);
  main.set_int("number_of_thread",o.threads);main.set_int("number_of_color",3);
  main.set_string("log_filename","stdout");main.set_string("ildg_log_filename","stdout");
  main.set_string("verbose_level","General");bridge_setup(main);
  const int volume=CommonParameters::Nvol(),nd=CommonParameters::Ndim(),nt=CommonParameters::Nt(),lt=o.lattice[3];
  Parameters pp,sp;
  pp.set_string("verbose_level","Crucial");
  sp.set_string("verbose_level","Crucial");sp.set_double("rho_uniform",o.rho);
  Field_G smeared(volume,nd);
  double raw,stout;
  std::string gauge_hash,coordinates;
  {
    Field_G U(volume,nd);
    if (o.gauge=="unit") U.set_unit(); else GaugeConfig("ILDG_Parallel").read(U,o.gauge);
    if(o.checkpoint) {
      gauge_hash=clover_checkpoint::sha256(U.ptr(0),std::size_t(U.size())*sizeof(double));
      std::vector<int> site(4);
      for(int d=0;d<4;++d) site[d]=Communicator::ipe(d);
      coordinates=solve_join(site);
    }
    Projection_Stout_SU3 proj(pp);Smear_APE smear(&proj,sp);smear.smear(smeared,U);
    Staple_lex staple;
    raw=staple.plaquette(U);
    stout=staple.plaquette(smeared);
  }
  const double kappa=1.0/(2.0*(4.0+o.mass));
  Parameters op;
  op.set_string("fermion_type","Clover");op.set_string("gamma_matrix_type","Dirac");
  op.set_double("hopping_parameter",kappa);op.set_double("clover_coefficient",o.csw);
  op.set_int_vector("boundary_condition",{1,1,1,-1});op.set_string("verbose_level","Crucial");
  std::unique_ptr<Fopr> scalar;
  auto make_full_operator=[&]() {
    scalar.reset(Fopr::New("Clover",op));
    if (!scalar) throw std::runtime_error("Bridge++ Clover factory is unavailable");
    scalar->set_config(&smeared);scalar->set_mode("D");
  };
  if (!o.compact_eo) make_full_operator();
  const int initial_volume=o.compact_eo ? 0 : volume;
  Field_F b(initial_volume,1),x(initial_volume,1),res(initial_volume,1);
  double setup_seconds=0;
  std::string operator_check="NA";
  std::string sap_block="NA";
  std::unique_ptr<Solver_BiCGStab_Cmplx> bicg;
  std::unique_ptr<EvenOddBiCGStab> eo;
  Parameters eo_parameters=op;
  if(o.compact_eo) eo_parameters.set_bool("compact_eo_inverse",true);
  auto make_eo_operator=[&]() {
    Timer setup("even_odd_setup");setup.start();
    eo.reset(new EvenOddBiCGStab(eo_parameters,smeared,o.tolerance,iterations,cycles));
    setup.stop();return setup.elapsed_sec();
  };
#ifdef CLOVER_HAS_QXS
  std::unique_ptr<QxsSystem> qsys;
  if (o.method=="mg") {
    if (o.precision=="single")
      qsys.reset(new QxsSystemImpl<FloatField>(o,op,smeared,*scalar,b,x,res));
    else
      qsys.reset(new QxsSystemImpl<DoubleField>(o,op,smeared,*scalar,b,x,res));
    setup_seconds=qsys->setup_seconds;
    sap_block=solve_join(qsys->block);
    std::ostringstream value;value<<std::setprecision(17)<<qsys->operator_error;operator_check=value.str();
  } else
#endif
  if (o.method=="eo-bicgstab") {
    if(!o.compact_eo) setup_seconds=make_eo_operator();
  } else {
    bicg.reset(new Solver_BiCGStab_Cmplx(scalar.get()));
    bicg->set_parameters(iterations,cycles,o.tolerance*o.tolerance,!o.resume_from.empty(),0.7);
  }
  Parameters source_params;source_params.set_int_vector("source_position",o.position);source_params.set_string("verbose_level","Crucial");
  std::unique_ptr<Source> src;
  if (o.source=="point") src.reset(new Source_Local(source_params)); else src.reset(new Source_Wall(source_params));
  const std::string dir=o.output+"/";
  std::ofstream meta,columns,times,field;
  if (!rank) {
    meta=open_output(dir+"run.tsv");columns=open_output(dir+"columns.tsv");times=open_output(dir+"timeslices.tsv");
    meta<<"key\tvalue\n"<<"gauge\t"<<o.gauge<<"\nlattice\t"<<solve_join(o.lattice)<<"\ngrid\t"<<solve_join(o.grid)
        <<"\nposition\t"<<solve_join(o.position)<<"\nsource\t"<<o.source<<"\nrhs\t"<<solve_join(o.rhs)
        <<"\nmass\t"<<o.mass<<"\nkappa\t"<<kappa<<"\ncsw\t"<<o.csw<<"\nrho\t"<<o.rho
        <<"\nrequested_relative_tolerance\t"<<o.tolerance<<"\nindependent_acceptance_strict\t"<<o.acceptance
        <<"\nscalar_qxs_relative_l2\t"<<operator_check<<"\nraw_plaquette\t"<<raw<<"\nstout_plaquette\t"<<stout
        <<"\nprecision\t"<<o.precision<<"\nrequested_method\t"<<requested_method<<"\nsolver\t"<<o.method<<"\nsap_block\t"<<sap_block<<"\nthreads_per_rank\t"<<o.threads<<"\nnvec\t"<<o.vectors<<"\nsetup_steps\t"<<o.setup_steps
        <<"\nsmoother_iters\t"<<o.smooth_steps<<"\nmg_max_iters\t"<<(o.mg_max_iters ? o.mg_max_iters : 5000)
        <<"\nsetup_seconds\t"<<setup_seconds<<"\nfields\t"<<o.fields
        <<"\nbasis\tBridge Dirac\nnormalization\t2*kappa\n";
    if (o.compact_eo) meta<<"compact_eo\t1\n";
    if (bicg || eo || o.compact_eo)
      meta<<"initial_guess\t"<<(!o.resume_from.empty() ? "checkpoint" : o.method=="eo-bicgstab" ? "zero" : "source")<<"\nbicg_max_iterations_per_cycle\t"
          <<iterations<<"\nbicg_max_cycles\t"<<cycles<<'\n';
    if(o.checkpoint) meta<<"checkpoint_format\tclover-bicgstab-checkpoint-v1\n";
    if(!o.resume_from.empty()) meta<<"resume_from\t"<<o.resume_from<<'\n';
    if (o.method=="eo-bicgstab")
      meta<<"schur_max_iterations_per_cycle\t"<<iterations
          <<"\nschur_max_cycles\t"<<cycles<<"\nschur_omega_tolerance\t0\n";
    columns<<"rhs\tspin\tcolor\tsource_norm2\tnconv\touter_loops\tsolver_relative\tindependent_relative\tsolve_seconds\taccepted\n";
    times<<"rhs\tt\tdensity\n";
    meta.flush();
    columns.flush();
    times.flush();
  }
  require_all(rank || (meta && columns && times),"cannot create result tables");
  if (o.fields=="full") {
    std::ostringstream name;name<<dir<<"field";
    for (int mu=0;mu<4;++mu) name<<'_'<<Communicator::ipe(mu);
    name<<".f64";
    field.open(name.str().c_str(),std::ios::binary);
    require_all(bool(field),"cannot create all rank-local field shards on the output filesystem");
  }
  std::vector<double> summed(lt,0.0);
  bool passed=true;
  auto set_source=[&](int rhs) {
    src->set(b,rhs);
    if (o.source=="wall") scal(b,double(o.lattice[0])*o.lattice[1]*o.lattice[2]);
  };
  for (int rhs:o.rhs) {
    double compact_setup=0.0;
    if (o.compact_eo) {
      if (o.resume_from.empty()) compact_setup+=make_eo_operator();
      else { make_full_operator();res.reset(volume,1); }
      b.reset(volume,1);x.reset(volume,1);
    }
    set_source(rhs);
    const double source_norm=b.norm2();
    if(!o.resume_from.empty()) {
      const double saved=clover_checkpoint::read(o,rhs,rank,ranks,coordinates,gauge_hash,x.ptr(0),x.size());
      scalar->mult(res,x);axpy(res,-1.0,b);
      const double initial=std::sqrt(res.norm2()/source_norm);
      // Grouping floating-point reductions differently may change the recomputed norm.
      const double allowance=128*std::numeric_limits<double>::epsilon() *
                             (1+std::sqrt(x.norm2()/source_norm)+saved);
      require_all(std::isfinite(initial) && std::isfinite(allowance) && std::abs(initial-saved)<=allowance,
                  "checkpoint residual differs from independently reapplying this full Dirac operator");
      if(!rank) {
        meta<<"resume_initial_relative_rhs_"<<rhs<<'\t'<<initial<<'\n';meta.flush();
        std::cout<<"RESUME_COLUMN "<<rhs<<" independent_relative="<<std::setprecision(17)<<initial<<std::endl;
      }
      if(o.compact_eo) {
        scalar.reset();res.reset(0,0);b.reset(0,0);
        compact_setup+=make_eo_operator();
        b.reset(volume,1);set_source(rhs);
      }
      if(eo) {
        eo->restore(x,b);
        if(o.compact_eo) {
          eo.reset();make_full_operator();res.reset(volume,1);
        }
        scalar->mult(res,x);axpy(res,-1.0,b);
        const double projected=std::sqrt(res.norm2()/source_norm);
        require_all(std::isfinite(projected),"non-finite full residual after reconstructing the EO initial guess");
        if(!rank) {meta<<"resume_eo_projected_relative_rhs_"<<rhs<<'\t'<<projected<<'\n';meta.flush();}
        if(o.compact_eo) {
          // Rebuild the Schur system after the full operator verified the projected field.
          scalar.reset();res.reset(0,0);b.reset(0,0);
          compact_setup+=make_eo_operator();
          b.reset(volume,1);set_source(rhs);
          eo->load_even(x,b);
        }
      }
    }
    if(o.compact_eo && !rank) {
      meta<<"compact_eo_setup_seconds_rhs_"<<rhs<<'\t'<<compact_setup<<'\n';meta.flush();
    }
    int nconv=0;double solver_error=0,seconds=0;
    bool solver_converged=true;
#ifdef CLOVER_HAS_QXS
    if (qsys) seconds=qsys->solve(x,b,nconv,solver_error,o.profile_mg,rhs);
    else
#endif
    {
      Timer solve("source_column");solve.start();
      if (eo) solver_converged=eo->solve(x,b,nconv,solver_error,!o.resume_from.empty());
      else {
#pragma omp parallel
        {
          const bool ok=bicg->try_solve(x,b,nconv,solver_error);
#pragma omp master
          { solver_converged=ok; }
        }
      }
      solve.stop();seconds=solve.elapsed_sec();
    }
    if(o.compact_eo) {
      eo.reset();make_full_operator();res.reset(volume,1);
    }
    scalar->mult(res,x);axpy(res,-1.0,b);
    const double true_error=std::sqrt(res.norm2()/source_norm);
    const bool accept=solver_converged && std::isfinite(true_error) && true_error<o.acceptance;
    passed=passed && accept;
    if(o.checkpoint && !accept && std::isfinite(true_error))
      clover_checkpoint::write(o,rhs,rank,ranks,coordinates,gauge_hash,x.ptr(0),x.size(),true_error);
    scal(x,2.0*kappa);
    if (o.fields=="full") {
      field.write(reinterpret_cast<const char*>(x.ptr(0)),std::streamsize(x.size())*sizeof(double));field.flush();
      require_all(bool(field),"error writing full field shards");
    }
    std::vector<double> local(lt,0.0),global(lt,0.0);
    const int slice_doubles=(volume/nt)*x.nin(),t0=Communicator::ipe(3)*nt;
#pragma omp parallel for schedule(static)
    for (int t=0;t<nt;++t) {
      const double* values=x.ptr(t*slice_doubles);double density=0;
      for (int i=0;i<slice_doubles;++i) density+=values[i]*values[i];
      local[t0+t]=density;
    }
    Communicator::reduce_sum(lt,global.data(),local.data());
    if (!rank) {
      const int outer_loops=o.method=="eo-bicgstab" ? (nconv/2-1)/2 : nconv/2;
      columns<<rhs<<'\t'<<rhs/3<<'\t'<<rhs%3<<'\t'<<source_norm<<'\t'<<nconv<<'\t'<<outer_loops<<'\t'<<solver_error<<'\t'<<true_error<<'\t'<<seconds<<'\t'<<accept<<'\n';columns.flush();
      for (int t=0;t<lt;++t) {times<<rhs<<'\t'<<t<<'\t'<<global[t]<<'\n';summed[t]+=global[t];}
      times.flush();
      std::cout<<"COLUMN "<<rhs<<" independent_relative="<<std::setprecision(17)<<true_error<<" accepted="<<accept<<std::endl;
    }
    if(o.compact_eo) {
      scalar.reset();b.reset(0,0);x.reset(0,0);res.reset(0,0);
    }
  }
  if (o.fields=="full") field.close();
  if (!rank) {
    for (int t=0;t<lt;++t) times<<"sum\t"<<t<<'\t'<<summed[t]<<'\n';
    meta.close();columns.close();times.close();
  }
  require_all((o.fields!="full" || bool(field)) && (rank || (meta && columns && times)),"error finalizing result files");
  if (!rank) {
    std::ofstream status=open_output(dir+(passed?"COMPLETE":"FAILED"));
    status<<o.rhs.size()<<" source columns\n";status.close();
    if (!status) throw std::runtime_error("cannot write result status");
  }
  return passed?0:1;
}
int main(int argc,char** argv) {
  if (argc==2 && std::string(argv[1])=="--help") {std::cout<<solve_usage();return 0;}
  SolveOptions options;
  try {options=solve_options(argc,argv);}
  catch (const std::exception& e) {std::cerr<<e.what()<<'\n'<<solve_usage();return 2;}
  bridge_initialize(&argc,&argv);
  int rank=0,ranks=1;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  int rc=0;
  try {rc=run(options,rank,ranks);}
  catch (const clover_checkpoint::CollectiveError& e) {
    if (rank==0) std::cerr<<"clover-solve: "<<e.what()<<std::endl;
    rc=2;
  }
  catch (const std::exception& e) {
    std::cerr<<"clover-solve rank "<<rank<<": "<<e.what()<<std::endl;MPI_Abort(MPI_COMM_WORLD,2);return 2;
  }
  bridge_finalize();
  return rc;
}
