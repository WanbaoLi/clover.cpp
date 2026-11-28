#pragma once
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
