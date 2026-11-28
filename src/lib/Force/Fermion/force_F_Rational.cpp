/*!
        @file    force_F_Rational.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Force/Fermion/force_F_Rational.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Tools/math_Rational.h"

const std::string Force_F_Rational::class_name = "Force_F_Rational";

//====================================================================
void Force_F_Rational::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  set_parameters_impl(params);

  setup(0);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Rational::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_Np = 0;

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Rational::setup(int Np_prev)
{
  ThreadManager::assert_single_thread(class_name);

  if(Np_prev != 0) tidyup();

  Parameters params_solver;
  params_solver.set_int("number_of_shifts", m_Np);
  params_solver.set_int("maximum_number_of_iteration", m_Niter);
  params_solver.set_double("convergence_criterion_squared", m_Stop_cond);
  params_solver.set_string("verbose_level", vout.get_verbose_level(m_vl));

  m_solver = new Shiftsolver_CG(m_fopr, params_solver);

  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  m_force1 = new Field_G(Nvol, Ndim);
  m_force2 = new Field_G(Nvol, Ndim);

  m_psi.resize(m_Np);

  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  for (int i = 0; i < m_Np; ++i) {
    m_psi[i].reset(NinF, NvolF, NexF);
  }

  m_eta = new Field(NinF, NvolF, NexF);

}

//====================================================================
void Force_F_Rational::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_solver;

  delete m_force1;
  delete m_force2;
  delete m_eta;

}


//====================================================================
void Force_F_Rational::set_parameters(const Parameters& params)
{

  int Np_prev = m_Np;

  set_parameters_impl(params);

  if(Np_prev != m_Np) setup(Np_prev);


}


//====================================================================
void Force_F_Rational::set_parameters_impl(const Parameters& params)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    std::string vlevel;
    if (!params.fetch_string("verbose_level", vlevel)) {
      m_vl = vout.set_verbose_level(vlevel);
    }
  }
#pragma omp barrier

  //- fetch and check input parameters
  int    Np, n_exp, d_exp;
  double x_min, x_max;
  int    Niter;
  double Stop_cond;

  int err = 0;
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_int("exponent_numerator", n_exp);
  err += params.fetch_int("exponent_denominator", d_exp);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(Np, n_exp, d_exp, x_min, x_max, Niter, Stop_cond);

}


//====================================================================
void Force_F_Rational::get_parameters(Parameters& params) const
{
  params.set_int("number_of_poles", m_Np);
  params.set_int("exponent_numerator", m_n_exp);
  params.set_int("exponent_denominator", m_d_exp);
  params.set_double("lower_bound", m_x_min);
  params.set_double("upper_bound", m_x_max);
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion_squared", m_Stop_cond);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Rational::set_parameters(
                            const int Np,
                            const int n_exp, const int d_exp,
                            const double x_min, const double x_max,
                            const int Niter, const double Stop_cond)
{
  int Np_prev = m_Np;

  set_parameters(Np, n_exp, d_exp, x_min, x_max, Niter, Stop_cond);

  if(Np_prev != m_Np) setup(Np_prev);

}


//====================================================================
void Force_F_Rational::set_parameters_impl(
                            const int Np,
                            const int n_exp, const int d_exp,
                            const double x_min, const double x_max,
                            const int Niter, const double Stop_cond)
{
#pragma omp barrier

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(Np);
  err += ParameterCheck::non_zero(n_exp);
  err += ParameterCheck::non_zero(d_exp);
                  // NB. x_min, x_max = 0.0 is allowed.
  err += ParameterCheck::non_zero(Niter);
  err += ParameterCheck::square_non_zero(Stop_cond);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_Np        = Np;
    m_n_exp     = n_exp;
    m_d_exp     = d_exp;
    m_x_min     = x_min;
    m_x_max     = x_max;
    m_Niter     = Niter;
    m_Stop_cond = Stop_cond;

    m_cl.resize(m_Np);
    m_bl.resize(m_Np);

    //- Rational approximation
    const double x_min2 = m_x_min * m_x_min;
    const double x_max2 = m_x_max * m_x_max;

    Math_Rational rational;
    rational.set_parameters(m_Np, m_n_exp, m_d_exp, x_min2, x_max2);
    rational.get_parameters(m_a0, m_bl, m_cl);
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: parameters:\n", class_name.c_str());
  vout.general(m_vl, "  Np        = %d\n", m_Np);
  vout.general(m_vl, "  n_exp     = %d\n", m_n_exp);
  vout.general(m_vl, "  d_exp     = %d\n", m_d_exp);
  vout.general(m_vl, "  x_min     = %12.8f\n", m_x_min);
  vout.general(m_vl, "  x_max     = %12.8f\n", m_x_max);
  vout.general(m_vl, "  Niter     = %d\n", m_Niter);
  vout.general(m_vl, "  Stop_cond = %8.2e\n", m_Stop_cond);
  vout.general(m_vl, "  a0 = %18.14f\n", m_a0);
  for (int i = 0; i < m_Np; i++) {
    vout.general(m_vl, "   bl[%d] = %18.14f  cl[%d] = %18.14f\n",
                 i, m_bl[i], i, m_cl[i]);
  }

#pragma omp barrier
}


//====================================================================
void Force_F_Rational::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr->set_config(U);
  m_force->set_config(U);

}

//====================================================================
void Force_F_Rational::force_udiv(Field& force_, const Field& eta_)
{
#pragma omp barrier

  copy(*m_eta, eta_);
#pragma omp barrier

  force_udiv_impl(*m_force2, *m_eta);

  copy(force_, *m_force2);
#pragma omp barrier

}


//====================================================================
void Force_F_Rational::force_udiv_impl(Field_G& force, const Field& eta)
{
#pragma omp barrier

  vout.general(m_vl, "  Shift solver in force calculation\n");
  vout.general(m_vl, "    Number of shift values = %d\n", m_cl.size());

  m_fopr->set_mode("DdagD");

  int            Nconv;
  double         diff;
  m_solver->solve(m_psi, m_cl, eta, Nconv, diff);
  vout.general(m_vl, "    diff(max) = %22.15e  \n", diff);

  force.set(0.0);
#pragma omp barrier

  for (int i = 0; i < m_Np; ++i) {
    m_force->force_udiv(*m_force1, m_psi[i]);
    scal(*m_force1, m_bl[i]);
    axpy(force, 1.0, *m_force1);
#pragma omp barrier
  }

}


//====================================================================
void Force_F_Rational::force_core1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_Rational::force_udiv1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//===========================================================END======
