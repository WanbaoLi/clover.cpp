/*!
        @file    action_F_Standard_SF.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-03-04 12:05:22 #$
        @version $LastChangedRevision: 2698 $
*/

#include "action_F_Standard_SF.h"

const std::string Action_F_Standard_SF::class_name = "Action_F_Standard_SF";

//====================================================================
void Action_F_Standard_SF::init(const Parameters& params)
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

  set_parameters(params);

  const int         Niter       = 100;
  const int         Nrestart    = 40;
  const double      Stop_cond   = 1.0e-24;
  const std::string solver_type = "CG";

  m_solver = Solver::New(solver_type, m_fopr);
  m_solver->set_parameters(Niter, Nrestart, Stop_cond);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Action_F_Standard_SF::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  const int         Niter       = 100;
  const int         Nrestart    = 40;
  const double      Stop_cond   = 1.0e-24;
  const std::string solver_type = "CG";

  m_solver = Solver::New(solver_type, m_fopr);
  m_solver->set_parameters(Niter, Nrestart, Stop_cond);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Action_F_Standard_SF::tidyup()
{
  delete m_solver;
}


//====================================================================
void Action_F_Standard_SF::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Standard_SF::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Standard_SF::set_config(Field *U)
{
  m_U = U;

#pragma omp parallel
 {
  m_fopr->set_config(U);
   m_forceF->set_config(U);
 }

}

//====================================================================
double Action_F_Standard_SF::langevin(RandomNumbers *rand)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr->field_nin();
  const int NvolF    = m_fopr->field_nvol();
  const int NexF     = m_fopr->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

#pragma omp parallel
 {
  m_fopr->set_config(m_U);
  m_fopr->set_mode("Ddag");
  m_fopr->mult(m_psf, xi);
 }

  Field_SF::set_boundary_zero(xi);

  const double xi2   = xi.norm();
  const double H_psf = xi2 * xi2;

  vout.decrease_indent();

  vout.general(m_vl, "  H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Standard_SF::calcH()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  const int NinF     = m_fopr->field_nin();
  const int NvolF    = m_fopr->field_nvol();
  const int NexF     = m_fopr->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  Field v1(NinF, NvolF, NexF);
  double H_psf = 0.0;

#pragma omp parallel
 {
  m_fopr->set_config(m_U);
  m_fopr->set_mode("DdagD");

  int    Nconv;
  double diff;
  m_solver->solve(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "  Nconv = %d  diff  = %.8e\n", Nconv, diff);

  Field_SF::set_boundary_zero(v1);
  Field_SF::set_boundary_zero(m_psf);

  const double H_psf1 = dot(v1, m_psf);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) H_psf = H_psf1;
 }

  vout.decrease_indent();
  vout.general(m_vl, "  H_Fstandard_SF = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof        = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Standard_SF::force(Field& force)
{
  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  Field   eta(NinF, NvolF, NexF);
  Field_G force1(Nvol, Nex);

#pragma omp parallel
 {
  m_fopr->set_config(m_U);
  m_fopr->set_mode("DdagD");

  int    Nconv;
  double diff;
  m_solver->solve(eta, m_psf, Nconv, diff);

  vout.general(m_vl, "  Solver: Nconv = %6d  diff  = %12.6e\n",
               Nconv, diff);

  m_forceF->set_config(m_U);

  m_forceF->force_core(force1, eta);

  Field_SF::set_boundary_spatial_link_zero(force1);

  copy(force, force1);
#pragma omp barrier

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
 }

  vout.decrease_indent();

}


//============================================================END=====
