/*!
        @file    action_F_Standard_lex.cpp
        @brief
        @author  UEDA, Satoru (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-03-04 12:05:22 #$
        @version $LastChangedRevision: 2698 $
*/

#include "lib/Action/Fermion/action_F_Standard_lex.h"
#include "lib/ResourceManager/threadManager.h"

const std::string Action_F_Standard_lex::class_name
                                          = "Action_F_Standard_lex";

//====================================================================
void Action_F_Standard_lex::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Standard_lex::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Standard_lex::set_config(Field *U)
{
  m_U = U;

#pragma omp parallel
 {
  m_fopr->set_config(U);
  m_fopr_force->set_config(U);
 }

}


//====================================================================
double Action_F_Standard_lex::langevin(RandomNumbers *rand)
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  m_timer.reset();
  m_timer.start();

  const int Nvol = CommonParameters::Nvol();

  const int NinF     = m_fopr->field_nin();
  const int NvolF    = m_fopr->field_nvol();
  const int NexF     = m_fopr->field_nex();
  const int NPE      = CommonParameters::NPE();
  double Fpsf = 1.0/(double(NinF * NvolF * NexF) * double(NPE));

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  double H_psf = 0.0;

#pragma omp parallel
 {
  m_fopr->set_config(m_U);
  m_fopr->set_mode("Ddag");

  m_fopr->mult(m_psf, xi);

  const double xi2   = xi.norm();
  const double H_psf1 = xi2 * xi2;

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();

  vout.decrease_indent();
  vout.general(m_vl, "  H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf * Fpsf);
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", elapsed_time);

  return H_psf;
}


//====================================================================
double Action_F_Standard_lex::calcH()
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  m_timer.reset();
  m_timer.start();

  const int NinF     = m_fopr->field_nin();
  const int NvolF    = m_fopr->field_nvol();
  const int NexF     = m_fopr->field_nex();
  const int NPE      = CommonParameters::NPE();
  double Fpsf = 1.0/(double(NinF * NvolF * NexF) * double(NPE));

  Field  v1(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_H->set_config(m_U);
  m_fprop_H->invert_DdagD(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "  Fprop_H: Nconv = %6d  diff  = %12.6e\n",
               Nconv, diff);

  double H_psf = 0.0;

#pragma omp parallel
 {  
  const double H_psf1 = dot(v1, m_psf);
  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();

  vout.decrease_indent();
  vout.general(m_vl, "  H_Fstandard  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf * Fpsf);
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", elapsed_time);

  return H_psf;
}


//====================================================================
void Action_F_Standard_lex::force(Field& force)
{
  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  m_timer.reset();
  m_timer.start();

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  const int NinF  = m_fopr->field_nin();
  const int NvolF = m_fopr->field_nvol();
  const int NexF  = m_fopr->field_nex();

  Field  eta(NinF, NvolF, NexF);
  int    Nconv;
  double diff;
  m_fprop_MD->set_config(m_U);
  m_fprop_MD->invert_DdagD(eta, m_psf, Nconv, diff);
  vout.general(m_vl, "  Fprop_MD: Nconv = %6d  diff  = %12.6e\n",
               Nconv, diff);

#pragma omp parallel
 {
  m_fopr_force->set_config(m_U);
  m_fopr_force->force_core(force, eta);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);

  vout.general(m_vl,"Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
 }

  vout.decrease_indent();

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", elapsed_time);

}


//============================================================END=====
