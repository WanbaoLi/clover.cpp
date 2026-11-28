/*!
        @file    action_F_Rational.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2026-03-04 12:05:22 #$

        @version $LastChangedRevision: 2698 $
*/

#include "lib/Action/Fermion/action_F_Rational.h"
#include "lib/ResourceManager/threadManager.h"

const std::string Action_F_Rational::class_name = "Action_F_Rational";

//====================================================================
void Action_F_Rational::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Rational::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Rational::set_config(Field *U)
{
  m_U = U;

#pragma omp parallel
 {
  m_fopr_langev->set_config(U);
  m_fopr_H->set_config(U);
  m_fopr_force_MD->set_config(U);
 }

}

//====================================================================
double Action_F_Rational::langevin(RandomNumbers *rand)
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  m_timer.reset();
  m_timer.start();

  const int NinF     = m_fopr_langev->field_nin();
  const int NvolF    = m_fopr_langev->field_nvol();
  const int NexF     = m_fopr_langev->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  m_psf.reset(NinF, NvolF, NexF);

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  double H_psf = 0.0;

#pragma omp parallel
 {
  m_fopr_langev->set_config(m_U);

  m_fopr_langev->mult(m_psf, xi);

  const double H_psf1 = xi.norm2();

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();

  vout.decrease_indent();
  vout.general(m_vl, "  H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", elapsed_time);

  return H_psf;
}


//====================================================================
double Action_F_Rational::calcH()
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  m_timer.reset();
  m_timer.start();

  const int NinF     = m_fopr_H->field_nin();
  const int NvolF    = m_fopr_H->field_nvol();
  const int NexF     = m_fopr_H->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  Field v1(NinF, NvolF, NexF);

  double H_psf = 0.0;

#pragma omp parallel
 {
  m_fopr_H->set_config(m_U);
  m_fopr_H->mult(v1, m_psf);

  const double H_psf1 = dot(v1, m_psf);
  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();

  vout.decrease_indent();
  vout.general(m_vl, "  H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", elapsed_time);

  return H_psf;
}


//====================================================================
void Action_F_Rational::force(Field& force)
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  m_timer.reset();
  m_timer.start();

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

#pragma omp parallel
 {
  m_fopr_force_MD->set_config(m_U);

  m_fopr_force_MD->force_core(force, m_psf);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);

  vout.general(m_vl, "Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
 }
  vout.decrease_indent();

  m_timer.stop();
  double elapsed_time = m_timer.elapsed_sec();
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", elapsed_time);

}

//============================================================END=====
