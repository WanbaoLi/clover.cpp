/*!
        @file    action_F_Rational_SF.cpp
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-03-04 12:05:22 #$
        @version $LastChangedRevision: 2698 $
*/

#include "lib/Action/Fermion/action_F_Rational_SF.h"
#include "lib/Fopr/fopr_thread-inc.h"

const std::string Action_F_Rational_SF::class_name = "Action_F_Rational_SF";

//====================================================================
void Action_F_Rational_SF::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  //vout.increase_indent();


  //vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Action_F_Rational_SF::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  set_parameters(params);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Action_F_Rational_SF::tidyup()
{
  // do nothing.
}

//====================================================================
void Action_F_Rational_SF::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Rational_SF::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Rational_SF::set_config(Field *U)
{
   m_U = U;

#pragma omp parallel
 {
  m_fopr_langev->set_config(U);
  m_fopr_H->set_config(U);
  m_forceF_MD->set_config(U);
 }

}


//====================================================================
double Action_F_Rational_SF::langevin(RandomNumbers *rand)
{
  const int NinF     = m_fopr_langev->field_nin();
  const int NvolF    = m_fopr_langev->field_nvol();
  const int NexF     = m_fopr_langev->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

  double H_psf = 0.0;

#pragma omp parallel
 {
  m_fopr_langev->set_config(m_U);
  m_fopr_langev->mult(m_psf, xi);

  Field_SF::set_boundary_zero(xi);

  const double xi2 = xi.norm();

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) H_psf = xi2 * xi2;
 }

  vout.general(m_vl, "  H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
double Action_F_Rational_SF::calcH()
{
  const int NinF     = m_fopr_H->field_nin();
  const int NvolF    = m_fopr_H->field_nvol();
  const int NexF     = m_fopr_H->field_nex();
  const int size_psf = NinF * NvolF * NexF * CommonParameters::NPE();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());

  Field v1(NinF, NvolF, NexF);

  double H_psf = 0.0;

#pragma omp parallel
 {
  m_fopr_H->set_config(m_U);
  m_fopr_H->mult(v1, m_psf);

  const double H_psf1 = dot(v1, m_psf);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) H_psf = H_psf1;
 }

  vout.general(m_vl, "  H_Frational  = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof      = %18.8f\n", H_psf / size_psf);

  return H_psf;
}


//====================================================================
void Action_F_Rational_SF::force(Field& force)
{
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  Field_G force1(Nvol, Nex);

#pragma omp parallel
 {
  m_forceF_MD->set_config(m_U);
  m_forceF_MD->force_core(force1, m_psf);

  Field_SF::set_boundary_spatial_link_zero(force1);

  copy(force, force1);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);
  vout.general(m_vl, "  Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
 }

  vout.decrease_indent();

}


//============================================================END=====
