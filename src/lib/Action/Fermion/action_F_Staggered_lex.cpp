/*!
        @file    action_F_Staggered_lex.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-03-04 12:05:22 #$
        @version $LastChangedRevision: 2698 $
*/

#include "Action/Fermion/action_F_Staggered_lex.h"
#include "lib/ResourceManager/threadManager.h"

const std::string Action_F_Staggered_lex::class_name
                                      = "Action_F_Staggered_lex";

//====================================================================
void Action_F_Staggered_lex::init()
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
void Action_F_Staggered_lex::tidyup()
{
  // do nothing.
}

//====================================================================
void Action_F_Staggered_lex::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }
}


//====================================================================
void Action_F_Staggered_lex::set_parameters()
{
  m_vl = CommonParameters::Vlevel();

  int Nc   = CommonParameters::Nc();
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();
  int NinG = 2 * Nc * Nc;

  vout.detailed(m_vl, "%s:\n", class_name.c_str());
}


//====================================================================
void Action_F_Staggered_lex::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Staggered_lex::set_config(Field *U)
{
  m_U = U;
#pragma omp parallel
 {
  m_fopr->set_config(U);
  m_forceF->set_config(U);
 }

}


//====================================================================
double Action_F_Staggered_lex::langevin(RandomNumbers *rand)
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();

  assert(NvolF == Nvol);
  m_psf.reset(NinF, NvolF, NexF);

  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());

  Field xi(NinF, NvolF, NexF);
  rand->gauss_lex_global(xi);

#pragma omp parallel
 {
  m_fopr->set_config(m_U);
  m_fopr->set_mode("Ddag");

  m_fopr->mult(m_psf, xi);

  // pseudo-fermion field on odd sites are discarded.
  m_fopr->mult_gm5(xi, m_psf);
  axpy(m_psf, 1.0, xi);
  scal(m_psf, 0.5);
 }

 // Hamiltonian must be determined by calling calc_H.
  double H_psf = calcH_impl();

  double Fpsf = 1.0/(double(NinF * NvolF/2 * NexF));
  Fpsf = Fpsf/double(CommonParameters::NPE());

  vout.decrease_indent();
  vout.general(m_vl, "  H_F     = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof = %18.8f\n", H_psf * Fpsf);

  return H_psf;
}


//====================================================================
double Action_F_Staggered_lex::calcH()
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();
  double Fpsf = 1.0/(double(NinF * NvolF/2 * NexF));
  Fpsf = Fpsf/double(CommonParameters::NPE());

  double H_psf = calcH_impl();

  vout.decrease_indent();
  vout.general(m_vl, "  H_F     = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof = %18.8f\n", H_psf * Fpsf);

  return H_psf;
}


//====================================================================
double Action_F_Staggered_lex::calcH_impl()
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF     = m_fopr->field_nin();
  int NvolF    = m_fopr->field_nvol();
  int NexF     = m_fopr->field_nex();

  Field v1(NinF, NvolF, NexF);

  m_fprop_H->set_config(m_U);

  int    Nconv;
  double diff;
  m_fprop_H->invert_DdagD(v1, m_psf, Nconv, diff);
  vout.general(m_vl, "  Fprop_MF: Nconv = %6d  diff  = %12.6e\n",
               Nconv, diff);

  double H_psf = 0.0;

#pragma omp parallel
 {
  double H_psf1 = dot(v1, m_psf);
 
  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  return H_psf;
}


//====================================================================
void Action_F_Staggered_lex::force(Field& force)
{
  vout.general(m_vl, "  %s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  int Nin  = m_U->nin();
  int Nvol = m_U->nvol();
  int Nex  = m_U->nex();
  int Nc   = CommonParameters::Nc();
  int Ndim = CommonParameters::Ndim();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  int   NinF  = m_fopr->field_nin();
  int   NvolF = m_fopr->field_nvol();
  int   NexF  = m_fopr->field_nex();
  Field eta(NinF, NvolF, NexF);


  int    Nconv;
  double diff;

  m_fprop_MD->set_config(m_U);
  m_fprop_MD->invert_DdagD(eta, m_psf, Nconv, diff);
  vout.general(m_vl, "  Fprop_MD:  %6d  %18.15e\n", Nconv, diff);

#pragma omp parallel
 {
  m_forceF->set_config(m_U);

  m_forceF->force_core(force, eta);

  double Fave, Fmax, Fdev;
  force.stat(Fave, Fmax, Fdev);

  vout.general(m_vl, "Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
               Fave, Fmax, Fdev);
 }

  vout.decrease_indent();

}

//============================================================END=====
