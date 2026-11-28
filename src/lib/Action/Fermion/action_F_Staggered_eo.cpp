/*!
        @file    action_F_Staggered_eo.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-03-04 12:05:22 #$
        @version $LastChangedRevision: 2698 $
*/

#include "lib/Action/Fermion/action_F_Staggered_eo.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Field/field_F_1spinor.h"

const std::string Action_F_Staggered_eo::class_name = "Action_F_Staggered_eo";

//====================================================================
void Action_F_Staggered_eo::init()
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
void Action_F_Staggered_eo::tidyup()
{
  // do nothing
}


//====================================================================
void Action_F_Staggered_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

}


//====================================================================
void Action_F_Staggered_eo::get_parameters(Parameters& params) const
{
  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Staggered_eo::set_config(Field *U)
{
  m_U = (Field_G *)U;

#pragma omp parallel
 {
  m_fopr->set_config(U);
 }

}

//====================================================================
double Action_F_Staggered_eo::langevin(RandomNumbers *rand)
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Ndim = CommonParameters::Ndim();

  const int Nvol  = CommonParameters::Nvol();
  const int Nvol2 = Nvol / 2;
  const int NPE   = CommonParameters::NPE();

  Field_F_1spinor xi;

  rand->gauss_lex_global(xi);

  m_psf.reset(2 * Nc, Nvol2, 1);

  Field_F_1spinor xi_e(Nvol2, 1);
  Field_F_1spinor xi_o(Nvol2, 1);

#pragma omp parallel
 {
  m_fopr->set_config(m_U);

  m_index_eo.convertField(xi_e, xi, 0);
  m_index_eo.convertField(xi_o, xi, 1);

  m_fopr->mult_dag(m_psf, xi_e, "Dee");
  m_fopr->mult_dag(xi_e, xi_o, "Deo");

  axpy(m_psf, 1.0, xi_e);
 }

 const double H_psf = calcH_impl();

  const double Fpsf = 1.0/(double(2*Nc * Nvol/2) * double(NPE));

  vout.decrease_indent();
  vout.general(m_vl, "  H_F     = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof = %18.8f\n", H_psf * Fpsf);

  return H_psf;
}


//====================================================================
double Action_F_Staggered_eo::calcH()
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  const int Nc   = CommonParameters::Nc();
  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const double H_psf = calcH_impl();

  const double Fpsf = 1.0/(double(2*Nc * Nvol/2) * double(NPE));

  vout.decrease_indent();
  vout.general(m_vl, "  H_F     = %18.8f\n", H_psf);
  vout.general(m_vl, "  H_F/dof = %18.8f\n", H_psf * Fpsf);

  return H_psf;
}


//====================================================================
double Action_F_Staggered_eo::calcH_impl()
{
  const int Nvol = CommonParameters::Nvol();

  const int    Niter     = 100;
  const int    Nrestart  = 40;
  const double Stop_cond = 1.0e-24;

  m_fopr->set_mode("D");

  Field v1(m_psf);
  assert(v1.nvol() == Nvol / 2);

  m_fprop_H->set_config(m_U);

  int    Nconv;
  double diff;
  m_fprop_H->invert_D(v1, m_psf, Nconv, diff);
  vout.general(m_vl, "  FproH: Nconv = %d  diff  = %12.6e\n",
               Nconv, diff);

  double H_psf = 0.0;

#pragma omp parallel
 {
  const double H_psf1 = dot(v1, m_psf);

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) H_psf = H_psf1;
 }

  return H_psf;

}


//====================================================================
void Action_F_Staggered_eo::force(Field& force)
{
  vout.general(m_vl, "%s: %s\n", class_name.c_str(), m_label.c_str());
  vout.increase_indent();

  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  const int    Niter     = 100;
  const int    Nrestart  = 40;
  const double Stop_cond = 1.0e-24;

  Field_F_1spinor eta(Nvol / 2, 1);

  int    Nconv;
  double diff;
  m_fprop_MD->set_config(m_U);
  m_fprop_MD->invert_D(eta, m_psf, Nconv, diff);

  vout.general(m_vl, "  Fprop_MD: Nconv = %6d  diff  = %12.6e\n",
               Nconv, diff);

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
