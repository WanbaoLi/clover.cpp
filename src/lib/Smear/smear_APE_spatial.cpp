/*!
        @file    smear_APE_spatial.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Smear/smear_APE_spatial.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE_spatial::register_factory();
}
#endif

const std::string Smear_APE_spatial::class_name = "Smear_APE_spatial";

//====================================================================
void Smear_APE_spatial::init(const Parameters& params)
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

  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  set_parameters(params);
 
  m_staple = new Staple_lex();

  m_ct  = new Field_G(m_Nvol, 1);
  m_ut1 = new Field_G(m_Nvol, 1);
  m_ut2 = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Smear_APE_spatial::tidyup()
{
  delete m_staple;

  delete m_ct;
  delete m_ut1;
  delete m_ut2;
}

//====================================================================
void Smear_APE_spatial::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double rho;

  int err = 0;
  err += params.fetch_double("rho", rho);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho);
}


//====================================================================
void Smear_APE_spatial::get_parameters(Parameters& params) const
{
  params.set_double("rho", m_rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Smear_APE_spatial::set_parameters(const double rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  rho = %10.6F\n", rho);

  //- range check
  // NB. rho == 0 is allowed.

  //- store values
  m_rho = rho;
}


//====================================================================
void Smear_APE_spatial::smear(Field_G& Usmear, const Field_G& U)
{
  assert(U.nvol() == m_Nvol);
  assert(U.nex()  == m_Ndim);
  assert(Usmear.nvol() == m_Nvol);
  assert(Usmear.nex()  == m_Ndim);

  const int Ndim_spc = m_Ndim - 1;

  //  Staple_lex staple;

  double plaq = m_staple->plaq_s(U);
  vout.general(m_vl, "  plaq_s(org)   = %12.8f\n", plaq);

  plaq = m_staple->plaq_t(U);
  vout.general(m_vl, "  plaq_t(org)   = %12.8f\n", plaq);

  Usmear.set(0.0);

  for (int mu = 0; mu < Ndim_spc; ++mu) {
    //Field_G c_tmp;
    m_ct->set(0.0);

    //Field_G u_tmp;
    copy(*m_ut1, 0, U, mu);

    //Field_G u_tmp2;

    for (int nu = 0; nu < Ndim_spc; ++nu) {
      if (nu != mu) {
        m_staple->upper(*m_ut2, U, mu, nu);
        axpy(*m_ct, 0, m_rho, *m_ut2, 0);

        m_staple->lower(*m_ut2, U, mu, nu);
        axpy(*m_ct, 0, m_rho, *m_ut2, 0);
      }
    }

    m_proj->project(*m_ut2, m_rho, *m_ct, *m_ut1);
    copy(Usmear, mu, *m_ut2, 0);
  }

  const int mu = m_Ndim - 1; // temporal link: unsmeared.
  copy(Usmear, mu, U, mu);

  plaq = m_staple->plaq_s(Usmear);
  vout.general(m_vl, "  plaq_s(smear) = %12.8f\n", plaq);

  plaq = m_staple->plaq_t(Usmear);
  vout.general(m_vl, "  plaq_t(smear) = %12.8f\n", plaq);
}


//============================================================END=====
