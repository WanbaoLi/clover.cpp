/*!
        @file    smear_APE.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Smear/smear_APE.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE::register_factory();
}
#endif

const std::string Smear_APE::class_name = "Smear_APE";

//====================================================================
void Smear_APE::init(const Parameters& params)
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

  m_rho.resize(m_Ndim * m_Ndim);
  for(int i = 0; i < m_Ndim * m_Ndim; ++i){
    m_rho[i] = 0.0;
  }

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
void Smear_APE::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n", class_name.c_str());
  vout.increase_indent();

  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  m_rho.resize(m_Ndim * m_Ndim);
  for(int i = 0; i < m_Ndim * m_Ndim; ++i){
    m_rho[i] = 0.0;
  }

  m_staple = new Staple_lex();

  m_ct  = new Field_G(m_Nvol, 1);
  m_ut1 = new Field_G(m_Nvol, 1);
  m_ut2 = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Smear_APE::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_staple;

  delete m_ct;
  delete m_ut1;
  delete m_ut2;
}


//====================================================================
void Smear_APE::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double rho1;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1);

}


//====================================================================
void Smear_APE::set_parameters(const double rho1)
{
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    for (int mu = 0; mu < m_Ndim; ++mu) {
      for (int nu = 0; nu < m_Ndim; ++nu) {
        m_rho[index(mu, nu)] = rho1;
      }
    }
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  rho (uniform) = %8.4f\n", rho1);

}


//====================================================================
void Smear_APE::set_parameters(const std::vector<double>& rho)
{
  assert(rho.size() == m_Ndim * m_Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    for (int mu = 0; mu < m_Ndim; ++mu) {
      for (int nu = 0; nu < m_Ndim; ++nu) {
        m_rho[index(mu, nu)] = rho[index(mu, nu)];
      }
    }
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      vout.general(m_vl, "  rho(%d,%d) = %8.4f\n",
                   mu, nu, rho[index(mu,nu)]);
    }
  }

}


//====================================================================
void Smear_APE::get_parameters(Parameters& params) const
{
  params.set_double_vector("rho", m_rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Smear_APE::smear(Field_G& Usmear, const Field_G& U)
{
#pragma omp barrier

  assert(U.nvol() == m_Nvol);
  assert(U.nex()  == m_Ndim);
  assert(Usmear.nvol() == m_Nvol);
  assert(Usmear.nex()  == m_Ndim);

  Usmear.set(0.0);
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {

    m_ct->set(0.0);
    copy(*m_ut1, 0, U, mu);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {

        double rho = m_rho[index(mu, nu)];
        m_staple->upper(*m_ut2, U, mu, nu);
        axpy(*m_ct, 0, rho, *m_ut2, 0);
#pragma omp barrier

        m_staple->lower(*m_ut2, U, mu, nu);
        axpy(*m_ct, 0, rho, *m_ut2, 0);
#pragma omp barrier
      }
    }

    double rho0 = m_rho[index(mu, mu)];
    m_proj->project(*m_ut2, rho0, *m_ct, *m_ut1);
    copy(Usmear, mu, *m_ut2, 0);
#pragma omp barrier
  }
}


//============================================================END=====
