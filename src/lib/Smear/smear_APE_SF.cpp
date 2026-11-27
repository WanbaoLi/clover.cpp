/*!
        @file    smear_APE_SF.cpp
        @brief
        @author  Yusuke Taniguchi (maintained by H. Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Smear/smear_APE_SF.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_APE_SF::register_factory();
}
#endif

const std::string Smear_APE_SF::class_name = "Smear_APE_SF";

//====================================================================
void Smear_APE_SF::init(const Parameters& params)
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
  const int Nc = CommonParameters::Nc();

  m_rho.resize(m_Ndim * m_Ndim);
  for(int i = 0; i < m_Ndim * m_Ndim; ++i){
    m_rho[i] = 0.0;
  }

  set_parameters(params);

  m_staple = new Staple_SF();
  m_staple->set_parameters(m_phi, m_phipr);


  m_ct  = new Field_G(m_Nvol, 1);
  m_ut1 = new Field_G(m_Nvol, 1);
  m_ut2 = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Smear_APE_SF::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_staple;

  delete m_ct;
  delete m_ut1;
  delete m_ut2;

}

//====================================================================
void Smear_APE_SF::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double              rho1;
  std::vector<double> phi, phipr;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho1);
  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1, phi, phipr);

}


//====================================================================
void Smear_APE_SF::get_parameters(Parameters& params) const
{
  std::vector<double> rho(m_rho.size());
  for (size_t i = 0; i < m_rho.size(); ++i) {
    rho[i] = m_rho[i];
  }
  params.set_double_vector("rho", rho);
  params.set_double_vector("phi", m_phi);
  params.set_double_vector("phipr", m_phipr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Smear_APE_SF::set_parameters(const double rho1,
                                  const std::vector<double>& phi,
                                  const std::vector<double>& phipr)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    for (int mu = 0; mu < m_Ndim; ++mu) {
      for (int nu = 0; nu < m_Ndim; ++nu) {
        m_rho[index(mu, nu)] = rho1;
      }
    }

    m_phi.resize(3);
    m_phipr.resize(3);
    for (int i = 0; i < 3; ++i) {
      m_phi[i]   = phi[i];
      m_phipr[i] = phipr[i];
    }

    Field_SF::set_boundary_matrix(m_wk, m_phi);
    Field_SF::set_boundary_matrix(m_wkpr, m_phipr);
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  rho = %8.4f\n", m_rho[0]);

  vout.general(m_vl, "  phi1  = %12.6f\n", m_phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", m_phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", m_phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", m_phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", m_phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", m_phipr[2]);

}


//====================================================================
void Smear_APE_SF::set_parameters(const std::vector<double>& rho,
                                  const std::vector<double>& phi,
                                  const std::vector<double>& phipr)
{
#pragma omp barrier

  assert(rho.size() == m_Ndim * m_Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    for (int mu = 0; mu < m_Ndim; ++mu) {
      for (int nu = 0; nu < m_Ndim; ++nu) {
        m_rho[index(mu, nu)] = rho[index(mu, nu)];
      }
    }

    m_phi.resize(3);
    m_phipr.resize(3);
    for (int i = 0; i < 3; ++i) {
      m_phi[i]   = phi[i];
      m_phipr[i] = phipr[i];
    }
    Field_SF::set_boundary_matrix(m_wk, m_phi);
    Field_SF::set_boundary_matrix(m_wkpr, m_phipr);
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      vout.general(m_vl, "  rho(%d,%d) = %8.4f\n",
                   mu, nu, m_rho[index(mu, nu)]);
    }
  }
  vout.general(m_vl, "  phi1  = %12.6f\n", m_phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", m_phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", m_phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", m_phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", m_phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", m_phipr[2]);

}


//====================================================================
void Smear_APE_SF::smear(Field_G& Usmear, const Field_G& U)
{
  const int Nvol = CommonParameters::Nvol();

  assert(U.nvol() == Nvol);
  assert(U.nex() == m_Ndim);
  assert(Usmear.nvol() == Nvol);
  assert(Usmear.nex() == m_Ndim);

  Usmear.set(0.0);
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {

    m_ct->set(0.0);
    copy(*m_ut1, 0, U, mu);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {
        double rho = m_rho[mu + m_Ndim * nu];
        m_staple->upper(*m_ut2, U, mu, nu);
        axpy(*m_ct, 0, rho, *m_ut2, 0);
#pragma omp barrier

        m_staple->lower(*m_ut2, U, mu, nu);
        axpy(*m_ct, 0, rho, *m_ut2, 0);
#pragma omp barrier
      }
    }

    double rho0 = m_rho[mu + m_Ndim * mu];
    m_proj->project(*m_ut2, rho0, *m_ct, *m_ut1);

    if (mu != 3) Field_SF::set_boundary_wk(*m_ut2, m_wk);
#pragma omp barrier

    copy(Usmear, mu, *m_ut2, 0);
#pragma omp barrier
  }
}

//============================================================END=====
