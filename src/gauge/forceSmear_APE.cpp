/*!
        @file    forceSmear_APE.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "gauge/forceSmear_APE.h"
#include "runtime/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ForceSmear_APE::register_factory();
}
#endif

const std::string ForceSmear_APE::class_name = "ForceSmear_APE";

//====================================================================
void ForceSmear_APE::init(const Parameters& params)
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
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      m_rho[index(mu, nu)] = 0.0;
    }
  }

  set_parameters(params);

  m_shift = new ShiftField_lex();

  m_U.resize(m_Ndim);
  m_iTheta.resize(m_Ndim);

  m_vt1 = new Field_G(m_Nvol, 1);
  m_vt2 = new Field_G(m_Nvol, 1);
  m_vt3 = new Field_G(m_Nvol, 1);
  m_v1  = new Field_G(m_Nvol, 1);
  m_v2  = new Field_G(m_Nvol, 1);

  m_ct1 = new Field_G(m_Nvol, 1);
  m_ct2 = new Field_G(m_Nvol, 1);
  m_xi  = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void ForceSmear_APE::tidyup()
{
  delete m_shift;

  delete m_vt1;
  delete m_vt2;
  delete m_vt3;
  delete m_v1;
  delete m_v2;

  delete m_ct1;
  delete m_ct2;
  delete m_xi;
}

//====================================================================
void ForceSmear_APE::set_parameters(const Parameters& params)
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
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(rho1);
}


//====================================================================
void ForceSmear_APE::set_parameters(const double rho1)
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
  vout.general(m_vl, "  rho = %8.4f\n", rho1);

}


//====================================================================
void ForceSmear_APE::set_parameters(const std::vector<double>& rho)
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
#pragma	omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  rho[%d] = %8.4f\n", mu, rho[mu]);
  }

}


//====================================================================
void ForceSmear_APE::get_parameters(Parameters& params) const
{
  params.set_double_vector("rho", m_rho);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void ForceSmear_APE::force_udiv(Field_G& Sigma,
                                const Field_G& Sigmap,
                                const Field_G& U)
{
  const int Nc = CommonParameters::Nc();

  assert(Sigmap.nin() == (2 * Nc * Nc));
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    copy(m_U[mu], 0, U, mu);
  }
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_ct1->set(0.0);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      double rho = m_rho[index(mu, nu)];

      staple(*m_ct2, m_U[mu], m_U[nu], mu, nu);
      axpy(*m_ct1, 0, rho, *m_ct2, 0);
#pragma omp barrier
    }

    copy(*m_ct2, 0, Sigmap, mu);
#pragma omp barrier

    double alpha = m_rho[index(mu, mu)];

    m_proj->force_recursive(*m_xi, m_iTheta[mu],
                            alpha, *m_ct2, *m_ct1, m_U[mu]);
    copy(Sigma, mu, *m_xi, 0);
#pragma omp barrier
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      double rho = m_rho[index(mu, nu)];

      force_each(*m_ct2, m_U[mu], m_U[nu],
                 m_iTheta[mu], m_iTheta[nu], mu, nu);
      axpy(Sigma, mu, rho, *m_ct2, 0);
#pragma omp barrier
    }
  }
}


//====================================================================
void ForceSmear_APE::force_each(Field_G& Sigma_mu,
                                const Field_G& V_mu,
                                const Field_G& V_nu,
                                const Field_G& iTheta_mu,
                                const Field_G& iTheta_nu,
                                const int mu, const int nu)
{
#pragma omp barrier

  Sigma_mu.set(0.0);
#pragma omp barrier

  m_shift->backward(*m_vt1, V_nu, mu);
  m_shift->backward(*m_vt2, V_mu, nu);
  mult_Field_Gnd(*m_vt3, 0, *m_vt1, 0, *m_vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, *m_vt3, 0, iTheta_nu, 0, 1.0);
#pragma omp barrier

  mult_Field_Gdn(*m_vt3, 0, iTheta_mu, 0, V_nu, 0);
  mult_Field_Gdn(*m_vt2, 0, *m_vt1, 0, *m_vt3, 0);
  m_shift->forward(*m_vt3, *m_vt2, nu);
  axpy(Sigma_mu, 1.0, *m_vt3);
#pragma omp barrier

  mult_Field_Gdn(*m_vt3, 0, V_mu, 0, iTheta_nu, 0);
  mult_Field_Gdn(*m_vt2, 0, *m_vt1, 0, *m_vt3, 0);
  m_shift->forward(*m_vt3, *m_vt2, nu);
  axpy(Sigma_mu, 1.0, *m_vt3);
#pragma omp barrier

  m_shift->backward(*m_vt1, iTheta_nu, mu);
  m_shift->backward(*m_vt2, V_mu, nu);
  mult_Field_Gnd(*m_vt3, 0, *m_vt1, 0, *m_vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, *m_vt3, 0, V_nu, 0, 1.0);
#pragma omp barrier

  mult_Field_Gdd(*m_vt2, 0, *m_vt1, 0, V_mu, 0);
  mult_Field_Gnn(*m_vt3, 0, *m_vt2, 0, V_nu, 0);
  m_shift->forward(*m_vt2, *m_vt3, nu);
  axpy(Sigma_mu, 1.0, *m_vt2);
#pragma omp barrier

  m_shift->backward(*m_vt1, V_nu, mu);
  m_shift->backward(*m_vt2, iTheta_mu, nu);
  mult_Field_Gnd(*m_vt3, 0, *m_vt1, 0, *m_vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, *m_vt3, 0, V_nu, 0, 1.0);
#pragma omp barrier

}


//====================================================================
void ForceSmear_APE::staple(Field_G& c,
                            const Field_G& u_mu, const Field_G& u_nu,
                            const int mu, const int nu)
{
  //- upper direction
  m_shift->backward(*m_v1, u_mu, nu);
  mult_Field_Gnn(*m_v2, 0, u_nu, 0, *m_v1, 0);
  m_shift->backward(*m_v1, u_nu, mu);
  mult_Field_Gnd(c, 0, *m_v2, 0, *m_v1, 0);
#pragma omp barrier

  //- lower direction
  m_shift->backward(*m_v2, u_nu, mu);
  mult_Field_Gnn(*m_v1, 0, u_mu, 0, *m_v2, 0);
  mult_Field_Gdn(*m_v2, 0, u_nu, 0, *m_v1, 0);
  m_shift->forward(*m_v1, *m_v2, nu);
  axpy(c, 0, 1.0, *m_v1, 0);
#pragma omp barrier
}


//============================================================END=====
