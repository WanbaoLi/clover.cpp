/*!
        @file    forceSmear_HYP.cpp
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Smear/forceSmear_HYP.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = ForceSmear_HYP::register_factory();
}
#endif

const std::string ForceSmear_HYP::class_name = "ForceSmear_HYP";

//====================================================================
void ForceSmear_HYP::init(const Parameters& params)
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

  m_shift = new ShiftField_lex();

  m_U.resize(m_Ndim);

  m_v1.resize(size1());
  m_v2.resize(size2());

  m_Sigma3.resize(size2());
  m_Sigma2.resize(size1b());

  m_iTheta3.resize(m_Ndim);
  m_iTheta2.resize(size2());
  m_iTheta1.resize(size1b());

  m_vt1 = new Field_G(m_Nvol, 1);
  m_vt2 = new Field_G(m_Nvol, 1);
  m_vt3 = new Field_G(m_Nvol, 1);
  m_ct1 = new Field_G(m_Nvol, 1);
  m_ct2 = new Field_G(m_Nvol, 1);
  m_ct3 = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void ForceSmear_HYP::tidyup()
{
  delete m_shift;

  delete m_vt1;
  delete m_vt2;
  delete m_vt3;
  delete m_ct1;
  delete m_ct2;
  delete m_ct3;
}

//====================================================================
void ForceSmear_HYP::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double alpha1, alpha2, alpha3;

  int err = 0;
  err += params.fetch_double("alpha1", alpha1);
  err += params.fetch_double("alpha2", alpha2);
  err += params.fetch_double("alpha3", alpha3);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(alpha1, alpha2, alpha3);
}


//====================================================================
void ForceSmear_HYP::set_parameters(const double alpha1,
                                    const double alpha2,
                                    const double alpha3)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_alpha1 = alpha1;
    m_alpha2 = alpha2;
    m_alpha3 = alpha3;
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  alpha1 = %10.6F\n", m_alpha1);
  vout.general(m_vl, "  alpha2 = %10.6F\n", m_alpha2);
  vout.general(m_vl, "  alpha3 = %10.6F\n", m_alpha3);

}


//====================================================================
void ForceSmear_HYP::get_parameters(Parameters& params) const
{
  params.set_double("alpha1", m_alpha1);
  params.set_double("alpha2", m_alpha2);
  params.set_double("alpha3", m_alpha3);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void ForceSmear_HYP::force_udiv(Field_G& Sigma,
                                const Field_G& Sigmap,
                                const Field_G& U)
{
#pragma omp barrier

  assert(U.nvol() == m_Nvol);
  assert(U.nex() == m_Ndim);
  assert(Sigmap.nvol() == m_Nvol);
  assert(Sigmap.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    copy(m_U[mu], 0, U, mu);
  }
#pragma omp barrier

  // vout.general(m_vl, " smearing step-1\n");
  smear_step1();
  // vout.general(m_vl, " smearing step-2\n");
  smear_step2();

  Sigma.set(0.0);
#pragma omp barrier

  // vout.general(m_vl, " smeared force step-3\n");
  force_step3(Sigma, Sigmap);

  // vout.general(m_vl, " smeared force step-2\n");
  force_step2(Sigma);

  // vout.general(m_vl, " smeared force step-1\n");
  force_step1(Sigma);

  // vout.general(m_vl, " smeared force finished\n");
}


//====================================================================
void ForceSmear_HYP::force_step3(Field_G& Sigma, const Field_G& Sigmap)
{
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {

    m_ct1->set(0.0);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      staple(*m_ct2, m_v2[idx2(mu, nu)], m_v2[idx2(nu, mu)], mu, nu);
      axpy(*m_ct1, 0, 1.0, *m_ct2, 0);
#pragma omp barrier
    }

    scal(*m_ct1, m_alpha1 / 6.0);
    copy(*m_ct2, 0, Sigmap, mu);
#pragma omp barrier

    m_proj->force_recursive(*m_ct3, m_iTheta3[mu],
                            m_alpha1, *m_ct2, *m_ct1, m_U[mu]);
    axpy(Sigma, mu, 1.0, *m_ct3, 0);
#pragma omp barrier
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      force_each(m_Sigma3[idx2(mu, nu)],
                 m_v2[idx2(mu, nu)], m_v2[idx2(nu, mu)],
                 m_iTheta3[mu], m_iTheta3[nu], mu, nu);

      scal(m_Sigma3[idx2(mu, nu)], m_alpha1/6.0);
#pragma omp barrier
    }
  }
}


//====================================================================
void ForceSmear_HYP::force_step2(Field_G& Sigma)
{
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      m_ct1->set(0.0);
#pragma omp barrier

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho == mu) || (rho == nu)) continue;

        staple(*m_ct2, m_v1[idx1(mu, nu, rho)],
               m_v1[idx1(rho, nu, mu)], mu, rho);
        axpy(*m_ct1, 0, 1.0, *m_ct2, 0);
#pragma omp barrier
      }

      scal(*m_ct1, m_alpha2/4.0);
#pragma omp barrier

      m_proj->force_recursive(*m_ct3, m_iTheta2[idx2(mu, nu)],
                              m_alpha2, m_Sigma3[idx2(mu, nu)],
                              *m_ct1, m_U[mu]);
      axpy(Sigma, mu, 1.0, *m_ct3, 0);
#pragma omp barrier
    }
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho == mu) || (rho == nu)) continue;
        force_each(m_Sigma2[idx1b(mu, nu, rho)],
                   m_v1[idx1(mu, nu, rho)], m_v1[idx1(rho, nu, mu)],
                   m_iTheta2[idx2(mu, nu)], m_iTheta2[idx2(rho, nu)],
                   mu, rho);

        scal(m_Sigma2[idx1b(mu, nu, rho)], m_alpha2/4.0);
#pragma omp barrier
      }
    }
  }
}


//====================================================================
void ForceSmear_HYP::force_step1(Field_G& Sigma)
{
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;
      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho == mu) || (rho == nu)) continue;

        int sig = 6 - mu - nu - rho;

        staple(*m_ct1, m_U[mu], m_U[sig], mu, sig);

        scal(*m_ct1, m_alpha3 / 2.0);
#pragma omp barrier

        m_proj->force_recursive(*m_ct2, m_iTheta1[idx1b(mu, nu, rho)],
                                m_alpha3, m_Sigma2[idx1b(mu, nu, rho)],
                                *m_ct1, m_U[mu]);
        axpy(Sigma, mu, 1.0, *m_ct2, 0);
#pragma omp barrier
      }
    }
  }

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho == mu) || (rho == nu)) continue;

        int sig = 6 - mu - nu - rho;

        force_each(*m_ct2, m_U[mu], m_U[sig],
                   m_iTheta1[idx1b(mu, nu, rho)],
                   m_iTheta1[idx1b(sig, nu, rho)],
                   mu, sig);

        scal(*m_ct2, m_alpha3/2.0);
        axpy(Sigma, mu, 1.0, *m_ct2, 0);
#pragma omp barrier
      }
    }
  }
}


//====================================================================
void ForceSmear_HYP::force_each(Field_G& Sigma_mu,
                                const Field_G& V_mu, const Field_G& V_nu,
                                const Field_G& iTheta_mu,
                                const Field_G& iTheta_nu,
                                const int mu, const int nu)
{
#pragma omp barrier

  Sigma_mu.set(0.0);
#pragma omp barrier

  //- The 1st block                                                            
  m_shift->backward(*m_vt1, V_nu, mu);
  m_shift->backward(*m_vt2, V_mu, nu);

  mult_Field_Gnd(*m_vt3, 0, *m_vt1, 0, *m_vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, *m_vt3, 0, iTheta_nu, 0, 1.0);
#pragma omp barrier

  //- The 2nd block
  mult_Field_Gdn(*m_vt3, 0, iTheta_mu, 0, V_nu, 0);
  mult_Field_Gdn(*m_vt2, 0, *m_vt1, 0, *m_vt3, 0);
  m_shift->forward(*m_vt3, *m_vt2, nu);
  axpy(Sigma_mu, 1.0, *m_vt3);
#pragma omp barrier

  //- The 3rd block
  mult_Field_Gdn(*m_vt3, 0, V_mu, 0, iTheta_nu, 0);
  mult_Field_Gdn(*m_vt2, 0, *m_vt1, 0, *m_vt3, 0);
  m_shift->forward(*m_vt3, *m_vt2, nu);
  axpy(Sigma_mu, 1.0, *m_vt3);
#pragma omp barrier

  //- The 4th bloc
  m_shift->backward(*m_vt1, iTheta_nu, mu);
  m_shift->backward(*m_vt2, V_mu, nu);
  mult_Field_Gnd(*m_vt3, 0, *m_vt1, 0, *m_vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, *m_vt3, 0, V_nu, 0, 1.0);
#pragma omp barrier

 //- The 5th block
  mult_Field_Gdd(*m_vt2, 0, *m_vt1, 0, V_mu, 0);
  mult_Field_Gnn(*m_vt3, 0, *m_vt2, 0, V_nu, 0);
  m_shift->forward(*m_vt2, *m_vt3, nu);
  axpy(Sigma_mu, 1.0, *m_vt2);
#pragma omp barrier

  //- The 6th block
  m_shift->backward(*m_vt1, V_nu, mu);
  m_shift->backward(*m_vt2, iTheta_mu, nu);
  mult_Field_Gnd(*m_vt3, 0, *m_vt1, 0, *m_vt2, 0);
  multadd_Field_Gnd(Sigma_mu, 0, *m_vt3, 0, V_nu, 0, 1.0);
#pragma omp barrier

}


//====================================================================
void ForceSmear_HYP::smear_step1()
{
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      for (int rho = nu + 1; rho < m_Ndim; ++rho) {
        if (rho == mu) continue;

        int sig = 6 - mu - nu - rho;

        staple(*m_ct1, m_U[mu], m_U[sig], mu, sig);
        scal(*m_ct1, m_alpha3/2.0);
#pragma omp barrier
        m_proj->project(m_v1[idx1(mu, nu, rho)],
                        m_alpha3, *m_ct1, m_U[mu]);
      }
    }
  }
}


//====================================================================
void ForceSmear_HYP::smear_step2()
{
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      m_ct2->set(0.0);
#pragma omp barrier

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho != mu) && (rho != nu)) {
          staple(*m_ct1, m_v1[idx1(mu, nu, rho)],
                        m_v1[idx1(rho, nu, mu)], mu, rho);
          axpy(*m_ct2, 0, 1.0, *m_ct1, 0);
#pragma omp barrier
        }
      }
      scal(*m_ct2, m_alpha2 / 4.0);
#pragma omp barrier
      m_proj->project(m_v2[idx2(mu, nu)], m_alpha2, *m_ct2, m_U[mu]);
    }
  }
}


//====================================================================
void ForceSmear_HYP::staple(Field_G& c,
                            const Field_G& u_mu, const Field_G& u_nu,
                            const int mu, const int nu)
{
#pragma omp barrier

  //- upper direction
  m_shift->backward(*m_vt1, u_mu, nu);
  mult_Field_Gnn(*m_vt2, 0, u_nu, 0, *m_vt1, 0);
#pragma omp barrier

  m_shift->backward(*m_vt1, u_nu, mu);
  mult_Field_Gnd(c, 0, *m_vt2, 0, *m_vt1, 0);
#pragma omp barrier

  //- lower direction
  m_shift->backward(*m_vt2, u_nu, mu);
  mult_Field_Gnn(*m_vt1, 0, u_mu, 0, *m_vt2, 0);
  mult_Field_Gdn(*m_vt2, 0, u_nu, 0, *m_vt1, 0);
  m_shift->forward(*m_vt1, *m_vt2, nu);
  axpy(c, 0, 1.0, *m_vt1, 0);
#pragma omp barrier
}


//============================================================END=====
