/*!
        @file    smear_HYP.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Smear/smear_HYP.h"
#include "lib/ResourceManager/threadManager.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Smear_HYP::register_factory();
}
#endif

const std::string Smear_HYP::class_name = "Smear_HYP";

//====================================================================
void Smear_HYP::init(const Parameters& params)
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
  m_v1.resize(size_v1());
  m_v2.resize(size_v2());
  
  m_c1 = new Field_G(m_Nvol, 1);
  m_c2 = new Field_G(m_Nvol, 1);
  m_c3 = new Field_G(m_Nvol, 1);
  m_ut = new Field_G(m_Nvol, 1);
  m_vt1 = new Field_G(m_Nvol, 1);
  m_vt2 = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Smear_HYP::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n", class_name.c_str());
  vout.increase_indent();

  m_Ndim = CommonParameters::Ndim();
  m_Nvol = CommonParameters::Nvol();

  m_shift = new ShiftField_lex();

  m_U.resize(m_Ndim);
  m_v1.resize(size_v1());
  m_v2.resize(size_v2());

  m_c1 = new Field_G(m_Nvol, 1);
  m_c2 = new Field_G(m_Nvol, 1);
  m_c3 = new Field_G(m_Nvol, 1);
  m_ut = new Field_G(m_Nvol, 1);
  m_vt1 = new Field_G(m_Nvol, 1);
  m_vt2 = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Smear_HYP::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_shift;

  delete m_c1;
  delete m_c2;
  delete m_c3;
  delete m_ut;
  delete m_vt1;
  delete m_vt2;

}

//====================================================================
void Smear_HYP::set_parameters(const Parameters& params)
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
void Smear_HYP::set_parameters(const double alpha1,
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

  //- print input parameters
  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  alpha1 = %10.6F\n", alpha1);
  vout.general(m_vl, "  alpha2 = %10.6F\n", alpha2);
  vout.general(m_vl, "  alpha3 = %10.6F\n", alpha3);
}


//====================================================================
void Smear_HYP::get_parameters(Parameters& params) const
{
  params.set_double("alpha1", m_alpha1);
  params.set_double("alpha2", m_alpha2);
  params.set_double("alpha3", m_alpha3);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Smear_HYP::smear(Field_G& Usmear, const Field_G& U)
{
  assert(U.nvol() == m_Nvol);
  assert(U.nex() == m_Ndim);

  assert(Usmear.nvol() == m_Nvol);
  assert(Usmear.nex() == m_Ndim);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    copy(m_U[mu], 0, U, mu);
  }
#pragma omp barrier

  step1();
  //  vout.general(m_vl, "level-1 step finished.\n");
  step2();
  //  vout.general(m_vl, "level-2 step finished.\n");
  step3(Usmear);
  //  vout.general(m_vl, "level-3 step finished.\n");
}


//====================================================================
void Smear_HYP::step1()
{
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      for (int rho = nu + 1; rho < m_Ndim; ++rho) {
        if (rho == mu) continue;

        int sig = 6 - mu - nu - rho;

        staple(*m_c1, m_U[mu], m_U[sig], mu, sig);

        scal(*m_c1, m_alpha3 * 0.5); // c1 *= m_alpha3 / 2.0;
#pragma omp barrier

        m_proj->project(m_v1[index_v1(mu, nu, rho)], m_alpha3,
                        *m_c1, m_U[mu]);
      }
    }
  }

#pragma omp barrier
}


//====================================================================
void Smear_HYP::step2()
{
  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      m_c2->set(0.0);
#pragma omp barrier

      for (int rho = 0; rho < m_Ndim; ++rho) {
        if ((rho != mu) && (rho != nu)) {

          staple(*m_ut, m_v1[index_v1(mu, nu, rho)],
                 m_v1[index_v1(rho, nu, mu)], mu, rho);

          axpy(*m_c2, 0, 1.0, *m_ut, 0);
#pragma omp barrier
        }
      }

      scal(*m_c2, m_alpha2 * 0.25); //c2 *= m_alpha2 / 4.0;
#pragma omp barrier

      m_proj->project(m_v2[index_v2(mu, nu)], m_alpha2, *m_c2, m_U[mu]);
    }
  }

#pragma omp barrier
}


//====================================================================
void Smear_HYP::step3(Field_G& v)
{
  for (int mu = 0; mu < m_Ndim; ++mu) {

    m_c3->set(0.0);
#pragma omp barrier

    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu != mu) {
        staple(*m_ut, m_v2[index_v2(mu, nu)],
                      m_v2[index_v2(nu, mu)], mu, nu);

        axpy(*m_c3, 0, 1.0, *m_ut, 0);
#pragma omp barrier
      }
    }

    scal(*m_c3, m_alpha1 / 6.0); //c3 *= m_alpha1 / 6.0;
#pragma omp barrier

    m_proj->project(*m_ut, m_alpha1, *m_c3, m_U[mu]);

    copy(v, mu, *m_ut, 0);
#pragma omp barrier
  }

}


//====================================================================
void Smear_HYP::staple(Field_G& c,
                       const Field_G& u_mu, const Field_G& u_nu,
                       const int mu, const int nu)
{
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
#pragma omp barrier

  mult_Field_Gdn(*m_vt2, 0, u_nu, 0, *m_vt1, 0);
#pragma omp barrier

  m_shift->forward(*m_vt1, *m_vt2, nu);

  axpy(c, 0, 1.0, *m_vt1, 0);
#pragma omp barrier

}


//============================================================END=====
