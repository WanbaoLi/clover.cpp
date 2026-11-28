/*!
        @file    force_F_Staggered_eo.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2026-03-17 18:06:41 #$
        @version $LastChangedRevision: 2699 $
*/

#include "lib/Force/Fermion/force_F_Staggered_eo.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Field/field_thread-inc.h"

const std::string Force_F_Staggered_eo::class_name
                                           = "Force_F_Staggered_eo";

//====================================================================
void Force_F_Staggered_eo::init(const Parameters& params)
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

  m_Ndim  = CommonParameters::Ndim();
  m_Nvol  = CommonParameters::Nvol();
  m_Nvol2 = m_Nvol/2;

  m_boundary.resize(m_Ndim);

  set_parameters_impl(params);
  
  m_fopr_ks = new Fopr_Staggered_eo(params);

  m_shift_eo = new ShiftField_eo();

  m_index_eo = new Index_eo();

  m_Ueo = new Field_G(m_Nvol, m_Ndim);

  m_force_e  = new Field_G(m_Nvol2, m_Ndim);
  m_force_o  = new Field_G(m_Nvol2, m_Ndim);
  m_force_eo = new Field_G(m_Nvol, m_Ndim);

  m_eta  = new Field_F_1spinor(m_Nvol2, 1);
  m_zeta = new Field_F_1spinor(m_Nvol2, 1);
  m_eta2 = new Field_F_1spinor(m_Nvol2, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Force_F_Staggered_eo::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n", class_name.c_str());
  vout.increase_indent();

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  m_fopr_ks = new Fopr_Staggered_eo();

  m_shift_eo = new ShiftField_eo();

  m_index_eo = new Index_eo();

  m_Ueo = new Field_G(m_Nvol, m_Ndim);

  m_force_e  = new Field_G(m_Nvol2, m_Ndim);
  m_force_o  = new Field_G(m_Nvol2, m_Ndim);
  m_force_eo = new Field_G(m_Nvol, m_Ndim);

  m_eta  = new Field_F_1spinor(m_Nvol2, 1);
  m_zeta = new Field_F_1spinor(m_Nvol2, 1);
  m_eta2 = new Field_F_1spinor(m_Nvol2, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Force_F_Staggered_eo::tidyup()
{

  delete m_fopr_ks;
  delete m_shift_eo;
  delete m_index_eo;

  delete m_Ueo;

  delete m_force_e;
  delete m_force_o;
  delete m_force_eo;

  delete m_eta;
  delete m_zeta;
  delete m_eta2;

}


//====================================================================
void Force_F_Staggered_eo::set_parameters(const Parameters& params)
{
  set_parameters_impl(params);

  m_fopr_ks->set_parameters(params);

}


//====================================================================
void Force_F_Staggered_eo::set_parameters_impl(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           mq;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(mq, bc);
}


//====================================================================
void Force_F_Staggered_eo::set_parameters(const double mq,
                                          const std::vector<int> bc)
{
  set_parameters_impl(mq, bc);

  m_fopr_ks->set_parameters(m_mq, m_boundary);

}


//====================================================================
void Force_F_Staggered_eo::set_parameters_impl(const double mq,
                                               const std::vector<int> bc)
{
#pragma omp barrier

  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_mq = mq;
    m_boundary = bc;
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %12.8f\n", mq);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

}


//====================================================================
void Force_F_Staggered_eo::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Staggered_eo::set_config(Field *U)
{
  m_index_eo->convertField(*m_Ueo, *U);
  // this method is multi-threaded with barriers.

  m_fopr_ks->set_config(U);
}


//====================================================================
void Force_F_Staggered_eo::force_core(Field& force_, const Field& eta_)
{
#pragma omp barrier

  const int Nc    = CommonParameters::Nc();

  assert(eta_.nin() == 2 * Nc);
  assert(eta_.nvol() == m_Nvol2);
  assert(eta_.nex() == 1);

  copy(*m_eta, eta_);
#pragma omp barrier

  m_fopr_ks->mult(*m_zeta, *m_eta, "Doe");

  force_udiv1(*m_force_e, *m_eta, *m_zeta, 0);

  force_udiv1(*m_force_o, *m_zeta, *m_eta, 1);
  scal(*m_force_o, -1.0);
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int site = 0; site < m_Nvol2; ++site) {
      int gx = m_index_eo->site(site, 0);

      Mat_SU_N ut(Nc);
      ut = m_Ueo->mat(gx, mu) * m_force_e->mat(site, mu);
      ut.at();
      ut *= -2.0;
      m_force_eo->set_mat(gx, mu, ut);

      gx = m_index_eo->site(site, 1);
      ut = m_Ueo->mat(gx, mu) * m_force_o->mat(site, mu);
      ut.at();
      ut *= -2.0;
      m_force_eo->set_mat(gx, mu, ut);
    }
  }
#pragma omp barrier

  m_index_eo->reverseField(force_, *m_force_eo);
}


//====================================================================
void Force_F_Staggered_eo::force_udiv(Field& force_, const Field& eta_)
{
#pragma omp barrier

  const int Nc    = CommonParameters::Nc();

  assert(eta_.nin() == 2 * Nc);
  assert(eta_.nvol() == m_Nvol2);
  assert(eta_.nex() == 1);

  copy(*m_eta, eta_);
#pragma omp barrier

  m_fopr_ks->mult(*m_zeta, *m_eta, "Doe");

  force_udiv1(*m_force_e, *m_eta, *m_zeta, 0);

  force_udiv1(*m_force_o, *m_zeta, *m_eta, 1);
  scal(*m_force_o, -1.0);
#pragma omp barrier

  m_index_eo->reverseField(force_, *m_force_e, 0);
  m_index_eo->reverseField(force_, *m_force_o, 1);
}


//====================================================================
void Force_F_Staggered_eo::force_udiv1(Field_G& force_h,
                                       const Field_F_1spinor& zeta,
                                       const Field_F_1spinor& eta,
                                       int ieo)
{
#pragma omp barrier

  const int Nc    = CommonParameters::Nc();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol2);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_shift_eo->backward_h(*m_eta2, eta, m_boundary[mu], mu, ieo);
    m_fopr_ks->mult_staggered_phase(*m_eta2, mu, ieo);

    for (int site = is; site < ns; ++site) {
      Mat_SU_N ut(Nc);

      for (int c1 = 0; c1 < Nc; ++c1) {
        for (int c2 = 0; c2 < Nc; ++c2) {
          double ut_r =  zeta.cmp_r(c2, site) * m_eta2->cmp_r(c1, site)
                       + zeta.cmp_i(c2, site) * m_eta2->cmp_i(c1, site);
          double ut_i =  zeta.cmp_r(c2, site) * m_eta2->cmp_i(c1, site)
                       - zeta.cmp_i(c2, site) * m_eta2->cmp_r(c1, site);
          // ut_r *= -0.5 / m_mq;  // hopping normalization
          // ut_i *= -0.5 / m_mq;
          ut_r *= -0.5;   // mass normalization
          ut_i *= -0.5;

          ut.set(c1, c2, ut_r, ut_i);
        }
      }

      force_h.set_mat(site, mu, ut);
    }
#pragma omp barrier
  }

}


//====================================================================
void Force_F_Staggered_eo::force_core1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: unimplemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_Staggered_eo::force_udiv1(Field&, const Field&, const Field&)
{
  vout.crucial(m_vl, "Error at %s: unimplemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//============================================================END=====
