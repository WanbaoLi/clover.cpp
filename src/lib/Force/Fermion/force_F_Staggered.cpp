/*!
        @file    force_F_Staggered.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2026-03-17 18:06:41 #$
        @version $LastChangedRevision: 2699 $
*/

#include "lib/Force/Fermion/force_F_Staggered.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Field/field_thread-inc.h"

using Bridge::vout;

const std::string Force_F_Staggered::class_name = "Force_F_Staggered";

//====================================================================
void Force_F_Staggered::init(const Parameters& params)
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

  m_boundary.resize(m_Ndim);

  set_parameters_impl(params);

  m_fopr_ks = new Fopr_Staggered(params);

  m_shift = new ShiftField_lex();

  m_eta  = new Field_F_1spinor(m_Nvol, 1);
  m_zeta = new Field_F_1spinor(m_Nvol, 1);
  m_eta2 = new Field_F_1spinor(m_Nvol, 1);

  m_force1 = new Field_G(m_Nvol, m_Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Force_F_Staggered::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n", class_name.c_str());
  vout.increase_indent();

  m_Ndim  = CommonParameters::Ndim();
  m_Nvol  = CommonParameters::Nvol();

  m_boundary.resize(m_Ndim);

  m_fopr_ks = new Fopr_Staggered();

  m_shift = new ShiftField_lex();

  m_eta  = new Field_F_1spinor(m_Nvol, 1);
  m_zeta = new Field_F_1spinor(m_Nvol, 1);
  m_eta2 = new Field_F_1spinor(m_Nvol, 1);

  m_force1 = new Field_G(m_Nvol, m_Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Force_F_Staggered::tidyup()
{
  delete m_fopr_ks;
  delete m_shift;

  delete m_eta;
  delete m_zeta;
  delete m_eta2;

  delete m_force1;

}


//====================================================================
void Force_F_Staggered::set_parameters(const Parameters& params)
{
  set_parameters_impl(params);

  m_fopr_ks->set_parameters(params);
}


//====================================================================
void Force_F_Staggered::set_parameters_impl(const Parameters& params)
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
    vout.crucial(m_vl, "%s: fetch error, input parameter not found.\n",
                 class_name.c_str());
    abort();
  }

  set_parameters_impl(mq, bc);
}


//====================================================================
void Force_F_Staggered::set_parameters(double mq, const std::vector<int> bc)
{
  set_parameters_impl(mq, bc);

  m_fopr_ks->set_parameters(m_mq, m_boundary);
}


//====================================================================
void Force_F_Staggered::set_parameters_impl(double mq,
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

  //- print input parameters
  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  mq   = %8.4f\n", mq);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

}


//====================================================================
void Force_F_Staggered::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
void Force_F_Staggered::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_ks->set_config(U);
}

//====================================================================
void Force_F_Staggered::force_udiv(Field& force_, const Field& eta_)
{
#pragma omp barrier

  int Nc   = CommonParameters::Nc();

  assert(eta_.nin() == 2 * Nc);
  assert(eta_.nvol() == m_Nvol);
  assert(eta_.nex() == 1);

  copy(*m_eta, eta_);
#pragma omp barrier

  m_fopr_ks->H(*m_zeta, *m_eta);

  force_udiv1_impl(*m_force1, *m_eta, *m_zeta);
  copy(force_, *m_force1);
#pragma omp barrier

  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);
  axpy(force_, 1.0, *m_force1);
#pragma omp barrier
}


//====================================================================
void Force_F_Staggered::force_udiv1(Field& force_,
                                    const Field& eta_,
                                    const Field& zeta_)
{
#pragma omp barrier

  int Nc   = CommonParameters::Nc();

  assert(eta_.nin() == 2 * Nc);
  assert(eta_.nvol() == m_Nvol);
  assert(eta_.nex() == 1);
  assert(zeta_.nin() == 2 * Nc);
  assert(zeta_.nvol() == m_Nvol);
  assert(zeta_.nex() == 1);

  copy(*m_zeta, zeta_);
  copy(*m_eta,  eta_);
#pragma omp barrier

  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);

  copy(force_, *m_force1);
#pragma omp barrier
}


//====================================================================
void Force_F_Staggered::force_udiv1_impl(Field_G& force,
                                         const Field_F_1spinor& zeta,
                                         const Field_F_1spinor& eta)
{
#pragma omp barrier

  int Nc   = CommonParameters::Nc();

  Mat_SU_N        ut(Nc);
  double          utr, uti;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  for (int mu = 0; mu < m_Ndim; ++mu) {
    m_shift->backward(*m_eta2, eta, m_boundary[mu], mu);
    m_fopr_ks->mult_staggered_phase(*m_eta2, mu);
    m_fopr_ks->mult_gm5(*m_eta2);

    for (int site = is; site < ns; ++site) {
      for (int c1 = 0; c1 < Nc; ++c1) {
        for (int c2 = 0; c2 < Nc; ++c2) {
          utr =   zeta.cmp_r(c2, site) * m_eta2->cmp_r(c1, site)
                + zeta.cmp_i(c2, site) * m_eta2->cmp_i(c1, site);
          uti =   zeta.cmp_r(c2, site) * m_eta2->cmp_i(c1, site)
                - zeta.cmp_i(c2, site) * m_eta2->cmp_r(c1, site);
          // utr *= 0.5/m_mq;  // hopping normalization
          // uti *= 0.5/m_mq;
          utr *= 0.5;          // mass normalization
          uti *= 0.5;
          ut.set(c1, c2, utr, uti);
        }
      }

      force.set_mat(site, mu, ut);
    }
#pragma omp barrier
  }

}

//============================================================END=====
