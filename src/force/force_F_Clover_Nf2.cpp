/*!
        @file    force_F_Clover_Nf2.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "force/force_F_Clover_Nf2.h"
#include "runtime/threadManager.h"

const std::string Force_F_Clover_Nf2::class_name = "Force_F_Clover_Nf2";

//====================================================================
void Force_F_Clover_Nf2::init(const Parameters& params)
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

  std::string repr;
  if (!params.fetch_string("gamma_matrix_type", repr)) {
    m_repr = repr;
  } else {
    m_repr = "Dirac"; // default gamma-matrix type                            
    vout.general(m_vl, "gamma_matrix_type is not given: defalt = %s\n",
                   m_repr.c_str());
  }

  if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
    vout.crucial("Error at %s: unsupported gamma-matrix type: %s\n",
                 class_name.c_str(), m_repr.c_str());
    exit(EXIT_FAILURE);
  }

  m_Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  m_fopr_c    = new Fopr_Clover(params);
  m_force_w   = new Force_F_Wilson_Nf2(params);
  m_force_csw = new Force_F_CloverTerm(params);

  m_boundary.resize(m_Ndim);

  set_parameters_impl(params);

  m_staple = new Staple_lex();

  m_force1 = new Field_G(Nvol, m_Ndim);
  m_force2 = new Field_G(Nvol, m_Ndim);

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Force_F_Clover_Nf2::init(std::string repr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n",
	       class_name.c_str());
  vout.increase_indent();

  m_repr = repr;

  const int Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_fopr_c    = new Fopr_Clover(repr);
  m_force_w   = new Force_F_Wilson_Nf2(repr);
  m_force_csw = new Force_F_CloverTerm(repr);

  m_boundary.resize(m_Ndim);

  m_staple = new Staple_lex();

  m_force1 = new Field_G(Nvol, m_Ndim);
  m_force2 = new Field_G(Nvol, m_Ndim);

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Force_F_Clover_Nf2::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_fopr_c;
  delete m_force_w;
  delete m_force_csw;

  delete m_staple;

  delete m_force1;
  delete m_force2;

  delete m_eta;
  delete m_zeta;
}


//====================================================================
void Force_F_Clover_Nf2::set_parameters(const Parameters& params)
{
  set_parameters_impl(params);

  m_fopr_c->set_parameters(params);
  m_force_w->set_parameters(params);
  m_force_csw->set_parameters(params);
}


//====================================================================
void Force_F_Clover_Nf2::set_parameters_impl(const Parameters& params)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    std::string vlevel;
    if (!params.fetch_string("verbose_level", vlevel)) {
      m_vl = vout.set_verbose_level(vlevel);
    }
  }
#pragma omp barrier

  //- fetch and check input parameters
  double           kappa, cSW;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(kappa, cSW, bc);
}


//====================================================================
void Force_F_Clover_Nf2::set_parameters(const double kappa,
                                        const double cSW,
                                        const std::vector<int> bc)
{
  set_parameters_impl(kappa, cSW, bc);

  //- propagate parameters
  Parameters params;
  get_parameters(params);

  m_fopr_c->set_parameters(params);
  m_force_w->set_parameters(params);
  m_force_csw->set_parameters(params);

}


//====================================================================
void Force_F_Clover_Nf2::set_parameters_impl(const double kappa,
                                             const double cSW,
                                             const std::vector<int> bc)
{
#pragma omp barrier

  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_kappa = kappa;
    m_cSW   = cSW;
    m_boundary = bc;
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: parameters:\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

}


//====================================================================
void Force_F_Clover_Nf2::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Clover_Nf2::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_c->set_config(U);
  m_force_w->set_config(U);
  m_force_csw->set_config(U);
}

//====================================================================
void Force_F_Clover_Nf2::force_udiv(Field& force_, const Field& eta_)
{
#pragma omp barrier

  copy(*m_eta, eta_);
#pragma omp barrier

  m_fopr_c->H(*m_zeta, *m_eta);

  force_udiv1_impl(*m_force1, *m_eta, *m_zeta);
  copy(force_, *m_force1);
#pragma omp barrier

  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);
  axpy(force_, 1.0, *m_force1);
#pragma omp barrier

}


//====================================================================
void Force_F_Clover_Nf2::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
{
#pragma omp barrier

  copy(*m_eta, eta_);
  copy(*m_zeta, zeta_);
#pragma omp barrier

  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);

  copy(force_, *m_force1);
#pragma omp barrier
}


//====================================================================
void Force_F_Clover_Nf2::force_udiv1_impl(Field_G& force,
                                          const Field_F& zeta,
                                          const Field_F& eta)
{
#pragma omp barrier

  force.set(0.0);
#pragma omp barrier

  m_force_w->set_mode("H");
  m_force_w->force_udiv1(force, zeta, eta);

  m_force_csw->force_udiv1(*m_force2, zeta, eta);

  axpy(force, 1.0, *m_force2);
#pragma omp barrier
}


//============================================================END=====
