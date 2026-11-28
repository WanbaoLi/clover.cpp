/*!
        @file    force_F_Wilson_SF.cpp
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Force/Fermion/force_F_Wilson_SF.h"
#include "lib/ResourceManager/threadManager.h"

const std::string Force_F_Wilson_SF::class_name = "Force_F_Wilson_SF";

//====================================================================
void Force_F_Wilson_SF::init(const Parameters& params)
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

  const int Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  set_parameters_impl(params);

  m_fopr_w = new Fopr_Wilson_SF(params);

  m_force1 = new Field_G(Nvol, m_Ndim);

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);

  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Wilson_SF::init()
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  m_repr = "Dirac"; // default gamma-matrix type

  const int Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();

  m_boundary.resize(m_Ndim);

  m_fopr_w = new Fopr_Wilson_SF();

  m_force1 = new Field_G(Nvol, m_Ndim);

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);

  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Wilson_SF::tidyup()
{
    delete m_fopr_w;

    delete m_force1;

    delete m_eta;
    delete m_zeta;
    delete m_eta2;
    delete m_eta3;
}

//====================================================================
void Force_F_Wilson_SF::set_parameters(const Parameters& params)
{
  set_parameters_impl(params);

  m_fopr_w->set_parameters(params);
}

//====================================================================
void Force_F_Wilson_SF::set_parameters_impl(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(kappa, bc);
}


//====================================================================
void Force_F_Wilson_SF::set_parameters(const double kappa,
                                       const std::vector<int> bc)
{
  set_parameters_impl(kappa, bc);

  //- propagate parameters
  m_fopr_w->set_parameters(m_kappa, m_boundary);
}


//====================================================================
void Force_F_Wilson_SF::set_parameters_impl(const double kappa,
                                            const std::vector<int> bc)
{
#pragma omp barrier

  assert(bc.size() == m_Ndim);
  
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_kappa = kappa;
    m_boundary = bc;
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}


//====================================================================
void Force_F_Wilson_SF::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Wilson_SF::set_config(Field *U)
{
  int ith = ThreadManager::get_thread_id();
  if (ith == 0)  m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_w->set_config(U);
}

//====================================================================
void Force_F_Wilson_SF::force_udiv(Field& force_, const Field& eta_)
{
#pragma omp barrier

  copy(*m_eta, eta_);
#pragma omp barrier

  m_fopr_w->set_mode("H");
  m_fopr_w->mult(*m_zeta, *m_eta);

  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);

  Field_SF::set_boundary_spatial_link_zero(*m_force1);
#pragma omp barrier

  copy(force_, *m_force1);
#pragma omp barrier

  force_udiv1_impl(*m_force1, *m_eta, *m_zeta);

  Field_SF::set_boundary_spatial_link_zero(*m_force1);
#pragma omp barrier

  axpy(force_, 1.0, *m_force1);
#pragma omp barrier

}


//====================================================================
void Force_F_Wilson_SF::force_udiv1(Field& force_,
                                    const Field& zeta_,
                                    const Field& eta_)
{
#pragma omp barrier

  copy(*m_zeta, zeta_);
  copy(*m_eta,  eta_);

  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);

  Field_SF::set_boundary_spatial_link_zero(*m_force1);
  // added [11 Jan 2026 H.Matsufuru]
#pragma omp barrier

  copy(force_, *m_force1);
#pragma omp barrier
}


//====================================================================
void Force_F_Wilson_SF::force_udiv1_impl(Field_G& force,
                                         const Field_F& zeta,
                                         const Field_F& eta)
{
#pragma omp barrier

  force.set(0.0);  // this may not be necessary
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {

    m_fopr_w->mult_gm5p(mu, *m_eta2, eta);

    mult_Field_Gd(*m_eta3, 0, *m_U, mu, *m_eta2, 0);

    tensorProd_Field_F(force, mu, *m_zeta, *m_eta3);
  }

  scal(force, -m_kappa);
#pragma omp barrier
}

//============================================================END=====
