/*!
        @file    force_F_Wilson_TwistedMass_Nf2.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Force/Fermion/force_F_Wilson_TwistedMass_Nf2.h"
#include "lib/ResourceManager/threadManager.h"

const std::string Force_F_Wilson_TwistedMass_Nf2::class_name
                                = "Force_F_Wilson_TwistedMass_Nf2";

//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::init(const Parameters& params)
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

  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  m_boundary.resize(Ndim);

  set_parameters_impl(params);

  m_fopr_tw = new Fopr_Wilson_TwistedMass(params);

  m_force1 = new Field_G(Nvol, Ndim);

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);
  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::init(const std::string repr)
{

  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  m_repr = repr;

  const int Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  m_boundary.resize(Ndim);

  m_fopr_tw = new Fopr_Wilson_TwistedMass(m_repr);

  m_force1 = new Field_G(Nvol, Ndim);

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);
  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::tidyup()
{
  delete m_fopr_tw;

  delete m_force1;
  delete m_eta;
  delete m_zeta;
  delete m_eta2;
  delete m_eta3;
}

//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::set_parameters(const Parameters& params)
{
  set_parameters_impl(params);

  m_fopr_tw->set_parameters(params);
}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::set_parameters_impl(
                                            const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa, tw_mass;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("twisted_mass", tw_mass);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(kappa, tw_mass, bc);
}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::set_parameters(
                                           const double kappa,
                                           const double tw_mass,
                                           const std::vector<int> bc)
{
  set_parameters_impl(kappa, tw_mass, bc);

  m_fopr_tw->set_parameters(m_kappa, m_tw_mass, m_boundary);
}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::set_parameters_impl(
                                           const double kappa,
                                           const double tw_mass,
                                           const std::vector<int> bc)
{
#pragma omp barrier

  const int Ndim = CommonParameters::Ndim();

  assert(bc.size() == Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_kappa    = kappa;
    m_tw_mass  = tw_mass;
    m_boundary = bc;
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma-matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa   = %12.8f\n", m_kappa);
  vout.general(m_vl, "  tw_mass = %12.8f\n", m_tw_mass);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("twisted_mass", m_tw_mass);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_tw->set_config(U);
}

//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::force_udiv(Field& force_,
                                                const Field& eta_)
{
#pragma omp barrier

  copy(*m_eta, eta_);
#pragma omp barrier

  m_fopr_tw->set_mode("H");
  m_fopr_tw->mult(*m_zeta, *m_eta);

  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);
  copy(force_, *m_force1);
#pragma omp barrier

  force_udiv1_impl(*m_force1, *m_eta, *m_zeta);
  axpy(force_, 1.0, *m_force1);
#pragma omp barrier

}


//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::force_udiv1(Field& force_,
                                                 const Field& zeta_,
                                                 const Field& eta_)
{
#pragma omp barrier

  copy(*m_zeta, zeta_);
  copy(*m_eta,  eta_);
#pragma omp barrier

  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);

  copy(force_, *m_force1);
#pragma omp barrier

}

//====================================================================
void Force_F_Wilson_TwistedMass_Nf2::force_udiv1_impl(
                                               Field_G& force,
                                               const Field_F& zeta,
                                               const Field_F& eta)
{
#pragma omp barrier

  const int Ndim = CommonParameters::Ndim();

  for (int mu = 0; mu < Ndim; ++mu) {
    m_fopr_tw->mult_gm5p(mu, *m_eta2, eta);

    mult_Field_Gd(*m_eta3, 0, *m_U, mu, *m_eta2, 0);
#pragma omp barrier

    tensorProd_Field_F(force, mu, zeta, *m_eta3);
  }

  scal(force, -m_kappa);
#pragma omp barrier
}


//============================================================END=====
