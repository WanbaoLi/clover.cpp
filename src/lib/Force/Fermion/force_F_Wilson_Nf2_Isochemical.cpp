/*!
        @file    force_F_Wilson_Nf2_Isochemical.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "force_F_Wilson_Nf2_Isochemical.h"
#include "lib/Force/Fermion/tensorProd.h"

const std::string Force_F_Wilson_Nf2_Isochemical::class_name
  = "Force_F_Wilson_Nf2_Isochemical";

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  //vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.crucial(m_vl, "%s: construction\n", class_name.c_str());
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

  m_boundary.resize(m_Ndim);

  set_parameters_impl(params);

  m_fopr_w = new Fopr_Wilson_Chemical(params);

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);
  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);

  m_force1 = new Field_G(Nvol, m_Ndim);

  vout.decrease_indent();
  //vout.general(m_vl, "%s: construction finished.\n",
  vout.crucial(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::init(const std::string repr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  m_repr = repr;

  m_fopr_w = new Fopr_Wilson_Chemical(m_repr);

  m_Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  m_boundary.resize(m_Ndim);

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);
  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);

  m_force1 = new Field_G(Nvol, m_Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_fopr_w;

  delete m_eta;
  delete m_zeta;
  delete m_eta2;
  delete m_eta3;

  delete m_force1;

}

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters(
                                           const Parameters& params)
{
  set_parameters_impl(params);

  vout.increase_indent();
  m_fopr_w->set_parameters(params);
  vout.decrease_indent();
}

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters_impl(
                                           const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           kappa, mu;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("chemical_potential", mu);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(kappa, mu, bc);

}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters(
                                         const double kappa,
                                         const double mu,
                                         const std::vector<int> bc)
{
  set_parameters_impl(kappa, mu, bc);

  Parameters params;
  get_parameters(params);
  m_fopr_w->set_parameters(params);

}

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_parameters_impl(
                                        const double kappa,
                                        const double mu,
                                        const std::vector<int> bc)
{
#pragma omp barrier

  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();

  if (ith == 0) {
    m_kappa  = kappa;
    m_mu     = mu;
    m_exp_mu = exp(mu);

    m_boundary.resize(m_Ndim);
    m_boundary = bc;
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  mu    = %12.8f\n", m_mu);
  for (int dir = 0; dir < m_Ndim; ++dir) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", dir, m_boundary[dir]);
  }

}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_config(Field *U)
{
  m_U = (Field_G *)U;
  m_fopr_w->set_config(U);
}

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::set_mode(const std::string& mode)
{
  m_mode = mode;
}

//====================================================================
void Force_F_Wilson_Nf2_Isochemical::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("chemical_potential", m_mu);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::force_udiv(Field& force_, const Field& eta_)
{
#pragma omp barrier

  copy(*m_eta, eta_);
#pragma omp barrier

  m_fopr_w->set_mode("H");
  m_fopr_w->mult(*m_zeta, *m_eta);

  set_mode("H");
  force_udiv1_impl(*m_force1, *m_zeta, *m_eta);
  copy(force_, *m_force1);
#pragma omp barrier

  set_mode("Hdag");
  force_udiv1_impl(*m_force1, *m_eta, *m_zeta);

  axpy(force_, 1.0, *m_force1);
#pragma omp barrier
}


//====================================================================
void Force_F_Wilson_Nf2_Isochemical::force_udiv1(Field& force_,
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
void Force_F_Wilson_Nf2_Isochemical::force_udiv1_impl(
                                                Field_G& force,
                                                const Field_F& zeta,
                                                const Field_F& eta)
{
#pragma omp barrier

  force.set(0.0);
#pragma omp barrier

  for (int dir = 0; dir < m_Ndim-1; ++dir) {
    m_fopr_w->mult_gm5p(dir, *m_eta2, eta);

    mult_Field_Gd(*m_eta3, 0, *m_U, dir, *m_eta2, 0);
    scal(*m_eta3, -m_kappa);
#pragma omp barrier

    tensorProd_Field_F(force, dir, zeta, *m_eta3);
  }

  const int dir = m_Ndim - 1;

  m_fopr_w->mult_gm5p(dir, *m_eta2, eta);

  mult_Field_Gd(*m_eta3, 0, *m_U, dir, *m_eta2, 0);

  if (m_mode == "H") {
    scal(*m_eta3, -(m_kappa * m_exp_mu));
  } else if (m_mode == "Hdag") {
    scal(*m_eta3, -(m_kappa / m_exp_mu));
  } else {
    vout.crucial(m_vl, "Error at %s: illegal mode.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  tensorProd_Field_F(force, dir, zeta, *m_eta3);

}

//============================================================END=====
