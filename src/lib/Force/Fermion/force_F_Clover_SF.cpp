/*!
        @file    force_F_Clover_SF.cpp
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Force/Fermion/force_F_Clover_SF.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Field/field_thread-inc.h"

const std::string Force_F_Clover_SF::class_name = "Force_F_Clover_SF";

//====================================================================
void Force_F_Clover_SF::init(const Parameters& params)
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

  m_fopr_c  = new Fopr_Clover_SF(params);

  m_force_w = new Force_F_Wilson_SF(params);

  m_staple = new Staple_SF();
  m_staple->set_parameters(m_phi, m_phipr);

  m_shift = new ShiftField_lex();

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);

  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);
  m_zeta1 = new Field_F(Nvol, 1);
  m_zeta2 = new Field_F(Nvol, 1);

  m_Cud = new Field_G(Nvol, m_Ndim * m_Ndim);

  m_force2 = new Field_G(Nvol, m_Ndim);

  m_force1 = new Field_G(Nvol, 1);

  m_Ut1 = new Field_G(Nvol, 1);
  m_Ut2 = new Field_G(Nvol, 1);
  m_Umu = new Field_G(Nvol, 1);
  m_Unu = new Field_G(Nvol, 1);

  m_vt1 = new Field_F(Nvol, 1);
  m_vt2 = new Field_F(Nvol, 1);
  m_vt3 = new Field_F(Nvol, 1);
  m_vt4 = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Clover_SF::init()
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

  m_fopr_c  = new Fopr_Clover_SF;

  m_force_w = new Force_F_Wilson_SF;

  m_staple = new Staple_SF();
  m_staple->set_parameters(m_phi, m_phipr);

  m_shift = new ShiftField_lex();

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);

  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);
  m_zeta1 = new Field_F(Nvol, 1);
  m_zeta2 = new Field_F(Nvol, 1);

  m_Cud = new Field_G(Nvol, m_Ndim * m_Ndim);

  m_force2 = new Field_G(Nvol, m_Ndim);

  m_force1 = new Field_G(Nvol, 1);

  m_Ut1 = new Field_G(Nvol, 1);
  m_Ut2 = new Field_G(Nvol, 1);
  m_Umu = new Field_G(Nvol, 1);
  m_Unu = new Field_G(Nvol, 1);

  m_vt1 = new Field_F(Nvol, 1);
  m_vt2 = new Field_F(Nvol, 1);
  m_vt3 = new Field_F(Nvol, 1);
  m_vt4 = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_Clover_SF::tidyup()
{
  delete m_force_w;
  delete m_fopr_c;

  delete m_staple;
  delete m_shift;

  delete m_eta;
  delete m_zeta;

  delete m_eta2;
  delete m_eta3;
  delete m_zeta1;
  delete m_zeta2;

  delete m_Cud;
  delete m_force1;
  delete m_force2;

  delete m_Ut1;
  delete m_Ut2;
  delete m_Umu;
  delete m_Unu;

  delete m_vt1;
  delete m_vt2;
  delete m_vt3;
  delete m_vt4;
}

//====================================================================
void Force_F_Clover_SF::set_parameters(const Parameters& params)
{
  set_parameters_impl(params);

  //- propagate parameters
  m_fopr_c->set_parameters(params);
  m_force_w->set_parameters(params);
}

//====================================================================
void Force_F_Clover_SF::set_parameters_impl(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double              kappa, cSW;
  std::vector<int>    bc;
  std::vector<double> phi, phipr;

  int err = 0;
  err += params.fetch_double("hopping_parameter", kappa);
  err += params.fetch_double("clover_coefficient", cSW);
  err += params.fetch_int_vector("boundary_condition", bc);
  err += params.fetch_double_vector("phi", phi);
  err += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(kappa, cSW, bc, phi, phipr);

}


//====================================================================
void Force_F_Clover_SF::set_parameters(const double kappa,
                                       const double cSW,
                                       const std::vector<int> bc,
                                       const std::vector<double> phi,
                                       const std::vector<double> phipr)
{
  set_parameters_impl(kappa, cSW, bc, phi, phipr);

  //- propagate parameters
  m_fopr_c->set_parameters(m_kappa, m_cSW, m_boundary, m_phi, m_phipr);
  m_force_w->set_parameters(m_kappa, m_boundary);

}

//====================================================================
void Force_F_Clover_SF::set_parameters_impl(
                                       const double kappa,
                                       const double cSW,
                                       const std::vector<int> bc,
                                       const std::vector<double> phi,
                                       const std::vector<double> phipr)
{
#pragma omp barrier

  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_kappa = kappa;
    m_cSW   = cSW;
    m_boundary = bc;

    m_phi.resize(3);
    m_phipr.resize(3);
    for (int i = 0; i < 3; ++i) {
      m_phi[i]   = phi[i];
      m_phipr[i] = phipr[i];
    }
#pragma omp barrier

    const int Lx     = CommonParameters::Lx();
    double    Lx_inv = 1.0 / double(Lx);
    double    c0r    = cos(phi[0] * Lx_inv);
    double    c0i    = sin(phi[0] * Lx_inv);
    double    c1r    = cos(phi[1] * Lx_inv);
    double    c1i    = sin(phi[1] * Lx_inv);
    double    c2r    = cos(phi[2] * Lx_inv);
    double    c2i    = sin(phi[2] * Lx_inv);
    m_wk.zero();
    m_wk.set(0, 0, c0r, c0i);
    m_wk.set(1, 1, c1r, c1i);
    m_wk.set(2, 2, c2r, c2i);

    c0r = cos(phipr[0] * Lx_inv);
    c0i = sin(phipr[0] * Lx_inv);
    c1r = cos(phipr[1] * Lx_inv);
    c1i = sin(phipr[1] * Lx_inv);
    c2r = cos(phipr[2] * Lx_inv);
    c2i = sin(phipr[2] * Lx_inv);

    m_wkpr.zero();
    m_wkpr.set(0, 0, c0r, c0i);
    m_wkpr.set(1, 1, c1r, c1i);
    m_wkpr.set(2, 2, c2r, c2i);
  }
#pragma omp barrier

  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", m_cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

  vout.general(m_vl, "  phi1  = %12.6f\n", m_phi[0]);
  vout.general(m_vl, "  phi2  = %12.6f\n", m_phi[1]);
  vout.general(m_vl, "  phi3  = %12.6f\n", m_phi[2]);
  vout.general(m_vl, "  phipr1= %12.6f\n", m_phipr[0]);
  vout.general(m_vl, "  phipr2= %12.6f\n", m_phipr[1]);
  vout.general(m_vl, "  phipr3= %12.6f\n", m_phipr[2]);

}

//====================================================================
void Force_F_Clover_SF::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_double_vector("phi", m_phi);
  params.set_double_vector("phipr", m_phipr);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Clover_SF::set_config(Field *U)
{
  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_c->set_config(m_U);
  m_force_w->set_config(U);
  set_component();

}

//====================================================================
void Force_F_Clover_SF::force_udiv(Field& force_, const Field& eta_)
{
#pragma omp barrier

  copy(*m_eta, eta_);
#pragma omp barrier

  m_fopr_c->set_mode("H");
  m_fopr_c->mult(*m_zeta, *m_eta);

  force_udiv1_impl(*m_force2, *m_eta, *m_zeta);

  Field_SF::set_boundary_spatial_link_zero(*m_force2);
#pragma omp barrier

  copy(force_, *m_force2);
#pragma omp barrier
  
  force_udiv1_impl(*m_force2, *m_zeta, *m_eta);

  Field_SF::set_boundary_spatial_link_zero(*m_force2);
#pragma omp barrier

  axpy(force_, 1.0, *m_force2);
#pragma omp barrier

}


//====================================================================
void Force_F_Clover_SF::force_udiv1(Field& force_, const Field& zeta_, const Field& eta_)
{
#pragma omp barrier

  copy(*m_zeta, zeta_);
  copy(*m_eta,  eta_);
#pragma omp barrier

  force_udiv1_impl(*m_force2, *m_zeta, *m_eta);

  copy(force_, *m_force2);
#pragma omp barrier
}


//====================================================================
void Force_F_Clover_SF::force_udiv1_impl(Field_G& force,
                                         const Field_F& zeta,
                                         const Field_F& eta)
{
#pragma omp barrier

  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  force.set(0.0);
#pragma omp barrier

  copy(*m_zeta1, zeta);
#pragma omp barrier

  set_boundary_zero(*m_zeta1);
#pragma omp barrier

  m_force_w->force_udiv1(force, zeta, eta);

  m_fopr_c->mult_gm5(*m_eta2, eta);
  set_boundary_zero(*m_eta2);
#pragma omp barrier

  double fac = -m_kappa * m_cSW / 8.0;

  for (int mu = 0; mu < Ndim; ++mu) {
    for (int nu = 0; nu < Ndim; ++nu) {
      if (nu == mu) continue;

      m_fopr_c->mult_isigma(*m_eta3, *m_eta2, mu, nu);

      copy(*m_Umu, 0, *m_U, mu);
      copy(*m_Unu, 0, *m_U, nu);
#pragma omp barrier

      if (mu != 3) Field_SF::set_boundary_wk(*m_Umu, m_wk);
      if (nu != 3) Field_SF::set_boundary_wk(*m_Unu, m_wk);
#pragma omp barrier

      // R(1) and R(5)
      mult_Field_Gd(*m_vt1, 0, *m_Cud, index_dir(mu, nu), *m_eta3, 0);
      tensorProd_Field_F(*m_force1, *m_zeta1, *m_vt1);
      axpy(force, mu, fac, *m_force1, 0);
#pragma omp barrier

      // R(2)
      mult_Field_Gd(*m_vt3, 0, *m_Umu, 0, *m_eta3, 0);

      m_shift->backward(*m_vt1, *m_vt3, nu);

      m_shift->backward(*m_vt2, *m_zeta1, nu);

      m_shift->backward(*m_Ut1, *m_Unu, mu);

      if (mu == 3){
	Field_SF::set_boundary_wkpr(*m_Ut1, m_wkpr);
#pragma omp barrier
      }

      mult_Field_Gn(*m_vt3, 0, *m_Ut1, 0, *m_vt1, 0);
      mult_Field_Gn(*m_vt4, 0, *m_Unu, 0, *m_vt2, 0);
      tensorProd_Field_F(*m_force1, *m_vt4, *m_vt3);
      axpy(force, mu, fac, *m_force1, 0);
#pragma omp barrier

      // R(4) and R(8)
      m_shift->backward(*m_vt1, *m_eta3, mu);
      m_shift->backward(*m_zeta2, *m_zeta1, mu);
      mult_Field_Gn(*m_vt4, 0, *m_Cud, index_dir(mu, nu), *m_zeta2, 0);
      tensorProd_Field_F(*m_force1, *m_vt4, *m_vt1);
      axpy(force, mu, fac, *m_force1, 0);
#pragma omp barrier

      // R(3)
      m_shift->backward(*m_vt1, *m_eta3, nu);
      mult_Field_Gn(*m_vt3, 0, *m_Unu, 0, *m_vt1, 0);
      mult_Field_Gn(*m_vt4, 0, *m_Umu, 0, *m_zeta2, 0);
      m_shift->backward(*m_vt1, *m_vt3, mu);
      m_shift->backward(*m_vt2, *m_vt4, nu);
      mult_Field_Gn(*m_vt4, 0, *m_Unu, 0, *m_vt2, 0);
      tensorProd_Field_F(*m_force1, *m_vt4, *m_vt1);
      axpy(force, mu, fac, *m_force1, 0);
#pragma omp barrier

      // R(6)
      m_shift->backward(*m_Ut1, *m_Unu, mu);
      if (mu == 3){
	Field_SF::set_boundary_wkpr(*m_Ut1, m_wkpr);
#pragma omp barrier
      }
      mult_Field_Gdd(*m_Ut2, 0, *m_Ut1, 0, *m_Umu, 0);
      mult_Field_Gn(*m_vt1, 0, *m_Ut2, 0, *m_eta3, 0);
      mult_Field_Gd(*m_vt2, 0, *m_Unu, 0, *m_zeta1, 0);
      m_shift->forward(*m_vt3, *m_vt1, nu);
      m_shift->forward(*m_vt4, *m_vt2, nu);
      tensorProd_Field_F(*m_force1, *m_vt4, *m_vt3);
      axpy(force, mu, -fac, *m_force1, 0);
#pragma omp barrier

      // R(7)
      mult_Field_Gd(*m_vt1, 0, *m_Unu, 0, *m_eta3, 0);
      mult_Field_Gn(*m_vt2, 0, *m_Umu, 0, *m_zeta2, 0);
      m_shift->backward(*m_vt3, *m_vt1, mu);
      m_shift->forward(*m_vt1, *m_vt3, nu);
      mult_Field_Gd(*m_vt4, 0, *m_Unu, 0, *m_vt2, 0);
      m_shift->forward(*m_vt2, *m_vt4, nu);
      tensorProd_Field_F(*m_force1, *m_vt2, *m_vt1);
      axpy(force, mu, -fac, *m_force1, 0);
#pragma omp barrier
    }
  }
}


//====================================================================
void Force_F_Clover_SF::set_component()
{
#pragma omp barrier

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      m_staple->upper(*m_Ut1, *m_U, mu, nu);
      copy(*m_Cud, index_dir(mu, nu), *m_Ut1, 0);
#pragma omp barrier

      m_staple->lower(*m_Ut1, *m_U, mu, nu);
      axpy(*m_Cud, index_dir(mu, nu), -1.0, *m_Ut1, 0);
#pragma omp barrier
    }
  }

}

//====================================================================
void Force_F_Clover_SF::set_boundary_zero(Field_F& f)
{
#pragma omp barrier

  int Nc   = CommonParameters::Nc();
  int Nd   = CommonParameters::Nd();
  int Nvol = CommonParameters::Nvol();
  int Svol = Nvol / CommonParameters::Nt();
  int Nc2  = 2 * Nc;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Svol);

  if (Communicator::ipe(3) == 0) {
    for (int site = is; site < ns; ++site) {
      for (int s = 0; s < Nd; ++s) {
        for (int cc = 0; cc < Nc2; ++cc) {
          f.set(cc + Nc2 * s, site, 0, 0.0);
        }
      }
    }
  }

#pragma omp barrier
}


//============================================================END=====
