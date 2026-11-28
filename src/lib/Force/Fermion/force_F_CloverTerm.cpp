/*!
        @file    force_F_CloverTerm.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: kanamori $
        @date    $LastChangedDate:: 2026-03-17 18:06:41 #$
        @version $LastChangedRevision: 2699 $
*/

#include "lib/Force/Fermion/force_F_CloverTerm.h"
#include "lib/ResourceManager/threadManager.h"

const std::string Force_F_CloverTerm::class_name = "Force_F_CloverTerm";

//====================================================================
void Force_F_CloverTerm::init(const Parameters& params)
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

  m_boundary.resize(m_Ndim);

  //m_fopr_csw = new Fopr_CloverTerm(m_repr);
  m_fopr_csw = new Fopr_CloverTerm(params);

  set_parameters_impl(params);

  m_shift  = new ShiftField_lex();
  m_staple = new Staple_lex();

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);

  m_force1 = new Field_G(Nvol, 1);
  m_force2 = new Field_G(Nvol, m_Ndim);

  m_Cud = new Field_G(Nvol, m_Ndim * m_Ndim);
  m_Ut1 = new Field_G(Nvol, 1);
  m_Ut2 = new Field_G(Nvol, 1);
  m_Unu = new Field_G(Nvol, 1);

  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);
  m_zeta_mu = new Field_F(Nvol, 1);

  m_vt1 = new Field_F(Nvol, 1);
  m_vt2 = new Field_F(Nvol, 1);
  m_vt3 = new Field_F(Nvol, 1);
  m_vt4 = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Force_F_CloverTerm::init(const std::string repr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n",
               class_name.c_str());
  vout.increase_indent();

  m_repr = repr;

  m_Ndim = CommonParameters::Ndim();
  const int Nvol = CommonParameters::Nvol();

  m_boundary.resize(m_Ndim);

  m_fopr_csw = new Fopr_CloverTerm(m_repr);

  m_shift    = new ShiftField_lex();
  m_staple   = new Staple_lex();

  m_eta  = new Field_F(Nvol, 1);
  m_zeta = new Field_F(Nvol, 1);

  m_force1 = new Field_G(Nvol, 1);
  m_force2 = new Field_G(Nvol, m_Ndim);

  m_Cud = new Field_G(Nvol, m_Ndim * m_Ndim);
  m_Ut1 = new Field_G(Nvol, 1);
  m_Ut2 = new Field_G(Nvol, 1);
  m_Unu = new Field_G(Nvol, 1);

  m_eta2 = new Field_F(Nvol, 1);
  m_eta3 = new Field_F(Nvol, 1);
  m_zeta_mu = new Field_F(Nvol, 1);

  m_vt1 = new Field_F(Nvol, 1);
  m_vt2 = new Field_F(Nvol, 1);
  m_vt3 = new Field_F(Nvol, 1);
  m_vt4 = new Field_F(Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());
}


//====================================================================
void Force_F_CloverTerm::tidyup()
{
  ThreadManager::assert_single_thread(class_name);

  delete m_fopr_csw;
  delete m_shift;
  delete m_staple;

  delete m_eta;
  delete m_zeta;

  delete m_Cud;
  delete m_Ut1;
  delete m_Ut2;
  delete m_Unu;

  delete m_force1;
  delete m_force2;

  delete m_eta2;
  delete m_eta3;
  delete m_zeta_mu;

  delete m_vt1;
  delete m_vt2;
  delete m_vt3;
  delete m_vt4;
}

//====================================================================
void Force_F_CloverTerm::set_parameters(const Parameters& params)
{
  set_parameters_impl(params);
  m_fopr_csw->set_parameters(params);
}


//====================================================================
void Force_F_CloverTerm::set_parameters_impl(const Parameters& params)
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
void Force_F_CloverTerm::set_parameters(const double kappa,
                                        const double cSW,
                                        const std::vector<int> bc)
{
  set_parameters_impl(kappa, cSW, bc);
  m_fopr_csw->set_parameters(m_kappa, m_cSW, m_boundary);
}

//====================================================================
void Force_F_CloverTerm::set_parameters_impl(const double kappa,
                                             const double cSW,
                                             const std::vector<int> bc)
{
#pragma omp barrier

  //- range check
  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_kappa = kappa;
    m_cSW   = cSW;
    m_boundary = bc;
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", m_cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }

}

//====================================================================
void Force_F_CloverTerm::get_parameters(Parameters& params) const
{
  params.set_double("hopping_parameter", m_kappa);
  params.set_double("clover_coefficient", m_cSW);
  params.set_int_vector("boundary_condition", m_boundary);
  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}

//====================================================================
void Force_F_CloverTerm::set_config(Field *U)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_csw->set_config(m_U);
  set_component();
}

//====================================================================
// override default force_core() in base class
void Force_F_CloverTerm::force_core(Field& force, const Field& eta)
{
  vout.crucial(m_vl, "Error at %s: force_core() is not available.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_CloverTerm::force_udiv(Field& force, const Field& eta)
{
  vout.crucial(m_vl, "Error at %s: force_udiv() is not available.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_CloverTerm::force_udiv1(Field& force_,
                                     const Field& zeta_, const Field& eta_)
{
#pragma omp barrier

  copy(*m_eta, eta_);
  copy(*m_zeta, zeta_);
#pragma omp barrier

  force_udiv1_impl(*m_force2, *m_zeta, *m_eta);

  copy(force_, *m_force2); // force_ = force;

#pragma omp barrier
}


//====================================================================
void Force_F_CloverTerm::force_udiv1_impl(Field_G& force,
                                          const Field_F& zeta,
                                          const Field_F& eta)
{
#pragma omp barrier

  m_force2->set(0.0);
#pragma omp barrier

  m_fopr_csw->mult_gm5(*m_eta2, eta);

  double fac = -m_kappa * m_cSW / 8.0;

  for (int mu = 0; mu < m_Ndim; ++mu) {
    for (int nu = 0; nu < m_Ndim; ++nu) {
      if (nu == mu) continue;

      m_fopr_csw->mult_isigma(*m_eta3, *m_eta2, mu, nu);

      copy(*m_Unu, 0, *m_U, nu);
#pragma omp barrier

      // R(1) and R(5)
      mult_Field_Gd(*m_vt1, 0, *m_Cud, index_dir(mu, nu), *m_eta3, 0);
      tensorProd_Field_F(*m_force1, zeta, *m_vt1);
      axpy(force, mu, 1.0, *m_force1, 0);
#pragma omp barrier

      // R(2)
      mult_Field_Gd(*m_vt3, 0, *m_U, mu, *m_eta3, 0);
      m_shift->backward(*m_vt2, zeta, nu);
      m_shift->backward(*m_vt1, *m_vt3, nu);
      m_shift->backward(*m_Ut1, *m_Unu, mu);
      mult_Field_Gn(*m_vt3, 0, *m_Ut1, 0, *m_vt1, 0);
      mult_Field_Gn(*m_vt4, 0, *m_U, nu, *m_vt2, 0);
      tensorProd_Field_F(*m_force1, *m_vt4, *m_vt3);
      axpy(force, mu, 1.0, *m_force1, 0);
#pragma omp barrier

      // R(4) and R(8)
      m_shift->backward(*m_vt1, *m_eta3, mu);
      m_shift->backward(*m_zeta_mu, zeta, mu);
      mult_Field_Gn(*m_vt4, 0, *m_Cud, index_dir(mu, nu), *m_zeta_mu, 0);
      tensorProd_Field_F(*m_force1, *m_vt4, *m_vt1);
      axpy(force, mu, 1.0, *m_force1, 0);
#pragma omp barrier

      // R(3)
      m_shift->backward(*m_vt1, *m_eta3, nu);
      mult_Field_Gn(*m_vt3, 0, *m_U, nu, *m_vt1, 0);
      m_shift->backward(*m_vt1, *m_vt3, mu);
      mult_Field_Gn(*m_vt4, 0, *m_U, mu, *m_zeta_mu, 0);
      m_shift->backward(*m_vt2, *m_vt4, nu);
      mult_Field_Gn(*m_vt4, 0, *m_U, nu, *m_vt2, 0);
      tensorProd_Field_F(*m_force1, *m_vt4, *m_vt1);
      axpy(force, mu, 1.0, *m_force1, 0);
#pragma omp barrier

      // R(6)
      m_shift->backward(*m_Ut1, *m_Unu, mu);
      mult_Field_Gdd(*m_Ut2, 0, *m_Ut1, 0, *m_U, mu);
      mult_Field_Gn(*m_vt1, 0, *m_Ut2, 0, *m_eta3, 0);
      mult_Field_Gd(*m_vt2, 0, *m_U, nu, zeta, 0);
      m_shift->forward(*m_vt3, *m_vt1, nu);
      m_shift->forward(*m_vt4, *m_vt2, nu);
      tensorProd_Field_F(*m_force1, *m_vt4, *m_vt3);
      axpy(force, mu, -1.0, *m_force1, 0);
#pragma omp barrier

      // R(7)
      mult_Field_Gd(*m_vt1, 0, *m_U, nu, *m_eta3, 0);
      mult_Field_Gn(*m_vt2, 0, *m_U, mu, *m_zeta_mu, 0);
      m_shift->backward(*m_vt3, *m_vt1, mu);
      m_shift->forward(*m_vt1, *m_vt3, nu);
      mult_Field_Gd(*m_vt4, 0, *m_U, nu, *m_vt2, 0);
      m_shift->forward(*m_vt2, *m_vt4, nu);
      tensorProd_Field_F(*m_force1, *m_vt2, *m_vt1);
      axpy(force, mu, -1.0, *m_force1, 0);
#pragma omp barrier

    } // nu-loop

    scal(force, mu, fac);
#pragma omp barrier

  } // mu-loop

}


//====================================================================
void Force_F_CloverTerm::set_component()
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


//============================================================END=====
