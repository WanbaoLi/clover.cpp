/*!
        @file    fopr_Clover_SF.cpp
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Fopr/fopr_Clover_SF.h"
#include "lib/Field/field_SF.h"
#include "lib/Fopr/fopr_thread-inc.h"


#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Fopr_Clover_SF::register_factory();
}
#endif

const std::string Fopr_Clover_SF::class_name = "Fopr_Clover_SF";

//====================================================================
namespace {
  inline double mult_uv_r(const double *g, const double *w)
  {
    return   g[0] * w[0] - g[1] * w[1]
           + g[2] * w[2] - g[3] * w[3]
           + g[4] * w[4] - g[5] * w[5];
  }


  inline double mult_uv_i(const double *g, const double *w)
  {
    return   g[0] * w[1] + g[1] * w[0]
           + g[2] * w[3] + g[3] * w[2]
           + g[4] * w[5] + g[5] * w[4];
  }
}

//====================================================================
void Fopr_Clover_SF::init(const Parameters& params)
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
    m_repr = "Dirac";  // default
    vout.general(m_vl, "gamma_matrix_type is not given: defalt = %s\n",
                 m_repr.c_str());
  }

  if(m_repr != "Dirac"){
    vout.crucial("%s: only Dirac repr. is available.\n");
    exit(EXIT_FAILURE);
  }

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;

  m_boundary.resize(m_Ndim);

  set_parameters_impl(params);

  setup_sg();

  m_U = 0;

  m_fopr_w = new Fopr_Wilson_SF(params);

  m_staple = new Staple_SF(params);

  m_shift = new ShiftField_lex();

  m_w1.reset(m_NinF, m_Nvol, 1);
  m_w2.reset(m_NinF, m_Nvol, 1);

  m_v1 = new Field_G(m_Nvol, 1);
  m_v2 = new Field_G(m_Nvol, 1);
  m_v3 = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Fopr_Clover_SF::init(const std::string repr)
{
  ThreadManager::assert_single_thread(class_name);

  m_vl = CommonParameters::Vlevel();

  vout.general(m_vl, "%s: construction (obsolete)\n", class_name.c_str());
  vout.increase_indent();

  m_repr = repr;
  if(m_repr != "Dirac"){
    vout.crucial("%s: only Dirac repr. is accepted/\n");
    exit(EXIT_FAILURE);
  }

  m_Nvol = CommonParameters::Nvol();
  m_Ndim = CommonParameters::Ndim();
  m_Nc   = CommonParameters::Nc();
  m_Nd   = CommonParameters::Nd();
  m_NinF = 2 * m_Nc * m_Nd;

  m_boundary.resize(m_Ndim);

  setup_sg();

  m_U = 0;

  m_fopr_w = new Fopr_Wilson_SF();

  m_staple = new Staple_SF();

  m_shift = new ShiftField_lex();

  m_w1.reset(m_NinF, m_Nvol, 1);
  m_w2.reset(m_NinF, m_Nvol, 1);

  m_v1 = new Field_G(m_Nvol, 1);
  m_v2 = new Field_G(m_Nvol, 1);
  m_v3 = new Field_G(m_Nvol, 1);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}


//====================================================================
void Fopr_Clover_SF::tidyup()
{
  delete m_fopr_w;
  delete m_staple;
  delete m_shift;

}


//====================================================================
void Fopr_Clover_SF::setup_sg()
{
  m_GM.resize(m_Ndim + 1);
  m_SG.resize(m_Ndim * m_Ndim);

  GammaMatrixSet *gmset = GammaMatrixSet::New(m_repr);

  m_GM[0] = gmset->get_GM(gmset->GAMMA1);
  m_GM[1] = gmset->get_GM(gmset->GAMMA2);
  m_GM[2] = gmset->get_GM(gmset->GAMMA3);
  m_GM[3] = gmset->get_GM(gmset->GAMMA4);
  m_GM[4] = gmset->get_GM(gmset->GAMMA5);

  m_SG[sg_index(0, 1)] = gmset->get_GM(gmset->SIGMA12);
  m_SG[sg_index(1, 2)] = gmset->get_GM(gmset->SIGMA23);
  m_SG[sg_index(2, 0)] = gmset->get_GM(gmset->SIGMA31);
  m_SG[sg_index(3, 0)] = gmset->get_GM(gmset->SIGMA41);
  m_SG[sg_index(3, 1)] = gmset->get_GM(gmset->SIGMA42);
  m_SG[sg_index(3, 2)] = gmset->get_GM(gmset->SIGMA43);

  m_SG[sg_index(1, 0)] = m_SG[sg_index(0, 1)].mult(-1);
  m_SG[sg_index(2, 1)] = m_SG[sg_index(1, 2)].mult(-1);
  m_SG[sg_index(0, 2)] = m_SG[sg_index(2, 0)].mult(-1);
  m_SG[sg_index(0, 3)] = m_SG[sg_index(3, 0)].mult(-1);
  m_SG[sg_index(1, 3)] = m_SG[sg_index(3, 1)].mult(-1);
  m_SG[sg_index(2, 3)] = m_SG[sg_index(3, 2)].mult(-1);

  m_SG[sg_index(0, 0)] = gmset->get_GM(gmset->UNITY);
  m_SG[sg_index(1, 1)] = gmset->get_GM(gmset->UNITY);
  m_SG[sg_index(2, 2)] = gmset->get_GM(gmset->UNITY);
  m_SG[sg_index(3, 3)] = gmset->get_GM(gmset->UNITY);
  // these 4 gamma matrices are actually not used.

  delete gmset;

}


//====================================================================
void Fopr_Clover_SF::set_parameters(const Parameters& params)
{
  set_parameters_impl(params);

  m_fopr_w->set_parameters(params);

  //m_staple->set_parameters(m_phi, m_phipr);
  m_staple->set_parameters(params);

}


//====================================================================
void Fopr_Clover_SF::set_parameters_impl(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double              kappa, cSW;
  std::vector<int>    bc;
  std::vector<double> phi, phipr;

  int err          = 0;
  int err_optional = 0;
  err          += params.fetch_double("hopping_parameter", kappa);
  err          += params.fetch_double("clover_coefficient", cSW);
  err          += params.fetch_int_vector("boundary_condition", bc);
  err          += params.fetch_double_vector("phi", phi);
  err          += params.fetch_double_vector("phipr", phipr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters_impl(kappa, cSW, bc, phi, phipr);

}


//====================================================================
void Fopr_Clover_SF::set_parameters(const double kappa,
                                    const double cSW,
                                    const std::vector<int> bc,
                                    const std::vector<double> phi,
                                    const std::vector<double> phipr)
{

  set_parameters_impl(kappa, cSW, bc, phi, phipr);

  //- propagate parameters
  Parameters params;
  get_parameters(params);

  m_fopr_w->set_parameters(params);

  m_staple->set_parameters(params);

}


//====================================================================
void Fopr_Clover_SF::set_parameters_impl(const double kappa,
                                         const double cSW,
                                         const std::vector<int> bc,
                                         const std::vector<double> phi,
                                         const std::vector<double> phipr)
{
#pragma omp barrier

  assert(bc.size() == m_Ndim);

  int ith = ThreadManager::get_thread_id();

  //- store values
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
  }
#pragma omp barrier

  //- print input parameters
  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
  vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
  vout.general(m_vl, "  cSW   = %12.8f\n", m_cSW);
  for (int mu = 0; mu < m_Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
  }
  vout.general(m_vl, "  phi1   = %12.8f\n", m_phi[0]);
  vout.general(m_vl, "  phi2   = %12.8f\n", m_phi[1]);
  vout.general(m_vl, "  phi3   = %12.8f\n", m_phi[2]);
  vout.general(m_vl, "  phipr1 = %12.8f\n", m_phipr[0]);
  vout.general(m_vl, "  phipr2 = %12.8f\n", m_phipr[1]);
  vout.general(m_vl, "  phipr3 = %12.8f\n", m_phipr[2]);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_SF::get_parameters(Parameters& params) const
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
void Fopr_Clover_SF::set_config(Field *U)
{
#pragma omp barrier

  // at this moment, this method is not multi-threaded yet.
  //  ThreadManager::assert_single_thread(class_name);

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_U = (Field_G *)U;
#pragma omp barrier

  m_fopr_w->set_config(U);
  set_csw();

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_SF::set_mode(const std::string mode)
{
#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) m_mode = mode;

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_SF::mult(Field& v, const Field& f)
{
  if (m_mode == "D") {
    D(v, f);
  } else if (m_mode == "DdagD") {
    DdagD(v, f);
  } else if (m_mode == "Ddag") {
    Ddag(v, f);
  } else if (m_mode == "H") {
    H(v, f);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Clover_SF::mult_dag(Field& v, const Field& f)
{
  if (m_mode == "D") {
    Ddag(v, f);
  } else if (m_mode == "DdagD") {
    DdagD(v, f);
  } else if (m_mode == "Ddag") {
    D(v, f);
  } else if (m_mode == "H") {
    H(v, f);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), m_mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Clover_SF::mult(Field& v, const Field& f,
                          const std::string mode)
{
  if (mode == "D") {
    D(v, f);
  } else if (mode == "DdagD") {
    DdagD(v, f);
  } else if (mode == "Ddag") {
    Ddag(v, f);
  } else if (mode == "H") {
    H(v, f);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Clover_SF::mult_dag(Field& v, const Field& f,
                              const std::string mode)
{
  if (mode == "D") {
    Ddag(v, f);
  } else if (mode == "DdagD") {
    DdagD(v, f);
  } else if (mode == "Ddag") {
    D(v, f);
  } else if (mode == "H") {
    H(v, f);
  } else {
    vout.crucial(m_vl, "Error at %s: undefined mode = %s.\n",
                 class_name.c_str(), mode.c_str());
    exit(EXIT_FAILURE);
  }
}


//====================================================================
void Fopr_Clover_SF::DdagD(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  D(m_w2, f);
  mult_gm5(w, m_w2);
  D(m_w2, w);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Clover_SF::Ddag(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  mult_gm5(w, f);
  D(m_w2, w);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Clover_SF::H(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  D(m_w2, f);
  mult_gm5(w, m_w2);
}


//====================================================================
void Fopr_Clover_SF::D(Field& w, const Field& f)
{
  assert(f.nex() == 1);

  m_fopr_w->D(w, f);

  mult_csw(m_w1, f);

  axpy(w, -1.0, m_w1);
#pragma omp barrier

  Field_SF::set_boundary_zero(w);

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_SF::mult_isigma(Field_F& v, const Field_F& w,
                                 const int mu, const int nu)
{
  assert(mu != nu);

  mult_iGM(v, m_SG[sg_index(mu, nu)], w);
}


//====================================================================
void Fopr_Clover_SF::mult_csw(Field& v, const Field& w)
{
  mult_csw_dirac(v, w);
}


//====================================================================
void Fopr_Clover_SF::mult_csw_dirac(Field& v, const Field& w)
{
#pragma omp barrier

  assert(w.nex() == 1);

  const int Nvc = 2 * m_Nc;
  const int Ndf = 2 * m_Nc * m_Nc;

  const int id1 = 0;
  const int id2 = Nvc;
  const int id3 = Nvc * 2;
  const int id4 = Nvc * 3;

  const double *w2 = w.ptr(0);
  double       *v2 = v.ptr(0);

  double *Bx = m_Bx.ptr(0);
  double *By = m_By.ptr(0);
  double *Bz = m_Bz.ptr(0);
  double *Ex = m_Ex.ptr(0);
  double *Ey = m_Ey.ptr(0);
  double *Ez = m_Ez.ptr(0);

  v.set(0.0);
#pragma omp barrier

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  for (int site = is; site < ns; ++site) {
    int iv = Nvc * m_Nd * site;
    int ig = Ndf * site;

    for (int ic = 0; ic < m_Nc; ++ic) {
      int ic_r = 2 * ic;
      int ic_i = 2 * ic + 1;
      int ic_g = ic * Nvc + ig;

      // isigma_23 * Bx
      v2[ic_r + id1 + iv] -= mult_uv_i(&Bx[ic_g], &w2[id2 + iv]);
      v2[ic_i + id1 + iv] += mult_uv_r(&Bx[ic_g], &w2[id2 + iv]);
      v2[ic_r + id2 + iv] -= mult_uv_i(&Bx[ic_g], &w2[id1 + iv]);
      v2[ic_i + id2 + iv] += mult_uv_r(&Bx[ic_g], &w2[id1 + iv]);

      v2[ic_r + id3 + iv] -= mult_uv_i(&Bx[ic_g], &w2[id4 + iv]);
      v2[ic_i + id3 + iv] += mult_uv_r(&Bx[ic_g], &w2[id4 + iv]);
      v2[ic_r + id4 + iv] -= mult_uv_i(&Bx[ic_g], &w2[id3 + iv]);
      v2[ic_i + id4 + iv] += mult_uv_r(&Bx[ic_g], &w2[id3 + iv]);

      // isigma_31 * By
      v2[ic_r + id1 + iv] += mult_uv_r(&By[ic_g], &w2[id2 + iv]);
      v2[ic_i + id1 + iv] += mult_uv_i(&By[ic_g], &w2[id2 + iv]);
      v2[ic_r + id2 + iv] -= mult_uv_r(&By[ic_g], &w2[id1 + iv]);
      v2[ic_i + id2 + iv] -= mult_uv_i(&By[ic_g], &w2[id1 + iv]);

      v2[ic_r + id3 + iv] += mult_uv_r(&By[ic_g], &w2[id4 + iv]);
      v2[ic_i + id3 + iv] += mult_uv_i(&By[ic_g], &w2[id4 + iv]);
      v2[ic_r + id4 + iv] -= mult_uv_r(&By[ic_g], &w2[id3 + iv]);
      v2[ic_i + id4 + iv] -= mult_uv_i(&By[ic_g], &w2[id3 + iv]);

      // isigma_12 * Bz
      v2[ic_r + id1 + iv] -= mult_uv_i(&Bz[ic_g], &w2[id1 + iv]);
      v2[ic_i + id1 + iv] += mult_uv_r(&Bz[ic_g], &w2[id1 + iv]);
      v2[ic_r + id2 + iv] += mult_uv_i(&Bz[ic_g], &w2[id2 + iv]);
      v2[ic_i + id2 + iv] -= mult_uv_r(&Bz[ic_g], &w2[id2 + iv]);

      v2[ic_r + id3 + iv] -= mult_uv_i(&Bz[ic_g], &w2[id3 + iv]);
      v2[ic_i + id3 + iv] += mult_uv_r(&Bz[ic_g], &w2[id3 + iv]);
      v2[ic_r + id4 + iv] += mult_uv_i(&Bz[ic_g], &w2[id4 + iv]);
      v2[ic_i + id4 + iv] -= mult_uv_r(&Bz[ic_g], &w2[id4 + iv]);

      // isigma_41 * Ex
      v2[ic_r + id1 + iv] += mult_uv_i(&Ex[ic_g], &w2[id4 + iv]);
      v2[ic_i + id1 + iv] -= mult_uv_r(&Ex[ic_g], &w2[id4 + iv]);
      v2[ic_r + id2 + iv] += mult_uv_i(&Ex[ic_g], &w2[id3 + iv]);
      v2[ic_i + id2 + iv] -= mult_uv_r(&Ex[ic_g], &w2[id3 + iv]);

      v2[ic_r + id3 + iv] += mult_uv_i(&Ex[ic_g], &w2[id2 + iv]);
      v2[ic_i + id3 + iv] -= mult_uv_r(&Ex[ic_g], &w2[id2 + iv]);
      v2[ic_r + id4 + iv] += mult_uv_i(&Ex[ic_g], &w2[id1 + iv]);
      v2[ic_i + id4 + iv] -= mult_uv_r(&Ex[ic_g], &w2[id1 + iv]);

      // isigma_42 * Ey
      v2[ic_r + id1 + iv] -= mult_uv_r(&Ey[ic_g], &w2[id4 + iv]);
      v2[ic_i + id1 + iv] -= mult_uv_i(&Ey[ic_g], &w2[id4 + iv]);
      v2[ic_r + id2 + iv] += mult_uv_r(&Ey[ic_g], &w2[id3 + iv]);
      v2[ic_i + id2 + iv] += mult_uv_i(&Ey[ic_g], &w2[id3 + iv]);

      v2[ic_r + id3 + iv] -= mult_uv_r(&Ey[ic_g], &w2[id2 + iv]);
      v2[ic_i + id3 + iv] -= mult_uv_i(&Ey[ic_g], &w2[id2 + iv]);
      v2[ic_r + id4 + iv] += mult_uv_r(&Ey[ic_g], &w2[id1 + iv]);
      v2[ic_i + id4 + iv] += mult_uv_i(&Ey[ic_g], &w2[id1 + iv]);

      // isigma_43 * Ez
      v2[ic_r + id1 + iv] += mult_uv_i(&Ez[ic_g], &w2[id3 + iv]);
      v2[ic_i + id1 + iv] -= mult_uv_r(&Ez[ic_g], &w2[id3 + iv]);
      v2[ic_r + id2 + iv] -= mult_uv_i(&Ez[ic_g], &w2[id4 + iv]);
      v2[ic_i + id2 + iv] += mult_uv_r(&Ez[ic_g], &w2[id4 + iv]);

      v2[ic_r + id3 + iv] += mult_uv_i(&Ez[ic_g], &w2[id1 + iv]);
      v2[ic_i + id3 + iv] -= mult_uv_r(&Ez[ic_g], &w2[id1 + iv]);
      v2[ic_r + id4 + iv] -= mult_uv_i(&Ez[ic_g], &w2[id2 + iv]);
      v2[ic_i + id4 + iv] += mult_uv_r(&Ez[ic_g], &w2[id2 + iv]);
    }
  }
#pragma omp barrier

  scal(v, m_kappa * m_cSW); // v *= m_kappa * m_cSW;

#pragma omp barrier
}


//====================================================================
void Fopr_Clover_SF::set_csw()
{
  set_fieldstrength(m_Bx, 1, 2);
  set_fieldstrength(m_By, 2, 0);
  set_fieldstrength(m_Bz, 0, 1);
  set_fieldstrength(m_Ex, 3, 0);
  set_fieldstrength(m_Ey, 3, 1);
  set_fieldstrength(m_Ez, 3, 2);
}


/*!
  The field strength defined by clover with the SF BC.
  <ul>
  <li>The field strength is set to zero at t=0 boundary.
  <li>This is performed automatically by a use of Staple_SF.
  </ul>
*/
//====================================================================
void Fopr_Clover_SF::set_fieldstrength(Field_G& Fst,
                                       const int mu, const int nu)
{
#pragma omp barrier

  m_staple->upper(*m_v2, *m_U, mu, nu);

  m_staple->lower(*m_v3, *m_U, mu, nu);

  mult_Field_Gnd(Fst, 0, *m_U, mu, *m_v2, 0);
  multadd_Field_Gnd(Fst, 0, *m_U, mu, *m_v3, 0, -1.0);

  mult_Field_Gdn(*m_v1, 0, *m_v2, 0, *m_U, mu);
  multadd_Field_Gdn(*m_v1, 0, *m_v3, 0, *m_U, mu, -1.0);
#pragma omp barrier

  m_shift->forward(*m_v2, *m_v1, mu);

  axpy(Fst, 1.0, *m_v2);

  ah_Field_G(Fst, 0);
  scal(Fst, 0.25);
#pragma omp barrier

}


//====================================================================
double Fopr_Clover_SF::flop_count()
{
  //- Counting of floating point operations in giga unit.
  //  not implemented, yet.

  vout.general(m_vl, "%s: flop_count() has not been implemented.\n",
               class_name.c_str());

  const double gflop = 0.0;

  return gflop;
}


//============================================================END=====
