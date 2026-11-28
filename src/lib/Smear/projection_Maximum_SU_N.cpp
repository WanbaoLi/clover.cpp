/*!
        @file    projection_Maximum_SU_N.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#include "lib/Smear/projection_Maximum_SU_N.h"
#include "lib/ResourceManager/threadManager.h"
#include "lib/Field/field_thread-inc.h"

#ifdef USE_FACTORY_AUTOREGISTER
namespace {
  bool init = Projection_Maximum_SU_N::register_factory();
}
#endif

const std::string Projection_Maximum_SU_N::class_name
                                      = "Projection_Maximum_SU_N";

//====================================================================
void Projection_Maximum_SU_N::init(const Parameters& params)
{
  ThreadManager::assert_single_thread(class_name);

  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  } else {
    m_vl = CommonParameters::Vlevel();
  }

  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  vout.general(m_vl, "%s: construction\n", class_name.c_str());
  vout.increase_indent();

  set_parameters(params);

  // to be modified Ndim -> 1.
  m_At  = new Field_G(Nvol, Ndim);
  m_Dlt = new Field_G(Nvol, Ndim);
  m_Ut  = new Field_G(Nvol, Ndim);
  m_vt  = new Field_G(Nvol, Ndim);
  m_wt  = new Field_G(Nvol, Ndim);

  vout.decrease_indent();
  vout.general(m_vl, "%s: construction finished.\n",
               class_name.c_str());

}

//====================================================================
void Projection_Maximum_SU_N::tidyup()
{
  delete m_At;
  delete m_Dlt;
  delete m_Ut;
  delete m_vt;
  delete m_wt;
}

//====================================================================
void Projection_Maximum_SU_N::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  int    Niter;
  double Enorm;

  int err = 0;
  err += params.fetch_int("maximum_number_of_iteration", Niter);
  err += params.fetch_double("convergence_criterion", Enorm);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  set_parameters(Niter, Enorm);
}


//====================================================================
void Projection_Maximum_SU_N::set_parameters(const int Niter,
                                             const double Enorm)
{
  //- range check
  int err = 0;
  err += ParameterCheck::non_negative(Niter);
  err += ParameterCheck::non_zero(Enorm);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    m_Niter = Niter;
    m_Enorm = Enorm;
  }

  //- print input parameters
  vout.general(m_vl, "%s: parameters\n", class_name.c_str());
  vout.general(m_vl, "  Niter  = %d\n", m_Niter);
  vout.general(m_vl, "  Enorm  = %12.4e\n", m_Enorm);

}


//====================================================================
void Projection_Maximum_SU_N::get_parameters(Parameters& params) const
{
  params.set_int("maximum_number_of_iteration", m_Niter);
  params.set_double("convergence_criterion", m_Enorm);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Projection_Maximum_SU_N::reset_field(const Field_G& Uref)
{
#pragma omp barrier

  const int Nex  = Uref.nex();
  const int Nvol = Uref.nvol();

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    delete m_At;
    delete m_Dlt;
    delete m_Ut;
    delete m_vt;
    delete m_wt;

    m_At  = new Field_G(Nvol, Nex);
    m_Dlt = new Field_G(Nvol, Nex);
    m_Ut  = new Field_G(Nvol, Nex);
    m_vt  = new Field_G(Nvol, Nex);
    m_wt  = new Field_G(Nvol, Nex);
  }
#pragma omp barrier

  vout.general(m_vl, "working fields reset.\n");

}

//====================================================================
void Projection_Maximum_SU_N::project(Field_G& U,
                                      const double alpha,
                                      const Field_G& Cst,
                                      const Field_G& Uorg)
{
  const int Nex  = Uorg.nex();
  const int Nvol = Uorg.nvol();
  const int Nc   = CommonParameters::Nc();

  assert(Cst.nex()  == Nex);
  assert(Cst.nvol() == Nvol);
  assert(U.nex()  == Nex);
  assert(U.nvol() == Nvol);

  if(Nex != m_Ut->nex() || Nvol != m_Ut->nvol()){
    reset_field(Uorg);
  }

  for (int ex = 0; ex < Nex; ++ex) {
    copy(*m_Ut, ex, Cst, ex);
    axpy(*m_Ut, ex, 1.0 - alpha, Uorg, ex);
  }
#pragma omp barrier

  maxTr(U, *m_Ut);

}


//====================================================================
void Projection_Maximum_SU_N::force_recursive(Field_G& Xi,
                                              Field_G& iTheta,
                                              const double alpha,
                                              const Field_G& Sigmap,
                                              const Field_G& Cst,
                                              const Field_G& Uorg)
{
  vout.crucial(m_vl, "Error at %s: force_recursive() is not available.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Projection_Maximum_SU_N::maxTr(Field_G& G0, const Field_G& Cst)
{
  const int Nc   = CommonParameters::Nc();
  const int Nvol = Cst.nvol();
  const int Nex  = Cst.nex();

  assert(Nvol == G0.nvol());
  assert(Nex  == G0.nex());

  const int Nmt = 1;  // number of subgroup maximization loop:
                      // seems not important because of outer iter-loop.

  Mat_SU_N unity(Nc);
  unity.unit();

  for (int ex = 0; ex < Nex; ++ex) {
    copy(*m_At, ex, Cst, ex);
  }
#pragma omp barrier

  vout.detailed(m_vl, "Maximum projection start.\n");

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      G0.set_mat(site, ex, unity);
    }
  }
#pragma omp barrier

  for (int iter = 0; iter < m_Niter; ++iter) {

    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        m_Dlt->set_mat(site, ex, unity);
      }
    }
#pragma omp barrier

    for (int imt = 0; imt < Nmt; ++imt) {
      for (int i1 = 0; i1 < Nc; ++i1) {
        int i2 = (i1 + 1) % Nc;
        maxTr_SU2(i1, i2, G0, *m_At, *m_Dlt);
      }
    }

    //- convergence test
    double retr1 = 0.0;
    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        for (int cc = 0; cc < Nc; ++cc) {
          retr1 += m_Dlt->cmp_r(cc * (1 + Nc), site, ex);
        }
      }
    }

    double retr = retr1;
    ThreadManager::reduce_sum_global(retr, ith, nth);

    const int Npe    = Communicator::size();
    double    deltaV = 1.0 - retr / (Nc * Nvol * Nex * Npe);
    vout.detailed(m_vl, "  iter = %d  deltaV = %12.4e\n", iter, deltaV);

    if (deltaV < m_Enorm) {
      for (int ex = 0; ex < Nex; ++ex) {
        for (int site = is; site < ns; ++site) {
          Mat_SU_N ut(Nc);
          G0.mat_dag(ut, site, ex);
          G0.set_mat(site, ex, ut);
        }
      }
#pragma omp barrier

      vout.detailed(m_vl, "Maximum projection converged.\n");

      return;
    }
  }

  vout.crucial(m_vl, "Error at %s: Maximum projection not converged.\n",
               class_name.c_str());
  exit(EXIT_FAILURE);
}


//====================================================================
void Projection_Maximum_SU_N::maxTr_SU2(const int i1, const int i2,
                                        Field_G& Gmax,
                                        Field_G& A, Field_G& Udelta)
{
  const int Nc   = CommonParameters::Nc();
  const int Nvol = A.nvol();
  const int Nex  = A.nex();

  if(Nex > m_vt->nex()){
    vout.crucial(m_vl, "too large matrix size in %s\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(i1 < Nc);
  assert(i2 < Nc);

  const int j1 = mindex(i1, i1, Nc);
  const int j2 = mindex(i2, i2, Nc);
  const int k1 = mindex(i1, i2, Nc);
  const int k2 = mindex(i2, i1, Nc);

  //----------[     | # # 0 | <i1  ]--------------------------
  //----------[ V = | # # 0 | <i2  ]--------------------------
  //----------[     | 0 0 1 |      ]--------------------------

  //Field_G v(Nvol, Nex);
  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {

      Mat_SU_N at(Nc);
      at = A.mat(site, ex);

      double xlamd =
          at.r(j1) * at.r(j1) + at.i(j1) * at.i(j1) + 2.0 * at.r(j1) * at.r(j2)
        + at.r(k1) * at.r(k1) + at.i(k1) * at.i(k1) - 2.0 * at.i(j1) * at.i(j2)
        + at.r(k2) * at.r(k2) + at.i(k2) * at.i(k2) - 2.0 * at.r(k1) * at.r(k2)
        + at.r(j2) * at.r(j2) + at.i(j2) * at.i(j2) + 2.0 * at.i(k1) * at.i(k2);
      xlamd = 1.0 / sqrt(xlamd);

      Mat_SU_N vt(Nc);
      vt.unit();
      vt.set(j1, (at.r(j1) + at.r(j2)) * xlamd, (-at.i(j1) + at.i(j2)) * xlamd);
      vt.set(k1, (at.r(k2) - at.r(k1)) * xlamd, (-at.i(k2) - at.i(k1)) * xlamd);
      vt.set(k2, (at.r(k1) - at.r(k2)) * xlamd, (-at.i(k1) - at.i(k2)) * xlamd);
      vt.set(j2, (at.r(j1) + at.r(j2)) * xlamd, (at.i(j1) - at.i(j2)) * xlamd);

      m_vt->set_mat(site, ex, vt);
    }
  }
#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    mult_Field_Gnn(*m_wt, ex, A, ex, *m_vt, ex);
    copy(A, ex, *m_wt, ex);
  }
#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    mult_Field_Gnn(*m_wt, ex, Gmax, ex, *m_vt, ex);
    copy(Gmax, ex, *m_wt, ex);
  }
#pragma omp barrier

  for (int ex = 0; ex < Nex; ++ex) {
    mult_Field_Gnn(*m_wt, ex, Udelta, ex, *m_vt, ex);
    copy(Udelta, ex, *m_wt, ex);
  }
#pragma omp barrier

}

//============================================================END=====
