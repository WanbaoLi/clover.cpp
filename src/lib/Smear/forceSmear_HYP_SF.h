/*!
        @file    forceSmear_HYP_SF.h
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCESMEAR_HYP_SF_INCLUDED
#define FORCESMEAR_HYP_SF_INCLUDED

#include "lib/Smear/forceSmear.h"
#include "lib/Field/field_SF.h"
#include "lib/Field/shiftField_lex.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Recursive calculation of HYP smeared fermion force with SF BC.

/*!
    Recursive calculation of HYP smeared fermion force with SF BC.
                                [26 May 2012 Y.Taniguchi]
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.    [10 Jan 2026 H.Matsufuru]
 */

class ForceSmear_HYP_SF : public ForceSmear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;         //!< spacetime dimension
  int m_Nvol;         //!< spacetime volume
  double m_alpha1;    //!< HYP smearing parameter (level-1)
  double m_alpha2;    //!< HYP smearing parameter (level-2)
  double m_alpha3;    //!< HYP smearing parameter (level-3)
  Bridge::VerboseLevel m_vl;

  Projection *m_proj;

  ShiftField_lex *m_shift;

  std::vector<Field_G> m_U;
  std::vector<Field_G> m_v1, m_v2;
  std::vector<Field_G> m_Sigma3, m_Sigma2;
  std::vector<Field_G> m_iTheta3, m_iTheta2, m_iTheta1;

  std::vector<double> m_phi;   //!< SF boundary condition at t=0
  std::vector<double> m_phipr; //!< SF boundary condition at t=Nt

  Mat_SU_N m_wk;               //!< SF boundary condition at t=0
  Mat_SU_N m_wkpr;             //!< SF boundary condition at t=Nt

  Field_G *m_vt1;
  Field_G *m_vt2;
  Field_G *m_vt3;
  Field_G *m_ct1;
  Field_G *m_ct2;
  Field_G *m_ct3;

 public:

  ForceSmear_HYP_SF(Projection *proj, const Parameters& params)
    : m_proj(proj),
      m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  { init(params);
  }

  ~ForceSmear_HYP_SF(){ };

  void set_parameters(const Parameters& params);
  void set_parameters(const double alpha1, const double alpha2,
                      const double alpha3,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void get_parameters(Parameters& params) const;

  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

 private:
  void init(const Parameters& params);

  void tidyup();

  void force_step3(Field_G&, const Field_G&);
  void force_step2(Field_G&);
  void force_step1(Field_G&);

  void force_each(Field_G&, const Field_G&, const Field_G&,
                  const Field_G&, const Field_G&, const int mu, const int nu);

  void smear_step1();
  void smear_step2();

  void staple(Field_G&, const Field_G&, const Field_G&,
              const int mu, const int nu);

  int idx1(const int mu, const int nu, const int rho)
  {
    int sig = 6 - mu - nu - rho;
    if (sig > mu) --sig;
    return mu + m_Ndim * sig;
  }

  int idx1b(const int mu, int nu, int rho)
  {
    if (nu > mu) --nu;
    if (rho > mu) --rho;
    if (rho > nu) --rho;
    return mu + m_Ndim * (nu + (m_Ndim - 1) * rho);
  }

  int idx2(const int mu, int nu)
  {
    if (nu > mu) --nu;
    return mu + m_Ndim * nu;
  }

  int size1(){ return m_Ndim * (m_Ndim - 1); }

  int size1b(){ return m_Ndim * (m_Ndim - 1) * (m_Ndim - 2); }

  int size2(){ return m_Ndim * (m_Ndim - 1); }

#ifdef USE_FACTORY
 private:

  static ForceSmear *create_object_with_params(Projection *proj,
                                               const Parameters& params)
  { return new ForceSmear_HYP_SF(proj, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= ForceSmear::Factory_params::Register("HYP_SF",
                                            create_object_with_params);
    return init;
  }
#endif

};
#endif
