/*!
        @file    forceSmear_APE_SF.h
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCESMEAR_APE_ALT_INCLUDED
#define FORCESMEAR_APE_ALT_INCLUDED

#include "lib/Smear/forceSmear.h"
#include "lib/Field/field_SF.h"
#include "lib/Field/shiftField_lex.h"
#include "lib/Smear/smear_APE_SF.h"
#include "lib/Smear/projection.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Recursive calculation for APE smeared fermion force.

/*!
    Recursive calculation for APE smeared fermion force.
                                [08 Apr 2012 H.Matsufuru]
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.    [10 Jun 2026 H.Matsufuru]
 */

class ForceSmear_APE_SF : public ForceSmear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim, m_Nvol;
  std::vector<double> m_rho;
  std::vector<double> m_phi;   //!< SF boundary condition at t=0
  std::vector<double> m_phipr; //!< SF boundary condition at t=Nt
  Bridge::VerboseLevel m_vl;

  Projection *m_proj;

  ShiftField_lex *m_shift;

  std::vector<Field_G> m_U;
  std::vector<Field_G> m_iTheta;

  Mat_SU_N m_wk;    //!< SF boundary condition at t=0
  Mat_SU_N m_wkpr;  //!< SF boundary condition at t=Nt

  Field_G *m_v1;
  Field_G *m_v2;
  Field_G *m_vt1;
  Field_G *m_vt2;
  Field_G *m_vt3;

  Field_G *m_ct1;
  Field_G *m_ct2;
  Field_G *m_xi;

 public:

  ForceSmear_APE_SF(Projection *proj, const Parameters& params)
    : m_proj(proj),
      m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  { init(params); }

  ~ForceSmear_APE_SF(){ tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double rho1,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void set_parameters(const std::vector<double>& rho,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void get_parameters(Parameters& params) const;

  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

 private:

  void init(const Parameters& params);

  void tidyup();

  int index(int mu, int nu){ return (mu + nu * m_Ndim); }

  void force_each(Field_G&,
                  const Field_G&, const Field_G&,
                  const Field_G&, const Field_G&,
                  const int mu, const int nu);

  void staple(Field_G&,
              const Field_G&, const Field_G&,
              const int mu, const int nu);

#ifdef USE_FACTORY
 private:

  static ForceSmear *create_object_with_params(Projection *proj,
                                               const Parameters& params)
  { return new ForceSmear_APE_SF(proj, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= ForceSmear::Factory_params::Register("APE_SF",
                                             create_object_with_params);
    return init;
  }
#endif

};
#endif
