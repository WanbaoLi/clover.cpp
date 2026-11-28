/*!
        @file    smear_HYP_SF.h
        @brief
        @author  Yusuke Tanigchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef SMEAR_HYP_SF_INCLUDED
#define SMEAR_HYP_SF_INCLUDED

#include "lib/Smear/smear.h"
#include "lib/Field/field_SF.h"
#include "lib/Field/shiftField_lex.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! HYP smearing of link variables with SF BC.

/*!
    HYP smearing of link variables with SF BC.
                                    [26 May 2012 Y.Taniguchi]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks.
                                    [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.       [09 Jan 2026 H.Matsufuru]
 */

class Smear_HYP_SF : public Smear
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Ndim, m_Nvol;
  double m_alpha1;              //!< HYP smearing parameters
  double m_alpha2;              //!< HYP smearing parameters
  double m_alpha3;              //!< HYP smearing parameters
  std::vector<double> m_phi;    //!< raw parameter values
  std::vector<double> m_phipr;  //!< raw parameter values

  Projection *m_proj;
  std::vector<Field_G> m_U;
  std::vector<Field_G> m_v1;
  std::vector<Field_G> m_v2;

  ShiftField_lex *m_shift;

  Mat_SU_N m_wk;        //!< SF boundary condition at t=0
  Mat_SU_N m_wkpr;      //!< SF boundary condition at t=Nt

  Field_G *m_c1;
  Field_G *m_c2;
  Field_G *m_c3;
  Field_G *m_ut;
  Field_G *m_vt1;
  Field_G *m_vt2;

 public:

  Smear_HYP_SF(Projection *proj, const Parameters& params)
    : m_proj(proj),
      m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  { init(params); }

  ~Smear_HYP_SF() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double alpha1,
                      const double alpha2,
                      const double alpha3,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void get_parameters(Parameters& params) const;

  void smear(Field_G& Usmear, const Field_G& U);

 private:

  void init(const Parameters& params);

  void tidyup();

  void staple(Field_G&, const Field_G&, const Field_G&,
              const int mu, const int nu);

  void step1();
  void step2();
  void step3(Field_G&);

  int index_v1(const int mu, const int nu, const int rho)
  {
    int sig = 6 - mu - nu - rho;
    if (sig > mu) --sig;
    return mu + m_Ndim * sig;
  }

  int index_v2(const int mu, int nu)
  {
    if (nu > mu) --nu;
    return mu + m_Ndim * nu;
  }

  int size_v1(){ return m_Ndim * (m_Ndim - 1); }

  int size_v2(){ return m_Ndim * (m_Ndim - 1); }


#ifdef USE_FACTORY
 private:

  static Smear *create_object_with_params(Projection *proj,
                                          const Parameters& params)
  { return new Smear_HYP_SF(proj, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Smear::Factory_params::Register("HYP_SF",
                                        create_object_with_params);
    return init;
  }
#endif

};
#endif
