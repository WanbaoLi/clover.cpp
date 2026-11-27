/*!
        @file    smear_APE_SF.h
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef SMEAR_APE_SF_INCLUDED
#define SMEAR_APE_SF_INCLUDED

#include "smear.h"

#include "lib/Measurements/Gauge/staple_SF.h"
#include "lib/Field/field_G.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! APE smearing of link variables for Schroedinger functional BC.

/*!
    This class is alternative to the Smear_APE class.
                                   [08 Apr 2012 H.Matsufuru]
    APE smearing of link variables for Schroedinger functional BC.
                                   [26 May 2012 Y.Taniguchi]                  
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks.
                                    [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.       [08 Jan 2026 H.Matsufuru]
 */

class Smear_APE_SF : public Smear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;
  int m_Nvol;
  std::vector<double> m_rho;
  std::vector<double> m_phi;   //!< SF boundary condition at t=0
  std::vector<double> m_phipr; //!< SF boundary condition at t=Nt
  Bridge::VerboseLevel m_vl;

  Projection *m_proj;

  Mat_SU_N m_wk;               //!< SF boundary condition at t=0
  Mat_SU_N m_wkpr;             //!< SF boundary condition at t=Nt

  Staple_SF *m_staple;

  Field_G *m_ct;
  Field_G *m_ut1;
  Field_G *m_ut2;

 public:

  Smear_APE_SF(Projection *proj, const Parameters& params)
    : m_proj(proj),
      m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  { init(params); }

  ~Smear_APE_SF() {}

  void set_parameters(const Parameters& params);

  void set_parameters(const double rho1,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void set_parameters(const std::vector<double>& rho,
                      const std::vector<double>& phi,
                      const std::vector<double>& phipr);

  void get_parameters(Parameters& params) const;

  void smear(Field_G& Usmear, const Field_G& U);

 private:
  void init(const Parameters& params);

  void tidyup();

  int index(const int mu, const int nu)
  { return mu + m_Ndim * nu; }

#ifdef USE_FACTORY
 private:

  static Smear *create_object_with_params(Projection *proj,
                                          const Parameters& params)
  { return new Smear_APE_SF(proj, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Smear::Factory_params::Register("APE_SF",
                                            create_object_with_params);
    return init;
  }
#endif
};

#endif
