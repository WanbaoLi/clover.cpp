/*!
        @file    force_F_Wilson_TwistedMass_Nf2.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_TMFWILSON_NF2_INCLUDED
#define FORCE_F_TMFWILSON_NF2_INCLUDED

#include "lib/Force/Fermion/force_F.h"
#include "lib/Force/Fermion/tensorProd.h"
#include "lib/Fopr/fopr_Wilson_TwistedMass.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force for the twisted-mass Wilson fermion operator

/*!
    This class calculates the force of the twisted-mass Wilson
    fermion.
    The gamma matrix representation is given as control string
    "Dirac"(default) or "Chiral" at the construction, which is
    used to construct the Fopr_Wilson instance.
                                     [23 Dec 2011 H.Matusfuru]
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.         [12 Jan 2026 H.Matsufuru]
 */


class Force_F_Wilson_TwistedMass_Nf2 : public Force
{
 public:
  static const std::string class_name;

 private:
  double m_kappa;
  double m_tw_mass;
  std::vector<int> m_boundary;
  std::string m_repr;
  Bridge::VerboseLevel m_vl;

  Fopr_Wilson_TwistedMass *m_fopr_tw;

  Field_G *m_force1;

  Field_F *m_eta;
  Field_F *m_zeta;
  Field_F *m_eta2;
  Field_F *m_eta3;

 public:

  DEPRECATED
  Force_F_Wilson_TwistedMass_Nf2()
  { init("Dirac"); }

  DEPRECATED
  Force_F_Wilson_TwistedMass_Nf2(const std::string repr)
  { init(repr); }

  Force_F_Wilson_TwistedMass_Nf2(const Parameters& params)
  { init(params); }

  ~Force_F_Wilson_TwistedMass_Nf2(){ tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa, const double tw_mass,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(Field& force, const Field& eta);
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:

  void init(const Parameters& params);

  void init(const std::string repr);

  void tidyup();

  void set_parameters_impl(const Parameters& params);

  void set_parameters_impl(const double kappa,
                           const double tw_mass,
                           const std::vector<int> bc);

  void force_udiv1_impl(Field_G& force,
                        const Field_F& zeta, const Field_F& eta);

};
#endif
