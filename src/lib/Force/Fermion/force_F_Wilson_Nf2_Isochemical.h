/*!
        @file    force_F_Wilson_Nf2_Isochemical.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_WILSON_NF2_ISOCHEMICAL_INCLUDED
#define FORCE_F_WILSON_NF2_ISOCHEMICAL_INCLUDED

#include <string>
#include <vector>

#include "lib/Force/Fermion/force_F.h"
#include "lib/Fopr/fopr_Wilson_Chemical.h"
#include "lib/Field/field_F.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force for the Wilson fermion operator with isospin chemical potential.

/*!
    This class calculates the force of the standard Wilson fermion with
    isospin chemical potential with two flavors.
                                             [24 Aug 2011 H.Matusfuru]
    YAML is implemented.                     [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.                 [05 Jan 2026 H.Matsufuru]
 */

class Force_F_Wilson_Nf2_Isochemical : public Force
{
 public:
  static const std::string class_name;

 private:
  int    m_Ndim;
  double m_kappa;                      //!< hopping parameter
  double m_mu;                         //!< isospin chemical potential
  double m_exp_mu;                     //!< exp(mu)
  std::vector<int> m_boundary;
  std::string m_repr;
  Bridge::VerboseLevel m_vl;

  Fopr_Wilson_Chemical *m_fopr_w;

  std::string m_mode;

  Field_F *m_eta;
  Field_F *m_zeta;
  Field_F *m_eta2;
  Field_F *m_eta3;

  Field_G *m_force1;

 public:
  DEPRECATED
  Force_F_Wilson_Nf2_Isochemical()
  { init("Dirac"); /* default gamma-matrix representation */ }

  DEPRECATED
  Force_F_Wilson_Nf2_Isochemical(const std::string repr)
  { init(repr); }

  Force_F_Wilson_Nf2_Isochemical(const Parameters& params)
  { init(params); }

  ~Force_F_Wilson_Nf2_Isochemical(){ tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa,
                      const double mu,
                      const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void set_mode(const std::string& mode);

  void force_udiv(Field& force, const Field& eta);
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  //! Obsolete initial setup.
  void init(const Parameters& params);

  //! Initial setup with parameters.
  void init(const std::string repr);

  //! Final clean-up.
  void tidyup();

  void set_parameters_impl(const Parameters& params);

  void set_parameters_impl(const double kappa,
                           const double mu,
                           const std::vector<int> bc);

  void force_udiv1_impl(Field_G& force,
                        const Field_F& zeta, const Field_F& eta);
};
#endif
