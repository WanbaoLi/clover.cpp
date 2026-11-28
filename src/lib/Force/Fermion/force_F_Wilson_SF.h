/*!
        @file    force_F_Wilson_SF.h
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_WILSON_SF_INCLUDED
#define FORCE_F_WILSON_SF_INCLUDED

#include "lib/Force/Fermion/force_F.h"
#include "lib/Force/Fermion/tensorProd.h"
#include "lib/Field/field_SF.h"
#include "lib/Fopr/fopr_Wilson_SF.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force for the Wilson quark action with SF BC.

/*!
    Force calculation for clover quark action with Schroedinger
    functional boundary condition.              [Y.Taniguchi]
                                     [23 Dec 2011 H.Matusfuru]
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.         [11 Jan 2026 H.Matsufuru]
 */

class Force_F_Wilson_SF : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;
  double m_kappa;
  std::vector<int> m_boundary;
  std::string m_repr;
  Bridge::VerboseLevel m_vl;

  Fopr_Wilson_SF *m_fopr_w;

  Field_G *m_force1;

  Field_F *m_eta;
  Field_F *m_zeta;

  Field_F *m_eta2;
  Field_F *m_eta3;

public:

  DEPRECATED
  Force_F_Wilson_SF()
  { init(); }

  Force_F_Wilson_SF(const Parameters& params)
  { init(params); }

  ~Force_F_Wilson_SF(){ tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double kappa, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(Field& force, const Field& eta);

  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:

  void init();

  void init(const Parameters& params);

  void tidyup();

  void set_parameters_impl(const Parameters& params);

  void set_parameters_impl(const double kappa, const std::vector<int> bc);

  void force_udiv1_impl(Field_G& force,
                        const Field_F& zeta, const Field_F& eta);

};
#endif
