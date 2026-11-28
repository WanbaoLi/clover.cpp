/*!
        @file    force_F_Staggered_eo.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_STAGGERED_EO_INCLUDED
#define FORCE_F_STAGGERED_EO_INCLUDED

#include "lib/Force/Fermion/force_F.h"
#include "lib/Fopr/fopr_Staggered_eo.h"
#include "lib/Field/field_F_1spinor.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for even-odd staggered fermion operator.

/*!
    Force calculation for even-odd staggered fermion operator.
                                      [28 Dec 2011 H.Matsufuru]
    YAML is implemented.              [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.          [20 Jan 2026 H.Matsufuru]
 */

class Force_F_Staggered_eo : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Nvol, m_Nvol2, m_Ndim;
  double m_mq;
  std::vector<int> m_boundary;
  Bridge::VerboseLevel m_vl;

  Fopr_Staggered_eo *m_fopr_ks;

  Index_eo *m_index_eo;

  ShiftField_eo *m_shift_eo;

  Field_G *m_Ueo;

  Field_G *m_force_e;
  Field_G *m_force_o;
  Field_G *m_force_eo;

  Field_F_1spinor *m_eta;
  Field_F_1spinor *m_zeta;
  Field_F_1spinor *m_eta2;


 public:

  DEPRECATED
  Force_F_Staggered_eo() : Force() { init(); }

  Force_F_Staggered_eo(const Parameters& params)
    : Force()
  { init(params); }

  ~Force_F_Staggered_eo() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double mq, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_core(Field& force, const Field& eta);
  void force_udiv(Field& force, const Field& eta);

  void force_core1(Field&, const Field&, const Field&);  // dummy entry
  void force_udiv1(Field&, const Field&, const Field&);  // dummy entry

  void force_udiv1(Field_G& force,
                   const Field_F_1spinor& zeta, const Field_F_1spinor& eta,
                   int ieo);

 private:
  //! Obsolete initial setup.
  void init();

  //! Initial setup with parameters.
  void init(const Parameters& params);

  //! final clean-up.
  void tidyup();

  void set_parameters_impl(const Parameters& params);

  void set_parameters_impl(const double mq,
                           const std::vector<int> bc);

};
#endif
