/*!
        @file    force_F_Staggered.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_STAGGERED_INCLUDED
#define FORCE_F_STAGGERED_INCLUDED

#include "lib/Force/Fermion/force_F.h"
#include "lib/Fopr/fopr_Staggered.h"
#include "lib/Field/field_F_1spinor.h"
#include "lib/Field/shiftField_lex.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for even-odd staggered fermion operator.

/*!
    Force calculation for lexical staggered fermion operator.
                                      [20 Sep 2015 H.Matsufuru]
    Multi-threading applied.          [20 Jan 2026 H.Matsufuru]                
 */

class Force_F_Staggered : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim, m_Nvol;
  double m_mq;                   //!< quark mass.
  std::vector<int> m_boundary;   //!< boundary conditions.
  Bridge::VerboseLevel m_vl;

  Fopr_Staggered *m_fopr_ks;

  ShiftField_lex *m_shift;

  Field_G *m_force1;

  Field_F_1spinor *m_eta;
  Field_F_1spinor *m_zeta;
  Field_F_1spinor *m_eta2;

 public:

  DEPRECATED
  Force_F_Staggered() : Force() { init(); }

  Force_F_Staggered(const Parameters& params)
  : Force()
  { init(params); }

  ~Force_F_Staggered(){ tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(double mq, const std::vector<int> bc);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(Field& force, const Field& eta);

  void force_udiv1(Field&, const Field&, const Field&);

 private:

  void init();

  void init(const Parameters& params);

  void tidyup();

  void set_parameters_impl(const Parameters& params);

  void set_parameters_impl(double mq, const std::vector<int> bc);

  void force_udiv1_impl(Field_G& force, const Field_F_1spinor& zeta,
                        const Field_F_1spinor& eta);

};
#endif
