/*!
        @file    force_F_Clover_Nf2.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_CLOVER_NF2_INCLUDED
#define FORCE_F_CLOVER_NF2_INCLUDED

#include "operator/fopr_Clover.h"
#include "force/force_F_Wilson_Nf2.h"
#include "force/force_F_CloverTerm.h"
#include "gauge/staple_lex.h"
#include "io/bridgeIO.h"
using Bridge::vout;

//! Force calculation for clover quark action.

/*!
    At present, only the Dirac representation for gamma-matrix
    is available.                     [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.              [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.          [03 Jan 2026 H.Matsufuru]
 */

class Force_F_Clover_Nf2 : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;                      //!< spacetime dimension
  double m_kappa;                  //!< hopping parameter
  double m_cSW;                    //!< clover coefficient
  std::vector<int> m_boundary;     //!< boundary conditions
  std::string m_repr;              //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;       //!< verbose level

  Fopr_Clover *m_fopr_c;           //!< fermion operator
  Force_F_Wilson_Nf2 *m_force_w;   //!< Wilson fermion force
  Force_F_CloverTerm *m_force_csw; //!< Clover term force

  Staple_lex *m_staple;

  Field_G *m_force1;
  Field_G *m_force2;
  Field_G *m_Cmu_ud;

  Field_F *m_zeta;
  Field_F *m_eta;

 public:
  DEPRECATED
  Force_F_Clover_Nf2()
  { init("Dirac");  /* default gamma matrix representation */ }

  DEPRECATED
  Force_F_Clover_Nf2(const std::string repr)
  { init(repr); }

  //! Constructor with parameters.
  Force_F_Clover_Nf2(const Parameters& params)
  { init(params); }

  //! Deconstructor
  ~Force_F_Clover_Nf2()
  { tidyup(); }

  //! Setting parameters of clover fermion force.
  void set_parameters(const Parameters& params);

  //! Setting parameters with clover fermion parameters.
  void set_parameters(const double kappa,
                      const double cSW,
                      const std::vector<int> bc);

  //! Getting parameters of clover fermion force.
  void get_parameters(Parameters& params) const;

  //! Setting gauge configuration
  void set_config(Field *U);

  //! For recursive calculation of smeared force.
  void force_udiv(Field& force, const Field& eta);

  //! For recursive calculation of smeared force.
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  //! Obsolete intial setup.
  void init(const std::string);

  //! Intial setup with parameters.
  void init(const Parameters& params);

  //! Final clean-up.
  void tidyup();

  //! Setting parameters of clover fermion force.
  void set_parameters_impl(const Parameters& params);

  //! Setting parameters with clover fermion parameters.
  void set_parameters_impl(const double kappa,
                           const double cSW,
                           const std::vector<int> bc);

  //! Core implemetation of clover force calculation.
  void force_udiv1_impl(Field_G& force,
                        const Field_F& zeta, const Field_F& eta);

};
#endif
