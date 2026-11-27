/*!
        @file    force_F_CloverTerm.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_CLOVERTERM_INCLUDED
#define FORCE_F_CLOVERTERM_INCLUDED

#include "lib/Force/Fermion/force_F.h"
#include "lib/Force/Fermion/tensorProd.h"
#include "lib/Fopr/fopr_CloverTerm.h"
#include "lib/Field/shiftField_lex.h"
#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for clover term of clover fermion.

/*!
    This class calculates contribution clover term to the clover
    fermion action.
                                [28 Dec 2011 H.Matsufuru]
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.    [03 Jan 2026 H.Matsufuru]
 */

class Force_F_CloverTerm : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;                   //!< spacetime dimension
  double m_kappa;               //!< hopping parameter
  double m_cSW;                 //!< clover coefficient
  std::vector<int> m_boundary;  //!< boundary conditions
  Bridge::VerboseLevel m_vl;    //!< verbose level

  std::string m_repr;           //!< gamma matrix representation

  Fopr_CloverTerm *m_fopr_csw;  //!< fermion operator

  Field_G *m_Cud;               //!< for force calculation

  ShiftField_lex *m_shift;      //!< Field shifter
  Staple_lex *m_staple;         //!< Staple constructor

  Field_F *m_eta;
  Field_F *m_zeta;

  Field_G *m_force1;
  Field_G *m_force2;

  Field_G *m_Ut1;
  Field_G *m_Ut2;
  Field_G *m_Unu;

  Field_F *m_eta2;
  Field_F *m_eta3;
  Field_F *m_zeta_mu;

  Field_F *m_vt1;
  Field_F *m_vt2;
  Field_F *m_vt3;
  Field_F *m_vt4;

 public:
  //! Constructor
  DEPRECATED
  Force_F_CloverTerm()
  { init("Dirac");  /* default gamma matrix representation */ }

  DEPRECATED
  Force_F_CloverTerm(const std::string repr)
  { init(repr); }

  Force_F_CloverTerm(const Parameters& params)
  { init(params);  }

  //! Deconstructor
  ~Force_F_CloverTerm(){ tidyup(); }

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

  //! Force determination for clover fermion.
  void force_core(Field& force, const Field& eta);  // override default

  //! For recursive calculation of smeared force.
  void force_udiv(Field& force, const Field& eta);

  //! For recursive calculation of smeared force.
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:
  //! Obsolete initial setup.
  void init(const std::string);

  //! Initial setup with parameters.
  void init(const Parameters& params);

  //! final clean-up.
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

  //! Set building components for force calculation.
  void set_component();

  int index_dir(const int mu, const int nu)
  { return mu + m_Ndim * nu; }

};
#endif
