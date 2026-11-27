/*!
        @file    force_F_Clover_Nf2_Isochemical.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_CLOVER_NF2_ISOCHEMICAL_INCLUDED
#define FORCE_F_CLOVER_NF2_ISOCHEMICAL_INCLUDED

#include "lib/Force/Fermion/force_F.h"
#include "lib/Fopr/fopr_Clover_Chemical.h"
#include "lib/Force/Fermion/force_F_Wilson_Nf2_Isochemical.h"
#include "lib/Force/Fermion/force_F_CloverTerm.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for clover quark action.

/*!
    At present, only the Dirac representation for gamma-matrix
    is available.                     [28 Dec 2011 H.Matsufuru]
    YAML is implemented.              [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.          [03 Jan 2026 H.Matsufuru]
 */

class Force_F_Clover_Nf2_Isochemical : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;                    //!< spacetime dimension
  double m_kappa;                //!< hopping parameter
  double m_cSW;                  //!< clover coefficient
  double m_mu;                   //!< Isospin chemical potential
  std::vector<int> m_boundary;   //!< boundary conditions
  Bridge::VerboseLevel m_vl;

  std::string m_repr;            //!< gamma matrix representation
  std::string m_mode;            //!< mult mode

  Fopr_Clover_Chemical *m_fopr_c;            //!< fermion operator
  Force_F_Wilson_Nf2_Isochemical *m_force_w; //!< Wilson fermion force
  Force_F_CloverTerm *m_force_csw;           //!< Clover term force

  Field_F *m_zeta;
  Field_F *m_eta;

  Field_G *m_force1;
  Field_G *m_force2;

 public:
  DEPRECATED
  Force_F_Clover_Nf2_Isochemical()
  { init("Dirac");  /* default gamma matrix representation */  }

  DEPRECATED
  Force_F_Clover_Nf2_Isochemical(const std::string repr)
  { init(repr); }

  Force_F_Clover_Nf2_Isochemical(const Parameters& params)
  { init(params); }

  //! Deconstructor
  ~Force_F_Clover_Nf2_Isochemical(){ tidyup(); }

  //! Setting parameters of clover fermion force.
  void set_parameters(const Parameters& params);

  //! Setting parameters with clover fermion parameters.
  void set_parameters(const double kappa,
                      const double cSW,
                      const double mu,
                      const std::vector<int> bc);

  //! Getting parameters of clover fermion force.
  void get_parameters(Parameters& params) const;

  void set_mode(const std::string& mode);

  //! Setting gauge configuration
  void set_config(Field *U);

  //! For recursive calculation of smeared force.
  void force_udiv(Field& force, const Field& eta);

  //! For recursive calculation of smeared force.
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);


 private:
  //! Obsolete initial setup.
  void init(const std::string repr);

  //! Initial setup with parameters.
  void init(const Parameters& params);

  //! Final clean-up.
  void tidyup();

  void set_parameters_impl(const Parameters& params);

  void set_parameters_impl(const double kappa,
                           const double cSW,
                           const double mu,
                           const std::vector<int> bc);

  //! Core implemetation of clover force calculation.
  void force_udiv1_impl(Field_G& force,
                        const Field_F& zeta, const Field_F& eta);

};
#endif
