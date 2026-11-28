/*!
        @file    force_F_Clover_SF.h
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_CLOVER_SF_INCLUDED
#define FORCE_F_CLOVER_SF_INCLUDED

#include "lib/Force/Fermion/force_F_Wilson_SF.h"
#include "lib/Force/Fermion/tensorProd.h"
#include "lib/Field/field_SF.h"
#include "lib/Fopr/fopr_Clover_SF.h"
#include "lib/Measurements/Gauge/staple_SF.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for clover quark action with SF BC.

/*!
    Force calculation for clover quark action with Schroedinger
    functional boundary condition.
                                            [Y.Taniguchi]
    At present, only the Dirac representation for gamma-matrix
    is available.
                                [28 Dec 2011 H.Matsufuru]
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.    [10 Jan 2026 H.Matsufuru]
 */

class Force_F_Clover_SF : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;
  double m_kappa;                //!< hopping parameter
  double m_cSW;                  //!< clover coefficient
  std::vector<int> m_boundary;   //!< boundary conditions
  std::string m_repr;            //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;     //!< verbose level

  Fopr_Clover_SF *m_fopr_c;
  Force_F_Wilson_SF *m_force_w;

  Staple_SF *m_staple;

  ShiftField_lex *m_shift;

  std::vector<double> m_phi;    //!< SF boundary condition at t=0
  std::vector<double> m_phipr;  //!< SF boundary condition at t=Nt

  Mat_SU_N m_wk;    //!< SF boundary condition at t=0
  Mat_SU_N m_wkpr;  //!< SF boundary condition at t=Nt

  Field_F *m_eta;
  Field_F *m_zeta;

  Field_F *m_eta2;
  Field_F *m_eta3;
  Field_F *m_zeta1;
  Field_F *m_zeta2;

  Field_G *m_force2;

  Field_G *m_Cud;      //!< for force calculation
  Field_G *m_force1;

  Field_G *m_Ut1;
  Field_G *m_Ut2;
  Field_G *m_Umu;
  Field_G *m_Unu;

  Field_F *m_vt1;
  Field_F *m_vt2;
  Field_F *m_vt3;
  Field_F *m_vt4;

 public:

  DEPRECATED
  Force_F_Clover_SF()
    : m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  { init(); }

  Force_F_Clover_SF(const Parameters& params)
    : m_wk(CommonParameters::Nc()), m_wkpr(CommonParameters::Nc())
  { init(params); }

  ~Force_F_Clover_SF(){ tidyup(); }

  void set_parameters(const Parameters& params);

  //! Setting parameters of clover fermion.
  void set_parameters(const double kappa,
                      const double cSW,
                      const std::vector<int> bc,
                      const std::vector<double> phi,
                      const std::vector<double> phipr);

  void get_parameters(Parameters& params) const;

  //! Setting gauge configuration
  void set_config(Field *U);

  //! For recursive calculation of smeared force.
  void force_udiv(Field& force, const Field& eta);

  //! For recursive calculation of smeared force.
  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:

  void init();

  void init(const Parameters& params);

  void tidyup();

  void set_parameters_impl(const Parameters& params);

  //! Setting parameters of clover fermion.
  void set_parameters_impl(const double kappa,
                           const double cSW,
                           const std::vector<int> bc,
                           const std::vector<double> phi,
                           const std::vector<double> phipr);

  //! Core implemetation of clover force calculation.
  void force_udiv1_impl(Field_G& force,
                        const Field_F& zeta, const Field_F& eta);

  //! Set building components for force calculation.
  void set_component();

  void set_boundary_zero(Field_F& f);

  int index_dir(const int mu, const int nu)
  { return mu + m_Ndim * nu; }
};

#endif
