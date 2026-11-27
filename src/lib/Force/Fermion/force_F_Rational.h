/*!
        @file    force_F_Rational.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_RATIONAL_INCLUDED
#define FORCE_F_RATIONAL_INCLUDED

#include "lib/Force/Fermion/force_F.h"
#include "lib/Solver/shiftsolver_CG.h"
#include "lib/Fopr/fopr.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for smeared fermion operators.

/*!
    This class determines the force of a rational approximation
    for a given fermion operator.
                                    [28 Dec 2011 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
 */

class Force_F_Rational : public Force
{
 public:
  static const std::string class_name;

 private:
  int m_Np;                  //!< number of poles in rational approx.
  int m_n_exp;               //!< numerator of the exponent
  int m_d_exp;               //!< denominator of the exponent
  double m_x_min;            //!< lower bound of approximate sign function
  double m_x_max;            //!< upper bound of approximate sign function
  int m_Niter;               //!< maximum iteration of shiftsolver
  double m_Stop_cond;        //!< stopping condition of shift solver
  Bridge::VerboseLevel m_vl; //!< verbose level

  Fopr *m_fopr;              //!< kernel fermion operator
  Force *m_force;            //!< kernel fermion force

  Shiftsolver_CG *m_solver;

  double m_a0;               //!< rational approx. coefficients
  std::vector<double> m_bl;  //!< rational approx. coefficients
  std::vector<double> m_cl;  //!< rational approx. coefficients

  std::vector<Field> m_psi;
  Field *m_eta;

  Field_G *m_force1;
  Field_G *m_force2;

 public:

  DEPRECATED
  Force_F_Rational(Fopr *fopr, Force *force)
    : m_fopr(fopr), m_force(force)
  { init(); }

  Force_F_Rational(Fopr *fopr, Force *force, const Parameters& params)
    : m_fopr(fopr), m_force(force)
  { init(params); }

  ~Force_F_Rational() { tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const int Np, const int n_exp, const int d_exp,
                      const double x_min, const double x_max,
                      const int Niter, const double Stop_cond);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  void force_udiv(Field&, const Field&);

  void force_core1(Field&, const Field&, const Field&);  // dummy entry
  void force_udiv1(Field&, const Field&, const Field&);  // dummy entry

 private:
  //! Initial setup with parameters
  void init(const Parameters& params);

  //! Obsolete initial setup with parameters
  void init();

  //! Finial clean-up.
  void tidyup();

  void setup(int Np_prev);

  void set_parameters_impl(const Parameters& params);

  void set_parameters_impl(const int Np,
                           const int n_exp, const int d_exp,
                           const double x_min, const double x_max,
                           const int Niter, const double Stop_cond);

  void force_udiv_impl(Field_G&, const Field&);

};
#endif
