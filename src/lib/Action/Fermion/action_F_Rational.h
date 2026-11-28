/*!
        @file    action_F_Rational.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$

        @version $LastChangedRevision: 2693 $
*/

#ifndef ACTION_F_RATIONAL_INCLUDED
#define ACTION_F_RATIONAL_INCLUDED

#include "lib/Action/action.h"
#include "lib/Fopr/fopr_Rational.h"
#include "lib/Force/Fermion/force_F_Rational.h"
#include "lib/Tools/timer.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! action class for RHMC, with externally constructed Fopr_Rational.

/*!
    For the class, Fopr and Force objects are instantiated outside
    the class and specified at the construction.
    This class just provides the framework of rational actions.
                                        [28 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                        [21 Mar 2015 Y.Namekawa]
    Multithreading applied.             [04 Jan 2026 H.Matsufuru]
 */


class Action_F_Rational : public Action
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  std::string m_label;    // label of action

  Fopr *m_fopr_langev;
  Fopr *m_fopr_H;
  Force *m_fopr_force_MD;

  Field *m_U;

  Field m_psf;

  Timer m_timer;

 public:
  //! constructor requires pointers to Fopr and Force instances.
  Action_F_Rational(
    Fopr *fopr_langev, Fopr *fopr_H, Force *fopr_force_MD)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr_langev(fopr_langev), m_fopr_H(fopr_H),
    m_fopr_force_MD(fopr_force_MD)
  { }

  Action_F_Rational(Fopr *fopr_langev, Fopr *fopr_H, Force *fopr_force_MD,
                    const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_fopr_langev(fopr_langev), m_fopr_H(fopr_H),
    m_fopr_force_MD(fopr_force_MD)
  {
    set_parameters(params);
  }

  //! destructor. constructed instances are deconstructed in tydyup().
  ~Action_F_Rational() {}

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  //! set the label of action.
  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  //! returns the label of action.
  std::string get_label()
  {
    return m_label;
  }

  //! setting gauge configuration.
  void set_config(Field *U);

  //! Langevin step called at the beginning of HMC.
  double langevin(RandomNumbers *);

  //! calculation of Hamiltonian.
  double calcH();

  //! returns the force for updating conjugate momentum.
  //const Field force();
  void force(Field&);

 private:
};
#endif
