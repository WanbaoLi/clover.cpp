/*!
        @file    action_F_Rational_SF.h
        @brief
        @author  Yusuke Taniguchi (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef ACTION_FRATIONAL_SF_INCLUDED
#define ACTION_FRATIONAL_SF_INCLUDED

#include "lib/Action/action.h"
#include "lib/Fopr/fopr_Rational.h"
#include "lib/Force/Fermion/force_F_Rational.h"
#include "lib/Field/field_SF.h"
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
    Multi-threading applied.             [01 Feb 2026 H.Matsufuru]
 */


class Action_F_Rational_SF : public Action {
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  std::string m_label;   // label of action

  Fopr *m_fopr_langev;
  Fopr *m_fopr_H;
  Force *m_forceF_MD;

  Field *m_U;

  Field m_psf;

 public:
  //! constructor requires pointers to Fopr and Force instances.
  Action_F_Rational_SF(Fopr *fopr_langev, Fopr *fopr_H,
                       Force *forceF_MD)
    : m_fopr_langev(fopr_langev), m_fopr_H(fopr_H),
      m_forceF_MD(forceF_MD)
  { init(); }

  Action_F_Rational_SF(Fopr *fopr_langev, Fopr *fopr_H,
                       Force *forceF_MD,
                       const Parameters& params)
    : m_fopr_langev(fopr_langev), m_fopr_H(fopr_H),
      m_forceF_MD(forceF_MD)
  { init(params); }

  //! destructor. constructed instances are deconstructed in tydyup().
  ~Action_F_Rational_SF()
  { tidyup(); }

  //! setting parameters and creating class instances.
  void set_parameters(const Parameters& params);

  //! getting parameters
  void get_parameters(Parameters& params) const;

  //! set the label of action.
  void set_label(const std::string label)
  {
    m_label = label;
    vout.detailed(m_vl, "  label: %s\n", m_label.c_str());
  }

  //! returns the label of action.
  std::string get_label() { return m_label; }

  //! setting gauge configuration.
  void set_config(Field *U);

  //! Langevin step called at the beginning of HMC.
  double langevin(RandomNumbers *);

  //! calculation of Hamiltonian.
  double calcH();

  //! returns force for updating conjugate momentum.
  void force(Field&);

 private:

  void init();

  void init(const Parameters& params);

  void tidyup();

};
#endif
