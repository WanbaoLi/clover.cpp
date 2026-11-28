/*!
        @file    action_F_Staggered_eo.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef ACTION_F_STAGGERED_EO_INCLUDED
#define ACTION_F_STAGGERED_EO_INCLUDED

#include "lib/Action/action.h"
#include "lib/Force/Fermion/force_F.h"
#include "lib/Field/index_eo.h"
#include "lib/Measurements/Fermion/fprop.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Action for staggered fermion: temporary implementation.

/*!
   This is Action class dedicated to staggered fermion in even-odd
   site index, as a temporary implementation.
   It may be not well organized and to be improved.
                                 [28 Dec 2011 H.Matsufuru]
    YAML is implemented.         [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                 [21 Mar 2015 Y.Namekawa]
    Multith-reading applied, and largely modified.
                                 [19 Jan 2026 H.Matsufuru]
 */

class Action_F_Staggered_eo : public Action {
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;
  std::string m_label;            //!< label of action

  Field_G *m_U;

  Fopr *m_fopr;
  Force *m_forceF;

  Fprop *m_fprop_MD;
  Fprop *m_fprop_H;

  Field m_psf;

  Index_eo m_index_eo;

 public:

  Action_F_Staggered_eo(
    Fopr *fopr, Force *forceF,
    Fprop *fprop_MD, Fprop *fprop_H)
    : Action(), m_fopr(fopr), m_forceF(forceF),
    m_fprop_MD(fprop_MD), m_fprop_H(fprop_H)
  { init(); }

  ~Action_F_Staggered_eo()
  { tidyup(); }

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  void set_config(Field *U);

  double langevin(RandomNumbers *);

  double calcH();

  void force(Field&);

 private:

  void init();

  void tidyup();

  double calcH_impl();

};
#endif
