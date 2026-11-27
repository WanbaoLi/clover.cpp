/*!
        @file    force_F_Smeared.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCE_F_SMEARED_INCLUDED
#define FORCE_F_SMEARED_INCLUDED

#include "lib/Force/Fermion/force_F.h"
#include "lib/Smear/forceSmear.h"
#include "lib/Fopr/fopr.h"
#include "lib/Smear/director_Smear.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! Force calculation for smeared fermion operators.

/*!
    This class determines the force of smeared fermion operator
    using smearing director (MultiSmear instance) and base
    fermion force instance.
                                      [28 Dec 2011 H.Matsufuru]
    Modified: set_mode() is added to incorporate non-hermitian H
                                      [21 Jan 2012 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                                      [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.          [11 Jan 2026 H.Matsufuru]
*/

class Force_F_Smeared : public Force
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Force *m_force;
  ForceSmear *m_force_smear;
  Director_Smear *m_director_smear;

  Field_G *m_force2;
  Field_G *m_force1;

 public:
  Force_F_Smeared(Force *force, ForceSmear *force_smear,
                  Director_Smear *director_smear)
    : m_force(force), m_force_smear(force_smear),
      m_director_smear(director_smear)
  { init(); }

  Force_F_Smeared(Force *force, ForceSmear *force_smear,
                  Director_Smear *director_smear,
                  const Parameters& params)
    : m_force(force), m_force_smear(force_smear),
      m_director_smear(director_smear)
  { init(params); }

  ~Force_F_Smeared(){ tidyup(); }

  void set_parameters(const Parameters&);

  void get_parameters(Parameters&) const;

  void set_config(Field *U);

  void set_mode(const std::string& mode)
  { m_force->set_mode(mode); }

  void force_udiv(Field& force, const Field& eta);

  void force_udiv1(Field& force, const Field& zeta, const Field& eta);

 private:

  void init(const Parameters& params);

  void init();

  void tidyup();

  void mult_jacobian(Field_G& force);

};
#endif
