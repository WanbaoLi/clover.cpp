/*!
        @file    director_Smear.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef DIRECTOR_SMEAR_INCLUDED
#define DIRECTOR_SMEAR_INCLUDED

#include <cassert>

#include "gauge/smear.h"
#include "runtime/director.h"
#include "field/field_G.h"
#include "io/bridgeIO.h"
using Bridge::vout;

//! Manager of smeared configurations.

/*!
    This director class handles smeared configurations.
                             [28 Dec 2011 H.Matsufuru]
    YAML is implemented.     [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                             [21 Mar 2015 Y.Namekawa]
    Multi-threading applied. [01 Feb 2026 H.Matsufuru]
 */

class Director_Smear : public Director
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  int m_Nsmear;                   //!< number of smearing to be applied
  Smear *m_smear;                 //!< smearing operator
  Field_G *m_U;                   //!< original thin link var.
  std::vector<Field_G> m_Usmear;  //!< smeared configs.
  int m_status_linkv;             //!< set to zero when link var. is updated

 public:
  //! constructor requires pointer to Smear object
  Director_Smear(Smear *smear)
    : m_smear(smear)
  { init(); }

  Director_Smear(Smear *smear, const Parameters& params)
    : m_smear(smear)
  { init(params); }

  //! set parameters, must be called before set_config
  void set_parameters(const Parameters& params);
  void set_parameters(const int Nsmear);

  //! get parameters
  void get_parameters(Parameters& params) const;

  //! get number of applied smearing operation
  int get_Nsmear() { return m_Nsmear; }

  //! get pointer to i-th smeared config (0th is original thin link)
  Field *getptr_smearedConfig(const int i_smear);

  Field_G *get_config();                  // smeared config
  Field_G *get_config(const int i_smear); // intermediate config

  //! set pointer to original thin link variable
  void set_config(Field *U);

  //! to be called when configuration is updated
  void notify_linkv() { m_status_linkv = 0; }

 private:

  void init(const Parameters& params);

  void init();

  //! smearing is performed by calling a function of Smear object
  void smear();
};
#endif
