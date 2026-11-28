/*!
        @file    smear_APE_spatial.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef SMEAR_APE_SPATIAL_INCLUDED
#define SMEAR_APE_SPATIAL_INCLUDED

#include "lib/Smear/smear.h"
#include "lib/Measurements/Gauge/staple_lex.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

//! APE type smearing of spatial link variables.

/*!
    This class smears spatial link variables with APE-type
    construction of smeared links with a given projection
    operator to SU(N) group element.
    Parameter is \rho, which specifies the mixing rate
    of original thin link and staples.
                               [09 Aug 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.       [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                               [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.   [09 Jan 2026 H.Matsufuru]
 */

class Smear_APE_spatial : public Smear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;                 //!< spacetime dimension
  int m_Nvol;                 //!< spacetime volume
  double m_rho;               //!< smearing parameter
  Bridge::VerboseLevel m_vl;  //!< verbose level

  Projection *m_proj;  //!< projector to group element.

  Staple_lex *m_staple;

  Field_G *m_ct;
  Field_G *m_ut1;
  Field_G *m_ut2;

 public:

  Smear_APE_spatial(Projection *proj, const Parameters& params)
    : m_proj(proj)
  { init(params); }

  //! Deconstructor
  ~Smear_APE_spatial() {tidyup(); }

  //! Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  //! Setting smearing parameter.
  void set_parameters(const double rho);

  //! Getting parameters by Parameters object.
  void get_parameters(Parameters& params) const;

  //! Smearing of a given gauge field.
  void smear(Field_G& Usmear, const Field_G& U);

 private:

  void init(const Parameters& params);

  void tidyup();

#ifdef USE_FACTORY
 private:

  static Smear *create_object_with_params(Projection *proj,
                                          const Parameters& params)
  { return new Smear_APE_spatial(proj, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Smear::Factory_params::Register("APE_spatial",
                                            create_object_with_params);
    return init;
  }
#endif

};
#endif
