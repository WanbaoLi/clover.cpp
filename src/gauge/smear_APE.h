/*!
        @file    smear_APE.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef SMEAR_APE_INCLUDED
#define SMEAR_APE_INCLUDED

#include "gauge/smear.h"
#include "gauge/staple_lex.h"
#include "io/bridgeIO.h"
using Bridge::vout;


//! APE type smearing of link variables.

/*!
    This class smears link variables with APE-type construction
    of smeared links with a given projection operator to SU(N)
    group element.
    Parameter is \rho(\mu,\nu), which in general depends on
    the directions of the modified link and the staple.
    By explicitly giving \rho(\mu,\nu) as std::vector object,
    anisotropic setup is possible, while isotropic setup
    requires only one double parameter, `rho_uniform'.
                       [08 Apr 2012/15 Jul 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.       [07 Jan 2026 H.Matsufuru]
 */

class Smear_APE : public Smear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;                   //!< spacetime dimension
  int m_Nvol;                   //!< spacetime volume
  std::vector<double> m_rho;    //!< smearing parameter
  Bridge::VerboseLevel m_vl;

  Projection *m_proj;           //!< projector to group element

  Staple_lex *m_staple;

  Field_G *m_ct;
  Field_G *m_ut1;
  Field_G *m_ut2;

 public:
  DEPRECATED
  Smear_APE(Projection *proj)
    : m_proj(proj)
  { init(); }

  //! Constructor with parameters.
  Smear_APE(Projection *proj, const Parameters& params)
    : m_proj(proj)
  { init(params); }

  //! Deconstructor
  ~Smear_APE() { tidyup(); }

  //! Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  //! Setting parameter with isotropic parameter.
  void set_parameters(const double rho1);

  //! Setting parameter with anisotropic parameter.
  void set_parameters(const std::vector<double>& rho);

  //! Getting parameters by Parameters object.
  void get_parameters(Parameters& params) const;

  //! Smearing of a given gauge field.
  void smear(Field_G& Usmear, const Field_G& U);

 private:
  //! Initial setup by parameters.
  void init(const Parameters& params);

  //! Obsolete initial setup.
  void init();

  //! Final clean-up.
  void tidyup();

  int index(const int mu, const int nu)
  { return mu + m_Ndim * nu; }

#ifdef USE_FACTORY
 private:
  static Smear *create_object(Projection *proj)
  { return new Smear_APE(proj); }

  static Smear *create_object_with_params(Projection *proj,
                                          const Parameters& params)
  { return new Smear_APE(proj, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Smear::Factory::Register("APE", create_object);
    init &= Smear::Factory_params::Register("APE", create_object_with_params);
    return init;
  }
#endif
};
#endif
