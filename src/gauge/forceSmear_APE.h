/*!
        @file    forceSmear_APE.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef FORCESMEAR_APE_INCLUDED
#define FORCESMEAR_APE_INCLUDED

#include "gauge/forceSmear.h"
#include "field/shiftField_lex.h"
#include "gauge/projection.h"
#include "gauge/smear_APE.h"
#include "io/bridgeIO.h"
using Bridge::vout;

//! Recursive calculation for APE smeared fermion force.

/*!
    Recursive calculation for APE smeared fermion force.
                                [08 Apr 2012 H.Matsufuru]
    YAML is implemented.        [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.    [10 Jun 2026 H.Matsufuru]
 */

class ForceSmear_APE : public ForceSmear
{
 public:
  static const std::string class_name;

 private:

  int m_Ndim, m_Nvol;         //!< spacetime dimension and volume
  std::vector<double> m_rho;
  Bridge::VerboseLevel m_vl;

  Projection *m_proj;

  ShiftField_lex *m_shift;

  std::vector<Field_G> m_U;
  std::vector<Field_G> m_iTheta;

  Field_G *m_v1;
  Field_G *m_v2;
  Field_G *m_vt1;
  Field_G *m_vt2;
  Field_G *m_vt3;

  Field_G *m_ct1;
  Field_G *m_ct2;
  Field_G *m_xi;

 public:

  ForceSmear_APE(Projection *proj, const Parameters& params)
    : m_proj(proj)
  { init(params); }

  ~ForceSmear_APE() { tidyup(); }

  // Setting parameters with Parameters object.
  void set_parameters(const Parameters& params);

  // Setting parameters with uniform smearing parameter.
  void set_parameters(const double rho1);

  // Setting parameters with anisotropic smearing parameter.
  void set_parameters(const std::vector<double>& rho);

  // Getting parameters by Parameters object.
  void get_parameters(Parameters& params) const;

  // Force computation.
  void force_udiv(Field_G& Sigma, const Field_G& Sigma_p, const Field_G& U);

 private:
  void init(const Parameters& params);

  void tidyup();

  int index(int mu, int nu){ return (mu + nu * m_Ndim); }

  void force_each(Field_G&,
                  const Field_G&, const Field_G&,
                  const Field_G&, const Field_G&, const int mu, const int nu);

  void staple(Field_G&,
              const Field_G&, const Field_G&,
              const int mu, const int nu);

#ifdef USE_FACTORY
 private:

  static ForceSmear *create_object_with_params(
                          Projection *proj, const Parameters& params)
  { return new ForceSmear_APE(proj, params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= ForceSmear::Factory_params::Register("APE",
                                           create_object_with_params);
    return init;
  }
#endif

};
#endif
