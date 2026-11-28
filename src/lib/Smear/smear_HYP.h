/*!
        @file    smear_HYP.h
        @brief
        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef SMEAR_HYP_INCLUDED
#define SMEAR_HYP_INCLUDED

#include "smear.h"

#include "Field/shiftField_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! HYP smearing of link variables.

/*!
    HYP smearing of link variables.
    Ported from the Fortran code.
                                    [18 Mar 2012 H.Matsufuru]
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    unique_ptr is introduced to avoid memory leaks
                                    [21 Mar 2015 Y.Namekawa]
    Multi-threading applied.       [07 Jan 2026 H.Matsufuru]
 */


class Smear_HYP : public Smear
{
 public:
  static const std::string class_name;

 private:
  int m_Ndim;         //!< spacetime dimension
  int m_Nvol;         //!< spacetime volume
  double m_alpha1;    //!< HYP smearing parameter (level-1)
  double m_alpha2;    //!< HYP smearing parameter (level-2)
  double m_alpha3;    //!< HYP smearing parameter (level-3)
  Bridge::VerboseLevel m_vl;

  Projection *m_proj;

  ShiftField_lex *m_shift;

  std::vector<Field_G> m_U;
  std::vector<Field_G> m_v1;
  std::vector<Field_G> m_v2;

  Field_G *m_c1;
  Field_G *m_c2;
  Field_G *m_c3;
  Field_G *m_ut;
  Field_G *m_vt1;
  Field_G *m_vt2;

 public:
  Smear_HYP(Projection *proj)
    : m_proj(proj)
  { init(); }

  Smear_HYP(Projection *proj, const Parameters& params)
    : m_proj(proj)
  { init(params); }

  ~Smear_HYP() {tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const double alpha1,
                      const double alpha2,
                      const double alpha3);

  void get_parameters(Parameters& params) const;

  void smear(Field_G& Usmear, const Field_G& U);

 private:
  void init();

  void init(const Parameters& params);

  void tidyup();

  void staple(Field_G&, const Field_G&, const Field_G&,
              const int mu, const int nu);

  void step1();
  void step2();
  void step3(Field_G&);

  int index_v1(const int mu, const int nu, const int rho)
  {
    int sig = 6 - mu - nu - rho;
    if (sig > mu) --sig;
    return mu + m_Ndim * sig;
  }

  int index_v2(const int mu, int nu)
  {
    if (nu > mu) --nu;
    return mu + m_Ndim * nu;
  }

  int size_v1(){ return m_Ndim * (m_Ndim - 1); }

  int size_v2(){ return m_Ndim * (m_Ndim - 1); }

#ifdef USE_FACTORY
 private:
  static Smear *create_object(Projection *proj)
  {
    return new Smear_HYP(proj);
  }

  static Smear *create_object_with_params(Projection *proj, const Parameters& params)
  {
    return new Smear_HYP(proj, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Smear::Factory::Register("HYP", create_object);
    init &= Smear::Factory_params::Register("HYP", create_object_with_params);
    return init;
  }
#endif
};
#endif
