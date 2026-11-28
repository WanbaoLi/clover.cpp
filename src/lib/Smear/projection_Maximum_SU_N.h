/*!
        @file    projection_Maximum_SU_N.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef PROJECTION_MAXIMUM_SU_N_INCLUDED
#define PROJECTION_MAXIMUM_SU_N_INCLUDED

#include "projection.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Maximum projection to SU(N) gauge group.

/*!
    Maximum projection for SU(N) matrix by Cabibbo-Marinari
    method, SU(2) subgroup transformation.
    The code was originally written by Takashi Umeda (1997) in
    Fortran by explicitely assuming SU(3) group.
    Genralization to SU((N) was done by H.M.
                                    [09 Aug 2012 H.Matsufuru]
    YAML is implemented.            [14 Nov 2012 Y.Namekawa]
    Multi-threading applied.        [09 Jan 2026 H.Matsufuru]
 */


class Projection_Maximum_SU_N : public Projection
{
 public:
  static const std::string class_name;

 private:
  int m_Niter;         //!< maximum iteration of maximization steps
  double m_Enorm;      //!< convergence criterion of maximization

  Bridge::VerboseLevel m_vl; //!< verbose level

  Field_G *m_At;
  Field_G *m_Dlt;
  Field_G *m_Ut;
  Field_G *m_vt;
  Field_G *m_wt;

 public:

  Projection_Maximum_SU_N(const Parameters& params)
  { init(params); }

  ~Projection_Maximum_SU_N() {tidyup(); }

  void set_parameters(const Parameters& params);

  void set_parameters(const int Niter, const double Enorm);

  void get_parameters(Parameters& params) const;

  //! projection U = P[alpha, C, Uorg]
  void project(Field_G& U,
               const double alpha,
               const Field_G& C, const Field_G& Uorg);

  //! force calculation: invalid in this class.
  void force_recursive(Field_G& Xi, Field_G& iTheta,
                       const double alpha, const Field_G& Sigmap,
                       const Field_G& C, const Field_G& U);

 private:

  void init(const Parameters& params);

  void tidyup();

  void reset_field(const Field_G& Uref);

  //! maximization of ReTr[U^\dag V].
  void maxTr(Field_G& U, const Field_G& V);

  //! maximization by SU(2) subgroup.
  void maxTr_SU2(const int, const int, Field_G&, Field_G&, Field_G&);

  //! matrix index for convenience.
  int mindex(const int i, const int j, const int Nc)
  { return i + j * Nc; }

  
#ifdef USE_FACTORY
 private:

  static Projection *create_object_with_params(const Parameters& params)
  { return new Projection_Maximum_SU_N(params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Projection::Factory_params::Register("Maximum_SU_N", create_object_with_params);
    return init;
  }
#endif
};
#endif
