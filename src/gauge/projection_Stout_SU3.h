/*!
        @file    projection_Stout_SU3.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2026-02-05 20:40:12 #$
        @version $LastChangedRevision: 2693 $
*/

#ifndef PROJECTION_STOUT_SU3_INCLUDED
#define PROJECTION_STOUT_SU3_INCLUDED

#include "gauge/projection.h"
#include "io/bridgeIO.h"
using Bridge::vout;

//! Stout(exponential)-type projection to SU(Nc=3) gauge group.

/*!
    Stout(exponential)-type projection to SU(Nc) gauge group.
    Present implementation applies to SU(3) case only.
    The SU(3) properties are explicitly used.
                                    [08 Apr 2012 H.Matsufuru]
    Add a strict check for Nc, complement to assert.
                                    [31 May 2021 Y.Namekawa]
 */

class Projection_Stout_SU3 : public Projection
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

 public:
  
  Projection_Stout_SU3() { init(); }

  Projection_Stout_SU3(const Parameters& params) { init(params); }

  ~Projection_Stout_SU3() {}

  void set_parameters(const Parameters& param);

  void get_parameters(Parameters& param) const;

  //! projection U = P[alpha, C, Uorg]
  void project(Field_G& U,
               const double alpha,
               const Field_G& C, const Field_G& Uorg);

  //! determination of fields for force calculation
  void force_recursive(Field_G& Xi, Field_G& iTheta,
                       const double alpha, const Field_G& Sigmap,
                       const Field_G& C, const Field_G& U);

  // void print_stat();

 private:

  void init();

  void init(const Parameters& params);

  void exp_iQ(Field_G& e_iQ, const Field_G& iQ);
  void exp_iQ_bf(Field_G& e_iQ, const Field_G& iQ);

  void set_uw(double& u, double& w,
              const Mat_SU_N& iQ2, const Mat_SU_N& iQ3);

  void set_fj(dcomplex& f0, dcomplex& f1, dcomplex& f2,
              const double& u, const double& w);

  double func_xi0(const double w);
  double func_xi1(const double w);

#ifdef USE_FACTORY
 private:

  static Projection *create_object()
  { return new Projection_Stout_SU3(); }

  static Projection *create_object_with_params(const Parameters& params)
  { return new Projection_Stout_SU3(params); }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Projection::Factory::Register("Stout_SU3", create_object);
    init &= Projection::Factory_params::Register("Stout_SU3",
                                           create_object_with_params);
    return init;
  }
#endif
};
#endif
