/*!
      @file    bridgeQXS_Wilson.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef BRIDGEQXS_WILSON_INCLUDED
#define BRIDGEQXS_WILSON_INCLUDED

namespace BridgeQXS {
  // Reconstructed EO kernels take the local boundary multiplier (one on
  // nonterminal ranks). Each SIMD tile packs two complex rows per column.
  enum { NDF12 = 12 };

  // real_t = double

  void mult_wilson_bulk_dirac(double *v2, double *up, double *v1,
                              double kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_wilson_1_dirac(double *buf_xp, double *buf_xm,
                           double *buf_yp, double *buf_ym,
                           double *buf_zp, double *buf_zm,
                           double *buf_tp, double *buf_tm,
                           double *up, double *v1,
                           int *bc, int *Nsize, int *do_comm);

  void mult_wilson_2_dirac(double *v2, double *up, double *v1,
                           double *buf_xp, double *buf_xm,
                           double *buf_yp, double *buf_ym,
                           double *buf_zp, double *buf_zm,
                           double *buf_tp, double *buf_tm,
                           double kappa, int *bc,
                           int *Nsize, int *do_comm);

  void mult_wilson_eo_bulk_dirac(double *v2, double *up, double *v1,
                                 double *xp,
                                 double kappa, int *bc,
                                 int *Nsize, int *do_comm, int *Leo,
                                 const int ieo, const int iflag);

  void mult_wilson_eo_bulk_recon12_dirac(double *v2, double *up, double *v1,
                                         double *xp,
                                         double kappa, int *bc,
                                         int *Nsize, int *do_comm, int *Leo,
                                         const int ieo, const int iflag);

  void mult_wilson_eo_1_dirac(double *buf_xp, double *buf_xm,
                              double *buf_yp, double *buf_ym,
                              double *buf_zp, double *buf_zm,
                              double *buf_tp, double *buf_tm,
                              double *up, double *v1, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_1_recon12_dirac(double *buf_xp, double *buf_xm,
                                      double *buf_yp, double *buf_ym,
                                      double *buf_zp, double *buf_zm,
                                      double *buf_tp, double *buf_tm,
                                      double *up, double *v1, int *bc,
                                      int *Nsize, int *do_comm, int *Leo,
                                      const int ieo, const int iflag);

  void mult_wilson_eo_2_dirac(double *v2, double *up, double *v1,
                              double *xp,
                              double *buf_xp, double *buf_xm,
                              double *buf_yp, double *buf_ym,
                              double *buf_zp, double *buf_zm,
                              double *buf_tp, double *buf_tm,
                              double kappa, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_2_recon12_dirac(double *v2, double *up, double *v1,
                                      double *xp,
                                      double *buf_xp, double *buf_xm,
                                      double *buf_yp, double *buf_ym,
                                      double *buf_zp, double *buf_zm,
                                      double *buf_tp, double *buf_tm,
                                      double kappa, int *bc,
                                      int *Nsize, int *do_comm, int *Leo,
                                      const int ieo, const int iflag);

  void mult_wilson_gm5_dirac(double *v2, double *v1, int *Nsize);

  // real_t = float

  void mult_wilson_bulk_dirac(float *v2, float *up, float *v1,
                              float kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_wilson_1_dirac(float *buf_xp, float *buf_xm,
                           float *buf_yp, float *buf_ym,
                           float *buf_zp, float *buf_zm,
                           float *buf_tp, float *buf_tm,
                           float *up, float *v1,
                           int *bc, int *Nsize, int *do_comm);

  void mult_wilson_2_dirac(float *v2, float *up, float *v1,
                           float *buf_xp, float *buf_xm,
                           float *buf_yp, float *buf_ym,
                           float *buf_zp, float *buf_zm,
                           float *buf_tp, float *buf_tm,
                           float kappa, int *bc,
                           int *Nsize, int *do_comm);

  void mult_wilson_eo_bulk_dirac(float *v2, float *up, float *v1,
                                 float *xp,
                                 float kappa, int *bc,
                                 int *Nsize, int *do_comm, int *Leo,
                                 const int ieo, const int iflag);

  void mult_wilson_eo_bulk_recon12_dirac(float *v2, float *up, float *v1,
                                         float *xp,
                                         float kappa, int *bc,
                                         int *Nsize, int *do_comm, int *Leo,
                                         const int ieo, const int iflag);

  void mult_wilson_eo_1_dirac(float *buf_xp, float *buf_xm,
                              float *buf_yp, float *buf_ym,
                              float *buf_zp, float *buf_zm,
                              float *buf_tp, float *buf_tm,
                              float *up, float *v1, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_1_recon12_dirac(float *buf_xp, float *buf_xm,
                                      float *buf_yp, float *buf_ym,
                                      float *buf_zp, float *buf_zm,
                                      float *buf_tp, float *buf_tm,
                                      float *up, float *v1, int *bc,
                                      int *Nsize, int *do_comm, int *Leo,
                                      const int ieo, const int iflag);

  void mult_wilson_eo_2_dirac(float *v2, float *up, float *v1,
                              float *xp,
                              float *buf_xp, float *buf_xm,
                              float *buf_yp, float *buf_ym,
                              float *buf_zp, float *buf_zm,
                              float *buf_tp, float *buf_tm,
                              float kappa, int *bc,
                              int *Nsize, int *do_comm, int *Leo,
                              const int ieo, const int iflag);

  void mult_wilson_eo_2_recon12_dirac(float *v2, float *up, float *v1,
                                      float *xp,
                                      float *buf_xp, float *buf_xm,
                                      float *buf_yp, float *buf_ym,
                                      float *buf_zp, float *buf_zm,
                                      float *buf_tp, float *buf_tm,
                                      float kappa, int *bc,
                                      int *Nsize, int *do_comm, int *Leo,
                                      const int ieo, const int iflag);

  void mult_wilson_gm5_dirac(float *v2, float *v1, int *Nsize);
}
#endif
