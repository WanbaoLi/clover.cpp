/*!
      @file    bridgeQXS_Clover.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef BRIDGEQXS_CLOVER_INCLUDED
#define BRIDGEQXS_CLOVER_INCLUDED

namespace BridgeQXS {
#ifdef CHIRAL_ROTATION
  void mult_clover_eo_bulk_cswinv_dirac_chrot(
    double *v2, double *v1_tmp, double *up, double *ct, double *v1, double *xp,
    double kappa, int *bc, int *Nsize, int *do_comm, int *Leo, const int ieo,
    const int iflag_aypx);

  void mult_clover_eo_bulk_cswinv_recon12_dirac_chrot(
    double *v2, double *v1_tmp, double *up, double *ct, double *v1, double *xp,
    double kappa, int *bc, int *Nsize, int *do_comm, int *Leo, const int ieo,
    const int iflag_aypx);

  void mult_clover_eo_2_cswinv_dirac_chrot(
    double *v2, double *v1_tmp, double *up, double *ct, double *v1, double *xp,
    double *buf_xp, double *buf_xm, double *buf_yp, double *buf_ym,
    double *buf_zp, double *buf_zm, double *buf_tp, double *buf_tm,
    double kappa, int *bc, int *Nsize, int *do_comm, int *Leo, const int ieo,
    const int iflag_aypx);

  void mult_clover_eo_2_cswinv_recon12_dirac_chrot(
    double *v2, double *v1_tmp, double *up, double *ct, double *v1, double *xp,
    double *buf_xp, double *buf_xm, double *buf_yp, double *buf_ym,
    double *buf_zp, double *buf_zm, double *buf_tp, double *buf_tm,
    double kappa, int *bc, int *Nsize, int *do_comm, int *Leo, const int ieo,
    const int iflag_aypx);

  void mult_clover_eo_bulk_cswinv_dirac_chrot(
    float *v2, float *v1_tmp, float *up, float *ct, float *v1, float *xp,
    float kappa, int *bc, int *Nsize, int *do_comm, int *Leo, const int ieo,
    const int iflag_aypx);

  void mult_clover_eo_bulk_cswinv_recon12_dirac_chrot(
    float *v2, float *v1_tmp, float *up, float *ct, float *v1, float *xp,
    float kappa, int *bc, int *Nsize, int *do_comm, int *Leo, const int ieo,
    const int iflag_aypx);

  void mult_clover_eo_2_cswinv_dirac_chrot(
    float *v2, float *v1_tmp, float *up, float *ct, float *v1, float *xp,
    float *buf_xp, float *buf_xm, float *buf_yp, float *buf_ym,
    float *buf_zp, float *buf_zm, float *buf_tp, float *buf_tm,
    float kappa, int *bc, int *Nsize, int *do_comm, int *Leo, const int ieo,
    const int iflag_aypx);

  void mult_clover_eo_2_cswinv_recon12_dirac_chrot(
    float *v2, float *v1_tmp, float *up, float *ct, float *v1, float *xp,
    float *buf_xp, float *buf_xm, float *buf_yp, float *buf_ym,
    float *buf_zp, float *buf_zm, float *buf_tp, float *buf_tm,
    float kappa, int *bc, int *Nsize, int *do_comm, int *Leo, const int ieo,
    const int iflag_aypx);

#endif

  void mult_clover_bulk_dirac(double *v2, double *up, double *ct,
                              double *v1, double kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_clover_bulk_dirac_chrot(
    double *v2, double *up, double *ct,
    double *v1, double kappa, int *bc,
    int *Nsize, int *do_comm);

  void mult_clover_dd_dirac(double *v2, double *up, double *ct,
                            double *v1, double kappa, int *bc,
                            int *Nsize, int *block_size, int ieo);

  void mult_clover_dd_dirac_chrot(
    double *v2, double *up, double *ct,
    double *v1, double kappa, int *bc,
    int *Nsize, int *block_size, int ieo);

  void mult_clover_cswinv_dirac_chrot(
    double *v2, double *ct, double *v1,
    int *Nsize);


  void mult_clover_bulk_dirac(float *v2, float *up, float *ct,
                              float *v1, float kappa, int *bc,
                              int *Nsize, int *do_comm);

  void mult_clover_bulk_dirac_chrot(
    float *v2, float *up, float *ct,
    float *v1, float kappa, int *bc,
    int *Nsize, int *do_comm);

  void mult_clover_dd_dirac(float *v2, float *up, float *ct,
                            float *v1, float kappa, int *bc,
                            int *Nsize, int *block_size, int ieo);

  void mult_clover_dd_dirac_chrot(
    float *v2, float *up, float *ct,
    float *v1, float kappa, int *bc,
    int *Nsize, int *block_size, int ieo);

  void mult_clover_cswinv_dirac_chrot(
    float *v2, float *ct, float *v1,
    int *Nsize);
}
#endif
