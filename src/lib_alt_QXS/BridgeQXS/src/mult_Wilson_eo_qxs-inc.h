/*!
      @file    mult_Wilson_eo_qxs-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef MULT_WILSON_EO_QXS_INCLUDED
#define MULT_WILSON_EO_QXS_INCLUDED

#include "mult_common_th-inc.h"

#include "mult_Wilson_eo_bulk_qxs-inc.h"

#define IMPLE 2

//====================================================================
template<bool Recon>
static void mult_wilson_eo_bulk_dirac_impl(real_t *v2, real_t *up,
                                           real_t *v1, real_t *xp,
                                           real_t kappa, int *bc,
                                           int *Nsize, int *do_comm,
                                           int *Leo, const int ieo,
                                           const int iflag)
{
  wilson_eo_bulk<Recon>(up, v1, nullptr, Nsize, do_comm, Leo, ieo, bc,
                 [&](int site, int, int, int, int, Vsimd_t *v2v) {
    svbool_t pg = set_predicate();
    if (iflag == 0) {
      real_t *vv2 = &v2[VLEN * NVCD * site];
      for (int i = 0; i < NVCD; ++i) {
        svreal_t v2t;
        load_vec(pg, v2t, &v2v[i].v[0]);
        scal_vec(pg, v2t, -kappa);
        save_vec(pg, &vv2[VLEN * i], v2t);
      }
    } else {
      mult_wilson_aypx_save(&v2[VLEN * NVCD * site],
                            kappa, v2v, &xp[VLEN * NVCD * site]);
    }
  });
}

void BridgeQXS::mult_wilson_eo_bulk_dirac(real_t *v2, real_t *up,
                                          real_t *v1, real_t *xp,
                                          real_t kappa, int *bc,
                                          int *Nsize, int *do_comm,
                                          int *Leo, const int ieo,
                                          const int iflag)
{
  mult_wilson_eo_bulk_dirac_impl<false>(v2, up, v1, xp, kappa, bc, Nsize, do_comm, Leo, ieo,
      iflag);
}

void BridgeQXS::mult_wilson_eo_bulk_recon12_dirac(real_t *v2, real_t *up,
                                                  real_t *v1, real_t *xp,
                                                  real_t kappa, int *bc,
                                                  int *Nsize, int *do_comm,
                                                  int *Leo, const int ieo,
                                                  const int iflag)
{
  mult_wilson_eo_bulk_dirac_impl<true>(v2, up, v1, xp, kappa, bc, Nsize, do_comm, Leo, ieo,
      iflag);
}


//====================================================================
template<bool Recon>
static void mult_wilson_eo_1_dirac_impl(real_t *buf_xp, real_t *buf_xm,
                                        real_t *buf_yp, real_t *buf_ym,
                                        real_t *buf_zp, real_t *buf_zm,
                                        real_t *buf_tp, real_t *buf_tm,
                                        real_t *up, real_t *v1, int *bc,
                                        int *Nsize, int *do_comm, int *Leo,
                                        const int ieo, const int iflag)
{
  int Nx2v = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nst2v = Nx2v * Ny * Nz * Nt;
  const EoGauge<Recon> gauge = {up, Nst2v, Nsize, bc, Leo};

  int Nxy2  = Nx2v * Ny;
  int Nxyz2 = Nx2v * Ny * Nz;

  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
#if IMPLE == 0
  svint_t svidx_xp, svidx_xm;
#else
  svuint_t svidx_xp, svidx_xm;
#endif
  set_index_xp_eo(svidx_xp);
  set_index_xm_eo(svidx_xm);

  int Nskipx = (VLENY + 1) / 2;


  if (do_comm[0] == 1) {
    int idir = 0;

    int Nyzt = Ny * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nyzt);

    for (int iyzt = is; iyzt < ns; ++iyzt) {
      int jeo = (ieo + Leo[VLENY * iyzt]) % 2;
      {
        int ix     = 0;
        int site   = ix + Nx2v * iyzt;
        int ibf_up = Nskipx * NVC * ND2 * iyzt;
        if (jeo == 0) {
#if VLENY > 1
#if IMPLE == 0
          set_index_xm_eo(svidx_xm);
          mult_wilson_eo_xp1(pg2o_xm, svidx_xm,
                             &buf_xp[ibf_up], &v1[VLEN * NVCD * site]);
#else
          mult_wilson_eo_xp1(pg2o_xm,
                             &buf_xp[ibf_up], &v1[VLEN * NVCD * site]);
#endif
#endif
        } else {
          if (VLENY == 1) ibf_up = Nskipx * NVC * ND2 * (iyzt / 2);
#if IMPLE == 0
          set_index_xm_eo(svidx_xm);
          mult_wilson_eo_xp1(pg2e_xm, svidx_xm,
                             &buf_xp[ibf_up], &v1[VLEN * NVCD * site]);
#else
          mult_wilson_eo_xp1(pg2e_xm,
                             &buf_xp[ibf_up], &v1[VLEN * NVCD * site]);
#endif
        }
      }

      {
        int    ix     = Nx2v - 1;
        int    site   = ix + Nx2v * iyzt;
        int    ibf_dn = Nskipx * NVC * ND2 * iyzt;
        const int block = 1 - ieo + 2 * idir;
        if (jeo == 0) {
          if (VLENY == 1) ibf_dn = Nskipx * NVC * ND2 * (iyzt / 2);
#if IMPLE == 0
          set_index_xp_eo(svidx_xp);
          mult_wilson_eo_xm1(pg2o_xp, svidx_xp, &buf_xm[ibf_dn],
                             gauge.link(block, site), &v1[VLEN * NVCD * site]);
#else
          mult_wilson_eo_xm1(pg2o_xp, &buf_xm[ibf_dn],
                             gauge.link(block, site), &v1[VLEN * NVCD * site]);
#endif
        } else {
#if VLENY > 1
#if IMPLE == 0
          set_index_xp_eo(svidx_xp);
          mult_wilson_eo_xm1(pg2e_xp, svidx_xp, &buf_xm[ibf_dn],
                             gauge.link(block, site), &v1[VLEN * NVCD * site]);
#else
          mult_wilson_eo_xm1(pg2e_xp, &buf_xm[ibf_dn],
                             gauge.link(block, site), &v1[VLEN * NVCD * site]);
#endif
#endif
        }
      }
    }
  }

  if (do_comm[1] > 0) {
    int idir = 1;

    int Nxzt = Nx2v * Nz * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxzt);

    for (int ixzt = is; ixzt < ns; ++ixzt) {
      int ix  = ixzt % Nx2v;
      int izt = ixzt / Nx2v;
      {
        int iy   = 0;
        int site = ix + Nx2v * (iy + Ny * izt);
        int ibf  = VLENX * NVC * ND2 * (ix + Nx2v * izt);
        mult_wilson_yp1(pg2_ym,
                        &buf_yp[ibf], &v1[VLEN * NVCD * site]);
      }
      {
        int    iy   = Ny - 1;
        int    site = ix + Nx2v * (iy + Ny * izt);
        const int block = (1 - ieo) + 2 * idir;
        int    ibf  = VLENX * NVC * ND2 * (ix + Nx2v * izt);
        mult_wilson_ym1(pg2_yp,
                        &buf_ym[ibf], gauge.link(block, site),
                        &v1[VLEN * NVCD * site]);
      }
    }
  }

  if (do_comm[2] > 0) {
    int idir  = 2;
    int Nxyt2 = Nxy2 * Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxyt2);

    for (int ixyt = is; ixyt < ns; ++ixyt) {
      int ixy = ixyt % Nxy2;
      int it  = ixyt / Nxy2;
      {
        int iz   = 0;
        int site = ixy + Nxy2 * (iz + Nz * it);
        int ibf  = VLEN * NVC * ND2 * (ixy + Nxy2 * it);
        mult_wilson_zp1(&buf_zp[ibf], &v1[VLEN * NVCD * site]);
      }
      {
        int    iz   = Nz - 1;
        int    site = ixy + Nxy2 * (iz + Nz * it);
        int    ibf  = VLEN * NVC * ND2 * (ixy + Nxy2 * it);
        const int block = 1 - ieo + 2 * idir;
        mult_wilson_zm1(&buf_zm[ibf], gauge.link(block, site),
                        &v1[VLEN * NVCD * site]);
      }
    }
  }

  if (do_comm[3] > 0) {
    int idir = 3;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Nxyz2);

    for (int ixyz = is; ixyz < ns; ++ixyz) {
      {
        int it   = 0;
        int site = ixyz + Nxyz2 * it;
        mult_wilson_tp1_dirac(&buf_tp[VLEN * NVC * ND2 * ixyz],
                              &v1[VLEN * NVCD * site]);
      }
      {
        int    it   = Nt - 1;
        int    site = ixyz + Nxyz2 * it;
        const int block = 1 - ieo + 2 * idir;
        mult_wilson_tm1_dirac(&buf_tm[VLEN * NVC * ND2 * ixyz],
                              gauge.link(block, site), &v1[VLEN * NVCD * site]);
      }
    }
  }
}

void BridgeQXS::mult_wilson_eo_1_dirac(real_t *buf_xp, real_t *buf_xm,
                                       real_t *buf_yp, real_t *buf_ym,
                                       real_t *buf_zp, real_t *buf_zm,
                                       real_t *buf_tp, real_t *buf_tm,
                                       real_t *up, real_t *v1, int *bc,
                                       int *Nsize, int *do_comm, int *Leo,
                                       const int ieo, const int iflag)
{
  mult_wilson_eo_1_dirac_impl<false>(buf_xp, buf_xm, buf_yp, buf_ym, buf_zp, buf_zm, buf_tp,
      buf_tm, up, v1, bc, Nsize, do_comm, Leo, ieo, iflag);
}

void BridgeQXS::mult_wilson_eo_1_recon12_dirac(real_t *buf_xp, real_t *buf_xm,
                                               real_t *buf_yp, real_t *buf_ym,
                                               real_t *buf_zp, real_t *buf_zm,
                                               real_t *buf_tp, real_t *buf_tm,
                                               real_t *up, real_t *v1, int *bc,
                                               int *Nsize, int *do_comm, int *Leo,
                                               const int ieo, const int iflag)
{
  mult_wilson_eo_1_dirac_impl<true>(buf_xp, buf_xm, buf_yp, buf_ym, buf_zp, buf_zm, buf_tp,
      buf_tm, up, v1, bc, Nsize, do_comm, Leo, ieo, iflag);
}


//====================================================================
template<bool Recon>
static void mult_wilson_eo_2_dirac_impl(real_t *v2, real_t *up, real_t *v1,
                                        real_t *xp,
                                        real_t *buf_xp, real_t *buf_xm,
                                        real_t *buf_yp, real_t *buf_ym,
                                        real_t *buf_zp, real_t *buf_zm,
                                        real_t *buf_tp, real_t *buf_tm,
                                        real_t kappa, int *bc,
                                        int *Nsize, int *do_comm, int *Leo,
                                        const int ieo, const int iflag)
{
  int Nx2v = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nst2v = Nx2v * Ny * Nz * Nt;
  const EoGauge<Recon> gauge = {up, Nst2v, Nsize, bc, Leo};

  int Nxy2  = Nx2v * Ny;
  int Nxyz2 = Nx2v * Ny * Nz;

  svbool_t pg1e_xp, pg2e_xp, pg3e_xp, pg1e_xm, pg2e_xm, pg3e_xm;
  svbool_t pg1o_xp, pg2o_xp, pg3o_xp, pg1o_xm, pg2o_xm, pg3o_xm;
  set_predicate_xp_eo(pg1e_xp, pg2e_xp, pg3e_xp, 0);
  set_predicate_xp_eo(pg1o_xp, pg2o_xp, pg3o_xp, 1);
  set_predicate_xm_eo(pg1e_xm, pg2e_xm, pg3e_xm, 0);
  set_predicate_xm_eo(pg1o_xm, pg2o_xm, pg3o_xm, 1);
  svbool_t pg1_yp, pg2_yp, pg1_ym, pg2_ym;
  set_predicate_yp(pg1_yp, pg2_yp);
  set_predicate_ym(pg1_ym, pg2_ym);
#if IMPLE == 0
  svint_t svidx_xp, svidx_xm;
#else
  svuint_t svidx_xp, svidx_xm;
#endif
  set_index_xp_eo(svidx_xp);
  set_index_xm_eo(svidx_xm);
  svbool_t pg = set_predicate();

  real_t kappa_eo = kappa;
  if (iflag == 0) {
    kappa_eo = -kappa;
  }

  int Nskipx = (VLENY + 1) / 2;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nst2v);

  for (int site = is; site < ns; ++site) {
    int ix   = site % Nx2v;
    int iyzt = site / Nx2v;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy2;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    if (!eo_halo_site(ix, iy, iz, it, Nsize, do_comm)) continue;
    int ixy  = ix + Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int jeo  = (ieo + Leo[VLENY * iyzt]) % 2;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    real_t zL[VLEN * NVCD], uL[VLEN * NDF];

    if ((ix == Nx2v - 1) && (do_comm[0] > 0)) {
      int ibf_up = Nskipx * NVC * ND2 * iyzt;
      if (VLENY == 1) ibf_up = Nskipx * NVC * ND2 * (iyzt / 2);
      const int block = ieo + 2 * 0;
      if (jeo == 0) {
        mult_wilson_eo_xp2(pg1e_xp, pg2e_xp, pg3e_xp, svidx_xp,
                           v2v, gauge.link(block, site),
                           &v1[VLEN * NVCD * site], &buf_xp[ibf_up]);
      } else {
        mult_wilson_eo_xp2(pg1o_xp, pg2o_xp, pg3o_xp, svidx_xp,
                           v2v, gauge.link(block, site),
                           &v1[VLEN * NVCD * site], &buf_xp[ibf_up]);
      }
    }

    if ((ix == 0) && (do_comm[0] > 0)) {
      int ibf_dn = Nskipx * NVC * ND2 * iyzt;
      if (VLENY == 1) ibf_dn = Nskipx * NVC * ND2 * (iyzt / 2);
      const int block = 1 - ieo + 2 * 0;
      if (jeo == 0) {
        mult_wilson_eo_xm2(pg1e_xm, pg2e_xm, pg3e_xm, svidx_xm,
                           v2v, gauge.link(block, site),
                           &v1[VLEN * NVCD * site], &buf_xm[ibf_dn]);
      } else {
        mult_wilson_eo_xm2(pg1o_xm, pg2o_xm, pg3o_xm, svidx_xm,
                           v2v, gauge.link(block, site),
                           &v1[VLEN * NVCD * site], &buf_xm[ibf_dn]);
      }
    }

    if ((iy == Ny - 1) && (do_comm[1] > 0)) {
      int    ibf = VLENX * NVC * ND2 * (ix + Nx2v * izt);
      const int block = ieo + 2 * 1;
      mult_wilson_yp2(pg1_yp, pg2_yp,
                      v2v, gauge.link(block, site),
                      &v1[VLEN * NVCD * site], &buf_yp[ibf]);
    }

    if ((iy == 0) && (do_comm[1] > 0)) {
      int    ibf = VLENX * NVC * ND2 * (ix + Nx2v * izt);
      const int block = 1 - ieo + 2 * 1;
      mult_wilson_ym2(pg1_ym, pg2_ym,
                      v2v, gauge.link(block, site),
                      &v1[VLEN * NVCD * site], &buf_ym[ibf]);
    }

    if ((iz == Nz - 1) && (do_comm[2] > 0)) {
      int    ibf = VLEN * NVC * ND2 * (ixy + Nxy2 * it);
      const int block = ieo + 2 * 2;
      mult_wilson_zp2(v2v, gauge.link(block, site), &buf_zp[ibf]);
    }

    if ((iz == 0) && (do_comm[2] > 0)) {
      int ibf = VLEN * NVC * ND2 * (ixy + Nxy2 * it);
      mult_wilson_zm2(v2v, &buf_zm[ibf]);
    }

    if ((it == Nt - 1) && (do_comm[3] > 0)) {
      const int block = ieo + 2 * 3;
      mult_wilson_tp2_dirac(v2v, gauge.link(block, site),
                            &buf_tp[VLEN * NVC * ND2 * ixyz]);
    }

    if ((it == 0) && (do_comm[3] > 0)) {
      mult_wilson_tm2_dirac(v2v, &buf_tm[VLEN * NVC * ND2 * ixyz]);
    }


    real_t *ww2 = &v2v[0].v[0];
    real_t *vv2 = &v2[VLEN * NVCD * site];
    for (int i = 0; i < NVCD; i += 8) {
      svreal_t wt1, wt2, wt3, wt4;
      svreal_t vt1, vt2, vt3, vt4;
      svreal_t wt5, wt6, wt7, wt8;
      svreal_t vt5, vt6, vt7, vt8;
      load_vec(pg, vt1, &vv2[VLEN * (i)]);
      load_vec(pg, vt2, &vv2[VLEN * (i + 1)]);
      load_vec(pg, vt3, &vv2[VLEN * (i + 2)]);
      load_vec(pg, vt4, &vv2[VLEN * (i + 3)]);

      load_vec(pg, wt1, &ww2[VLEN * (i)]);
      load_vec(pg, wt2, &ww2[VLEN * (i + 1)]);
      load_vec(pg, wt3, &ww2[VLEN * (i + 2)]);
      load_vec(pg, wt4, &ww2[VLEN * (i + 3)]);

      axpy_vec(pg, vt1, kappa_eo, wt1);
      axpy_vec(pg, vt2, kappa_eo, wt2);
      axpy_vec(pg, vt3, kappa_eo, wt3);

      load_vec(pg, vt5, &vv2[VLEN * (i + 4)]);
      load_vec(pg, vt6, &vv2[VLEN * (i + 5)]);
      load_vec(pg, vt7, &vv2[VLEN * (i + 6)]);
      load_vec(pg, vt8, &vv2[VLEN * (i + 7)]);

      load_vec(pg, wt5, &ww2[VLEN * (i + 4)]);
      load_vec(pg, wt6, &ww2[VLEN * (i + 5)]);
      load_vec(pg, wt7, &ww2[VLEN * (i + 6)]);
      load_vec(pg, wt8, &ww2[VLEN * (i + 7)]);

      axpy_vec(pg, vt4, kappa_eo, wt4);
      axpy_vec(pg, vt5, kappa_eo, wt5);
      axpy_vec(pg, vt6, kappa_eo, wt6);
      axpy_vec(pg, vt7, kappa_eo, wt7);
      axpy_vec(pg, vt8, kappa_eo, wt8);

      save_vec(pg, &vv2[VLEN * (i)], vt1);
      save_vec(pg, &vv2[VLEN * (i + 1)], vt2);
      save_vec(pg, &vv2[VLEN * (i + 2)], vt3);
      save_vec(pg, &vv2[VLEN * (i + 3)], vt4);
      save_vec(pg, &vv2[VLEN * (i + 4)], vt5);
      save_vec(pg, &vv2[VLEN * (i + 5)], vt6);
      save_vec(pg, &vv2[VLEN * (i + 6)], vt7);
      save_vec(pg, &vv2[VLEN * (i + 7)], vt8);
    }
  }
}

void BridgeQXS::mult_wilson_eo_2_dirac(real_t *v2, real_t *up, real_t *v1,
                                       real_t *xp,
                                       real_t *buf_xp, real_t *buf_xm,
                                       real_t *buf_yp, real_t *buf_ym,
                                       real_t *buf_zp, real_t *buf_zm,
                                       real_t *buf_tp, real_t *buf_tm,
                                       real_t kappa, int *bc,
                                       int *Nsize, int *do_comm, int *Leo,
                                       const int ieo, const int iflag)
{
  mult_wilson_eo_2_dirac_impl<false>(v2, up, v1, xp, buf_xp, buf_xm, buf_yp, buf_ym, buf_zp,
      buf_zm, buf_tp, buf_tm, kappa, bc, Nsize, do_comm, Leo, ieo, iflag);
}

void BridgeQXS::mult_wilson_eo_2_recon12_dirac(real_t *v2, real_t *up, real_t *v1,
                                               real_t *xp,
                                               real_t *buf_xp, real_t *buf_xm,
                                               real_t *buf_yp, real_t *buf_ym,
                                               real_t *buf_zp, real_t *buf_zm,
                                               real_t *buf_tp, real_t *buf_tm,
                                               real_t kappa, int *bc,
                                               int *Nsize, int *do_comm, int *Leo,
                                               const int ieo, const int iflag)
{
  mult_wilson_eo_2_dirac_impl<true>(v2, up, v1, xp, buf_xp, buf_xm, buf_yp, buf_ym, buf_zp,
      buf_zm, buf_tp, buf_tm, kappa, bc, Nsize, do_comm, Leo, ieo, iflag);
}



#endif
//============================================================END=====
