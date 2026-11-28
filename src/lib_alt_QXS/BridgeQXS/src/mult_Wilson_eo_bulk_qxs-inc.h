// Shared EO stencil. The caller consumes each site before its accumulator expires.
#ifndef MULT_WILSON_EO_BULK_QXS_INCLUDED
#define MULT_WILSON_EO_BULK_QXS_INCLUDED

namespace {
inline bool eo_halo_site(int x, int y, int z, int t, const int *n, const int *comm)
{
  return (comm[0] && (x == 0 || x == n[0] - 1))
      || (comm[1] && (y == 0 || y == n[1] - 1))
      || (comm[2] && (z == 0 || z == n[2] - 1))
      || (comm[3] && (t == 0 || t == n[3] - 1));
}

// A temporary owns the expanded SIMD tile until the stencil call returns.
template<bool Recon> struct EoLink;
template<> struct EoLink<false> {
  real_t *get(real_t *u, int site, int block, int volume,
              const int *, const int *, const int *) {
    return u + NDF * VLEN * (site + block * volume);
  }
};

template<> struct EoLink<true> {
  alignas(64) real_t data[NDF * VLEN];

  real_t *get(real_t *u, int site, int block, int volume,
              const int *n, const int *bc, const int *leo) {
    const int dir = block / 2;
    const real_t *src = u + BridgeQXS::NDF12 * VLEN * (site + block * volume);
    for (int c = 0; c < 3; ++c) {
      // QXS component order is column, row, real/imaginary, SIMD lane.
      real_t *out = data + 6 * VLEN * c;
      for (int i = 0; i < 4 * VLEN; ++i) out[i] = src[4 * VLEN * c + i];
      const real_t *a = src + 4 * VLEN * ((c + 1) % 3);
      const real_t *b = src + 4 * VLEN * ((c + 2) % 3);
      for (int k = 0; k < VLEN; ++k) {
        // Column c of conj(row0 x row1).
        out[4 * VLEN + k] =
          ((a[k] * b[2 * VLEN + k] - a[VLEN + k] * b[3 * VLEN + k])
         - (b[k] * a[2 * VLEN + k] - b[VLEN + k] * a[3 * VLEN + k]));
        out[5 * VLEN + k] =
          ((b[k] * a[3 * VLEN + k] + b[VLEN + k] * a[2 * VLEN + k])
         - (a[k] * b[3 * VLEN + k] + a[VLEN + k] * b[2 * VLEN + k]));
      }
    }
    // The cross product contains the boundary phase twice. Correct row 2
    // at the origin of this link, irrespective of the hopping direction.
    if (bc[dir] != 1 && bc[dir] != 0) {
      const int yzt = site / n[0];
      const int tile[] = {site % n[0], yzt % n[1], yzt / n[1] % n[2],
                          yzt / (n[1] * n[2])};
      if (tile[dir] == n[dir] - 1) {
        const real_t scale = real_t(1) / bc[dir];
        for (int k = 0; k < VLEN; ++k) {
          if (dir == 0 && (k % VLENX != VLENX - 1 ||
              (block + leo[VLENY * yzt + k / VLENX]) % 2 != 1)) continue;
          if (dir == 1 && k / VLENX != VLENY - 1) continue;
          for (int c = 0; c < 3; ++c) {
            data[(6 * c + 4) * VLEN + k] *= scale;
            data[(6 * c + 5) * VLEN + k] *= scale;
          }
        }
      }
    }
    return data;
  }
};

template<bool Recon> struct EoGauge {
  real_t *links;
  int volume;
  const int *size, *boundary, *leo;

  // Each argument to a shifted load gets its own scratch tile. The
  // temporary lasts through the enclosing call, including two-tile x/y loads.
  real_t *link(int block, int site, EoLink<Recon>&& scratch = {}) const {
    return scratch.get(links, site, block, volume, size, boundary, leo);
  }
};

template<bool Recon, class Store>
inline void wilson_eo_bulk(real_t *up, real_t *v1, const real_t *ct,
                           int *Nsize, int *do_comm, int *Leo, int ieo,
                           const int *bc, const Store& store)
{
  int Nx2v = Nsize[0];
  int Ny   = Nsize[1];
  int Nz   = Nsize[2];
  int Nt   = Nsize[3];

  int Nst2v = Nx2v * Ny * Nz * Nt;
  const EoGauge<Recon> gauge = {up, Nst2v, Nsize, bc, Leo};

  int Nxy2  = Nx2v * Ny;
  int Nxyz2 = Nx2v * Ny * Nz;

  svbool_t pg1e_xp, pg1e_xm;
  svbool_t pg1o_xp, pg1o_xm;
  svbool_t pg1_yp, pg1_ym;

  svuint_t idx1e_xp, idx1o_xp, idx1e_xm, idx1o_xm;
  svuint_t idx1_yp, idx1_ym;
  set_idx_predicate_xp_eo(pg1e_xp, idx1e_xp, 0);
  set_idx_predicate_xp_eo(pg1o_xp, idx1o_xp, 1);
  set_idx_predicate_xm_eo(pg1e_xm, idx1e_xm, 0);
  set_idx_predicate_xm_eo(pg1o_xm, idx1o_xm, 1);
  set_idx_predicate_yp(pg1_yp, idx1_yp);
  set_idx_predicate_ym(pg1_ym, idx1_ym);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nst2v);

  for (int site = is; site < ns; ++site) {
    if (ct && site + 1 < Nst2v) {
      __builtin_prefetch(v1 + VLEN * NVCD * (site + 1));
      __builtin_prefetch(ct + VLEN * NDF * ND * (site + 1));
    }
    int ix   = site % Nx2v;
    int iyzt = site / Nx2v;
    int iy   = iyzt % Ny;
    int izt  = site / Nxy2;
    int iz   = izt % Nz;
    int it   = izt / Nz;
    int ixy  = ix + Nx2v * iy;
    int ixyz = ixy + Nxy2 * iz;
    int jeo  = (ieo + Leo[VLENY * iyzt]) % 2;

    Vsimd_t v2v[NVCD];
    clear_vec(v2v, NVCD);

    if ((ix < Nx2v - 1) || (do_comm[0] == 0)) {
      int nei = ix + 1 + Nx2v * iyzt;
      if (ix == Nx2v - 1) nei = 0 + Nx2v * iyzt;
      const int block = ieo + 2 * 0;

      if (jeo == 0) {
        mult_wilson_eo_xpb(pg1e_xp, idx1e_xp,
                           v2v, gauge.link(block, site),
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {
        mult_wilson_eo_xpb(pg1o_xp, idx1o_xp,
                           v2v, gauge.link(block, site),
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }
    }

    if ((ix > 0) || (do_comm[0] == 0)) {
      int nei = ix - 1 + Nx2v * iyzt;
      if (ix == 0) nei = Nx2v - 1 + Nx2v * iyzt;
      const int block = 1 - ieo + 2 * 0;
      if (jeo == 0) {
        mult_wilson_eo_xmb(pg1e_xm, idx1e_xm, v2v,
                           gauge.link(block, site), gauge.link(block, nei),
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      } else {
        mult_wilson_eo_xmb(pg1o_xm, idx1o_xm, v2v,
                           gauge.link(block, site), gauge.link(block, nei),
                           &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
      }
    }

    if ((iy < Ny - 1) || (do_comm[1] == 0)) {
      int    iy2 = (iy + 1) % Ny;
      int    nei = ix + Nx2v * (iy2 + Ny * izt);
      const int block = ieo + 2 * 1;
      mult_wilson_ypb(pg1_yp, idx1_yp, v2v,
                      gauge.link(block, site),
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if ((iy > 0) || (do_comm[1] == 0)) {
      int    iy2 = (iy - 1 + Ny) % Ny;
      int    nei = ix + Nx2v * (iy2 + Ny * izt);
      const int block = 1 - ieo + 2 * 1;
      mult_wilson_ymb(pg1_ym, idx1_ym, v2v,
                      gauge.link(block, site), gauge.link(block, nei),
                      &v1[VLEN * NVCD * site], &v1[VLEN * NVCD * nei]);
    }

    if ((iz < Nz - 1) || (do_comm[2] == 0)) {
      int    iz2 = (iz + 1) % Nz;
      int    nei = ixy + Nxy2 * (iz2 + Nz * it);
      const int block = ieo + 2 * 2;
      mult_wilson_zpb(v2v, gauge.link(block, site), &v1[VLEN * NVCD * nei]);
    }

    if ((iz > 0) || (do_comm[2] == 0)) {
      int    iz2 = (iz - 1 + Nz) % Nz;
      int    nei = ixy + Nxy2 * (iz2 + Nz * it);
      const int block = 1 - ieo + 2 * 2;
      mult_wilson_zmb(v2v, gauge.link(block, nei), &v1[VLEN * NVCD * nei]);
    }

    if ((it < Nt - 1) || (do_comm[3] == 0)) {
      int    it2 = (it + 1) % Nt;
      int    nei = ixyz + Nxyz2 * it2;
      const int block = ieo + 2 * 3;
      mult_wilson_tpb_dirac(v2v, gauge.link(block, site), &v1[VLEN * NVCD * nei]);
    }

    if ((it > 0) || (do_comm[3] == 0)) {
      int    it2 = (it - 1 + Nt) % Nt;
      int    nei = ixyz + Nxyz2 * it2;
      const int block = 1 - ieo + 2 * 3;
      mult_wilson_tmb_dirac(v2v, gauge.link(block, nei), &v1[VLEN * NVCD * nei]);
    }

    store(site, ix, iy, iz, it, v2v);
  }
}
} // namespace
#endif
