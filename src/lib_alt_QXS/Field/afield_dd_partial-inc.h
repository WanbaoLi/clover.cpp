// Per-site block algebra for SAP blocks that split a QXS register.
// The full-pack kernels remain available when both block lengths are aligned.
template<typename INDEX>
bool partial_xy_block(const INDEX& index)
{
  return index.block_size(0) % VLENX || index.block_size(1) % VLENY;
}

template<typename INDEX, typename ACTION>
void for_partial_blocks(const INDEX& index, const ACTION& action)
{
  int thread, threads, first, end;
  set_threadtask(thread, threads, first, end, index.coarse_nvol());
  for (int block = first; block < end; ++block) action(block);
}

// The offset of an imaginary part is one vector beyond its real part.
template<typename INDEX, typename AFIELD, typename ACTION>
void for_block_components(const INDEX& index, const AFIELD& field,
                          int block, const ACTION& action)
{
  int bx, by, bz, bt;
  index.get_coarse_coord(bx, by, bz, bt, block);
  const int sx = index.block_size(0), sy = index.block_size(1);
  const int sz = index.block_size(2), st = index.block_size(3);
  const int nxv = index.fine_lattice_size(0) / VLENX;
  const int nyv = index.fine_lattice_size(1) / VLENY;
  const int nz = index.fine_lattice_size(2);
  const int packs = field.nvol() / VLEN;
  for (int t = bt*st; t < (bt+1)*st; ++t)
    for (int z = bz*sz; z < (bz+1)*sz; ++z)
      for (int y = by*sy; y < (by+1)*sy; ++y)
        for (int x = bx*sx; x < (bx+1)*sx; ++x) {
          int pack = x/VLENX + nxv*(y/VLENY + nyv*(z + nz*t));
          int lane = x%VLENX + VLENX*(y%VLENY);
          for (int ex = 0; ex < field.nex(); ++ex)
            for (int k = 0; k < field.nin(); k += 2) {
              int i = lane + VLEN*(k + field.nin()*(pack + packs*ex));
              action(i, i + VLEN);
            }
        }
}

template<typename INDEX, typename AFIELD, typename COEFF>
void update_partial_blocks(AFIELD& v, const COEFF *a, const AFIELD& w,
                           int parity, typename AFIELD::real_t fac,
                           const INDEX& index, bool add)
{
  using Real = typename AFIELD::real_t;
  using Complex = typename AFIELD::complex_t;
  Real *vp = v.ptr(0);
  const Real *wp = const_cast<AFIELD&>(w).ptr(0);
  for_partial_blocks(index, [&](int block) {
    if (index.block_eo(block) == 1-parity) return;
    const Complex alpha = fac * Complex(a[block]);
    for_block_components(index, v, block, [&](int ir, int ii) {
      const Complex value = alpha * Complex(wp[ir], wp[ii]);
      if (add) {vp[ir] += real(value); vp[ii] += imag(value);}
      else {vp[ir] = real(value); vp[ii] = imag(value);}
    });
  });
}
