# Old multigrid.cpp ports · node5001

The reference is cb57afc. These runs compare GCR before and after the
QXS residual-and-norm fusion, on cfg20000 with the same 16 MPI × 16 OpenMP
layout, the same wall source, and the established FP32/FP64 residual gates.
One run was measured for each case. GCR is now the default for the
clover-solve MG path. FBiCGStab remains selectable.

| Precision | Columns | Outer iterations | Max true relative residual | MG setup (s) | Solve (s) |
|---|---:|---:|---:|---:|---:|
| single | 1 | 31 → 31 | 3.256009e-07 → 3.256009e-07 | 25.782 → 25.870 | 3.854 → 3.833 (-0.5%) |
| single | 12 | 382 → 382 | 4.369791e-07 → 4.369791e-07 | 26.099 → 25.992 | 48.031 → 48.465 (+0.9%) |
| double | 1 | 48 → 48 | 8.776504e-10 → 8.776504e-10 | 25.872 → 26.274 | 6.367 → 6.367 (+0.0%) |
| double | 12 | 588 → 588 | 9.837571e-10 → 9.837571e-10 | 26.106 → 26.194 | 77.020 → 77.173 (+0.2%) |

The old code at 797ac9c supplied the single-pass reliable-residual pattern.
HMC ports reuse Schur scratch (fc5c530), apply the Schur normal operator in
one OpenMP team (91265d3), and combine shifted-CG curvature norms in a real
MPI reduction (97eb38c). RHMC also reuses the per-pole verification fields.

Two one-trajectory runs per version were measured on 4⁴, using one MPI
rank and four OpenMP threads. Gauge fields, RNG state, acceptance and
ΔH were identical in all four runs. Median trajectory time before checkpoint
I/O was 3.038 → 3.044 s. This is a small-grid smoke timing.

The old fine/coarse halo batching and setup-vector batching use a native
layout and batched MG cycle. They need Bridge++ batch plumbing and their own
accuracy and memory checks. The selected endpoint energy is already cached
by the current HMC controller. The independent pure-gauge force path in the
old code uses the Wilson plaquette action, whereas this HMC uses the
tree-level improved gauge action.
