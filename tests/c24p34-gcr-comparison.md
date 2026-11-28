# C24P34 cfg20000 · MG outer-solver comparison

Measured on node5001 with 16 MPI ranks and 16 OpenMP threads per rank.
Both use the 24³×64 C24P34 field, a t=0 wall source, a 6×6×6×4 SAP block,
24 setup vectors, four setup refinements and four smoothing sweeps.
Each algorithm and precision was run once for RHS 0 and once for all
twelve columns. Twelve columns are solved sequentially with a shared
hierarchy. Times exclude file output and independent FP64 verification.
The solver tolerance is 1e-6 in FP32 and 1e-9 in FP64.

GCR uses the paired direction/image basis and true-residual refresh from
the old multigrid.cpp native solver. The 16-direction restart follows its
bounded-workspace configuration. The implementation reuses the QXS fused
dot product and norm reduction. Both algorithms use the same initial
guess x=b and the same MG cycle. For FP32 GCR, the recursive cutoff is
2.5e-7, with a refreshed FP32 early-stop target of 8e-7. The final independent FP64 gate is 1e-6 for both solvers.

| Outer precision | Columns | Iterations, FBiCGStab → GCR | Max true relative residual | MG setup (s) | Solve (s) |
|---|---:|---:|---:|---:|---:|
| single | 1 | 23 → 31 | 4.586913e-07 → 3.256009e-07 | 26.572 → 26.360 | 5.954 → 3.879 (-34.9%) |
| single | 12 | 314 → 382 | 6.351899e-07 → 4.369791e-07 | 26.709 → 25.911 | 80.667 → 48.389 (-40.0%) |
| double | 1 | 38 → 48 | 7.533936e-10 → 8.776504e-10 | 26.639 → 26.318 | 10.399 → 6.397 (-38.5%) |
| double | 12 | 516 → 588 | 9.827129e-10 → 9.837571e-10 | 26.729 → 26.056 | 136.575 → 77.166 (-43.5%) |

Raw run directories: mg-c24p34.sMLDom (FBiCGStab) and mg-c24p34-gcr.RCRxAF (GCR).
The companion JSON contains per-column data and phase timings.

Old repository reviewed at all 15 branch tips. The scalar reference is
\`src/core/gcr.cpp\`, the production solver is \`src/solver/native/solver.cpp\`,
and the independent batched update is in \`src/solver/native/batched_solver.cpp\`.
Coarse-halo batching (71160fa), SVE neighbor-projection reuse (ea30336),
and adaptive coarse-iteration escalation (5b99879) are separate candidates.
This comparison uses the current Bridge++ MG cycle for both outer solvers.
