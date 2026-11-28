# clover.cpp

Bridge++ 2.1.3 with the Wilson–Clover propagator driver from clover.cpp develop.
Build on ITP node5001 with Open MPI 4.1.4, GCC, CMake and Ninja.

```bash
./build.sh
```

## cfg3000 point-source propagator

From the repository root on node5001, run this command for source column 0
(spin 0, color 0). Use a new output directory for each run.

```bash
export PATH="/public2/software/mpi/openmpi/4.1.4/gnu/bin:$PATH"
export OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=false
cfg="$HOME/Configure_clover/beta6.20_mu-0.2790_ms-0.2400_L24x64/beta6.20_mu-0.2790_ms-0.2400_L24x64_cfg_3000.lime"

mpirun -np 16 --map-by ppr:2:numa:PE=16 --bind-to core --mca pml ob1 --mca btl self,vader \
  build/clover-solve --method mg --lattice 24,24,24,64 --grid 1,1,2,8 \
  --threads 16 --gauge "$cfg" --source point --position 0,0,0,0 --rhs 0 \
  --mass -0.279 --csw 1.160678675188 --rho 0.125 \
  --precision double --tol 1e-9 --accept 1e-9 --fields full --output build/prop-cfg3000
```

Use `--rhs all` to compute all 12 columns. Check `columns.tsv` for
`independent_relative` and `accepted`. A successful output includes
`COMPLETE` and the per-rank FP64 field files.

## MG timing

### C24P34 cfg20000 wall-source benchmark

Run `tests/benchmark_mg_c24p34.sh` on node5001. The gauge is stored in
`$HOME/Configure_clover/beta6.20_mu-0.2770_ms-0.2310_L24x64/` and was copied
from the Kunshan ensemble. SHA-256 is
`47969c7fe9b4a4716d5eab40e5a1dbdaef09abdb055934000e034a99d20d2d30`.

The light-quark run follows the C24P34 entry in the supplied PyQUDA input:
24³×64, mass -0.277, Clover coefficient 1/0.951479³, one stout step at
rho 0.125, antiperiodic time, t=0 wall source, at most 1000 outer iterations. The script uses SAP blocks 6×6×6×4,
24 Bridge++ setup vectors, four setup refinement steps and four smoothing
sweeps. Bridge++ uses CPU FBiCGStab and solves the twelve columns in order
with a shared hierarchy. Results are under `build/mg-c24p34.*/`.

The benchmark reports one and twelve columns for both precisions. Pass single
or double to run one precision. The FP32 and FP64 true residual limits are 1e-6
and 1e-9. The independent check uses FP64 in both cases. For FP32 the recursive
stop threshold is 2.5e-7 to allow for rounding error.

## Wilson–Clover 2+1 ensemble generation

The new build-tree launcher is build/bin/clover-hmc. It needs Python 3.9+
with NumPy and uses the MPI launcher available on PATH. Set CLOVER_HMC_PYTHON
to select an interpreter. To run one small test trajectory:

    build/bin/clover-hmc run --dims 2 2 2 2 --seed 1703 --tau 0.12 \
      --outer-steps 1 --inner-steps 1 --output build/hmc-smoke

The resume and export commands share the same checkpoint format. See
build/bin/clover-hmc run --help for the archived Chroma restart option.
The original Wilson Nf=2 generator remains available as build/bridge-hmc.

The default QXS build uses the generic implementation. On AArch64,
CMake accepts -DCLOVER_QXS=sve512 with a supporting SVE2 toolchain.

## MG outer solver

The default is flexible BiCGStab. Select restarted flexible GCR with
`--method mg --mg-outer gcr`. The default restart length is 16. The fine
operator, MG setup, coarse solver and SAP smoother are shared by both modes.
`outer_loops` counts the GCR steps, while `nconv` includes fine-operator
applications used to refresh the GCR residual. The independently checked
full-field residual remains the acceptance criterion.

Run `tests/benchmark_mg_c24p34.sh all gcr` for the one- and twelve-column
measurements. Omit the second argument to benchmark the default solver.
