#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."

export PATH="/public2/software/mpi/openmpi/4.1.4/gnu/bin:$PATH"
export OMP_NUM_THREADS=16 OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=false
cfg="$HOME/Configure_clover/beta6.20_mu-0.2770_ms-0.2310_L24x64/beta6.20_mu-0.2770_ms-0.2310_L24x64_cfg_20000.lime"
test -x build/clover-solve || { echo 'Run ./build.sh first.' >&2; exit 1; }
test -r "$cfg" || { echo "Missing C24P34 cfg20000: $cfg" >&2; exit 1; }

case ${1:-all} in
  single|double) precisions=("$1") ;;
  all) precisions=(single double) ;;
  *) echo "Usage: $0 [single|double|all]" >&2; exit 2 ;;
esac
out=$(mktemp -d build/mg-c24p34.XXXXXX)
printf 'precision\tcolumns\touter_iters\tmax_true_residual\tmg_setup_s\tsolve_s\n' | tee "$out/timing.tsv"

for precision in "${precisions[@]}"; do
  if [[ $precision == single ]]; then tol=1e-6; else tol=1e-9; fi
  for rhs in 0 all; do
    result="$out/$precision-rhs-$rhs"
    expected=1
    if [[ $rhs == all ]]; then expected=12; fi
    if ! mpirun -np 16 --map-by ppr:2:numa:PE=16 --bind-to core \
      --mca pml ob1 --mca btl self,vader \
      build/clover-solve --method mg --lattice 24,24,24,64 --grid 1,1,2,8 \
      --threads "$OMP_NUM_THREADS" --gauge "$cfg" --source wall --position 0,0,0,0 \
      --rhs "$rhs" --mass -0.277 --csw 1.16092022578464 --rho 0.125 \
      --sap-block 6,6,6,4 --mg-max-iters 1000 --profile-mg --nvec 24 --setup-steps 4 --smooth-iters 4 \
      --precision "$precision" --tol "$tol" --accept "$tol" --fields none --output "$result" > "$result.log" 2>&1 ||
      [[ ! -f "$result/COMPLETE" ]]; then
      echo "Solve failed. See $result.log" >&2
      exit 1
    fi

    awk -F '\t' -v precision="$precision" -v expected="$expected" '
      NR == FNR { if ($1 == "setup_seconds") { setup = $2; have_setup = 1 }; next }
      FNR == 1 { next }
      {
        if ($10 != 1) failed = 1
        if ($8 + 0 > residual) residual = $8 + 0
        iterations += $6
        total += $9
        count++
      }
      END {
        if (failed || !have_setup || count != expected) exit 1
        printf "%s\t%d\t%d\t%.6e\t%.3f\t%.3f\n", precision, count, iterations, residual, setup, total
      }
    ' "$result/run.tsv" "$result/columns.tsv" | tee -a "$out/timing.tsv"
  done
done
