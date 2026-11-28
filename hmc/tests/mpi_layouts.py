"""Compare native MPI decompositions against two ordinary serial HMC updates."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime",type=Path,required=True,help="installed libexec/clover-hmc directory")
    parser.add_argument("--report",type=Path)
    args=parser.parse_args()
    runtime=args.runtime.resolve()
    sys.path.insert(0,str(runtime))
    from generator import read_state,statepath
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",CLOVER_HMC_FERMION_THREADS="2")
    with tempfile.TemporaryDirectory(prefix="clover-hmc-grid-test-") as folder:
        root=Path(folder)
        def run(name,grid,command):
            start=command=="run"
            options=["--seed","1703","--tau",".12","--outer-steps","1","--inner-steps","1"] if start else []
            subprocess.run([sys.executable,runtime/"generator.py",command,"--output",root/name,*options],
                           env=dict(env,CLOVER_HMC_GRID=grid),check=True)
        results=[]
        for name,grids in (("serial",["1 1 1 1"]*2),("distributed",["1 1 1 2","2 1 1 1"])):
            run(name,grids[0],"run")
            run(name,grids[1],"resume")
        for n,grid in enumerate(([1,1,1,2],[2,1,1,1]),1):
            u,state=read_state(statepath(root/"distributed",n))
            v,reference=read_state(statepath(root/"serial",n))
            relative=float(np.linalg.norm(u-v)/np.linalg.norm(v))
            record=state["step"]
            expected=reference["step"]
            if not (relative<1e-12 and record["accepted"]==expected["accepted"]
                    and state["rng"]==reference["rng"] and record["native_mpi_grid"]==grid):
                raise RuntimeError(f"update {n} failed for grid {grid}, field relative L2 {relative}")
            delta=record["delta_h"]-expected["delta_h"]
            if abs(delta)>1e-8:raise RuntimeError(f"Hamiltonian change differs by {delta}")
            results.append(dict(update=n,grid=grid,field_relative_l2=relative,delta_h_difference=delta,
                                accepted=record["accepted"],same_rng=True))
    text=json.dumps(results,indent=2,allow_nan=False)
    if args.report:args.report.write_text(text+"\n")
    print(text)


if __name__=="__main__":
    main()
