"""Run real native HMC twice from the same random stream and verify trajectory-boundary restart."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import tempfile
import numpy as np

def state(folder,n):
    with np.load(folder/f"trajectory-{n:08d}.npz",allow_pickle=False) as f:
        return f["gauge"],json.loads(str(f["metadata"]))
def main():
    p=argparse.ArgumentParser();p.add_argument("--clover-hmc",required=True,type=Path);args=p.parse_args()
    exe=str(args.clover_hmc.resolve())
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    def run(*argv):subprocess.run([exe,*map(str,argv)],check=True,env=env)
    with tempfile.TemporaryDirectory(prefix="clover-hmc-test-") as tmp:
        base=Path(tmp);split=base/"split";plain=base/"plain"
        options=("--dims",2,2,2,2,"--seed",1703,"--tau",.12,"--outer-steps",1,"--inner-steps",1)
        run("run","--output",split,*options)
        run("resume","--output",split)
        run("run","--output",plain,*options,"--trajectories",2)
        for n in (1,2):
            a,am=state(split,n);b,bm=state(plain,n)
            if am!=bm or not np.array_equal(a,b):
                raise RuntimeError(f"restart differs from uninterrupted trajectory {n}")
        lime=base/"update.lime";run("export","--output",split,"--lime",lime)
        # A fresh chain starting from the exported field must store identical trajectory zero.
        run("run","--output",base/"roundtrip",*options,"--gauge",lime)
        loaded,_=state(base/"roundtrip",0)
        if not np.array_equal(loaded,a):raise RuntimeError("ILDG readback mismatch")
        print("PASS: both native trajectories, restored RNG and ILDG readback agree.")
if __name__=="__main__":main()
