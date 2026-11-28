"""Compare native CG and QXS MG on the same nontrivial links and odd RHS."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np

def relative(left,right):
    return float(np.linalg.norm(left-right)/np.linalg.norm(right))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker",type=Path,required=True)
    parser.add_argument("--mpiexec",default="mpiexec")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--dims",type=int,nargs=4,default=(16,16,4,4))
    parser.add_argument("--grid", type=int, nargs=4, default=(1,1,1,1))
    args=parser.parse_args()
    exe=args.worker.resolve()
    dims=tuple(args.dims)
    grid=tuple(args.grid)
    if any(p<1 or n%p for n,p in zip(dims,grid)):
        parser.error("MPI grid must divide each lattice extent")
    count=int(np.prod(grid))
    if count>1 and not args.mpiexec:
        parser.error("multiple ranks require an MPI launcher")
    subprocess.run([exe,"--check-mg",*(str(n//p) for n,p in zip(dims,grid))],check=True)
    rng=np.random.default_rng(701130)
    vol=int(np.prod(dims))
    noise=rng.normal(size=(vol,4,3,3))+1j*rng.normal(size=(vol,4,3,3))
    hermitian=.025*(noise+noise.conj().swapaxes(-1,-2))
    hermitian-=np.trace(hermitian,axis1=-2,axis2=-1)[...,None,None]*np.eye(3)/3
    vals,vecs=np.linalg.eigh(hermitian)
    gauge=np.ascontiguousarray((vecs*np.exp(1j*vals[...,None,:]))@vecs.conj().swapaxes(-1,-2))
    source=rng.normal(size=(vol,4,3))+1j*rng.normal(size=(vol,4,3))
    records={}
    with tempfile.TemporaryDirectory() as tmp:
        for method in ("cg","mg"):
            path=Path(tmp)/method
            path.mkdir()
            gauge.tofile(path/"u.bin");source.tofile(path/"phi.bin")
            env=dict(os.environ,CLOVER_HMC_LIGHT_SOLVER=method,OMP_NUM_THREADS="1")
            cmd=([*shlex.split(args.mpiexec),"-n",str(count)] if args.mpiexec else [])+[exe,*map(str,dims)]
            if count>1:
                cmd.extend(map(str,grid))
            proc=subprocess.run(cmd,cwd=path,env=env,input="2 -0.279 -0.2785 1e-12 1000\n0 -0.279 -0.2785 1e-12 1000\n",
                                text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            if proc.returncode: raise RuntimeError(proc.stdout[-3000:])
            rows=[line.split() for line in proc.stdout.splitlines() if line.split()]
            refresh=next(r for r in rows if r[0]=="REFRESH")
            action=next(r for r in rows if r[0]=="RATIO")
            calls=[r for r in rows if r[0]=="MG_NORMAL"]
            records[method]=dict(heatbath_residual=float(refresh[1]),normal_residual=float(action[4]),
                                 action=float(action[2]),polishing_iterations=int(action[5]),
                                 mg=[dict(stage=r[1],operator_count=int(r[2]),normal_residual_before_polish=float(r[3])) for r in calls])
            assert records[method]["heatbath_residual"]<=1.5e-12 and records[method]["normal_residual"]<=1.5e-12
            if method=="mg": assert len(calls)==2 and all(int(r[2])>0 and float(r[3])<1e-11 for r in calls)
        for file in ("refreshed.bin","force.bin"):
            arrays=[np.fromfile(Path(tmp)/method/file,dtype="=c16") for method in ("mg","cg")]
            error=relative(*arrays)
            records[file]=dict(relative_l2=error)
            assert error<(2e-8 if file=="force.bin" else 1e-10),(file,error)
    records["grid"]=grid
    records["dims"]=dims
    args.output.write_text(json.dumps(records,indent=2)+"\n")
    print(json.dumps(records))

if __name__=="__main__":
    main()
