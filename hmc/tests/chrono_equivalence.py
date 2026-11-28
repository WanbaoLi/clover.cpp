"""Compare every light action and force under cold and chronological starts on nearby fixed SU(3) fields."""
import argparse
import json
import math
import sys
import tempfile
from contextlib import closing
from pathlib import Path
import numpy as np

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--runtime",type=Path,required=True,help="installed libexec/clover-hmc")
    parser.add_argument("--report",type=Path)
    args=parser.parse_args();sys.path.insert(0,str(args.runtime.resolve()))
    from engine import Native, MASSES,SOLVERS,project
    rng=np.random.Generator(np.random.PCG64(24011996));dims=(4,4,4,4);shape=(256,4,3,3)
    def algebra():
        return project(rng.standard_normal(shape)+1j*rng.standard_normal(shape))
    def exponential(a):
        vals,evecs=np.linalg.eigh(1j*a)
        return (evecs*np.exp(-1j*vals)[...,None,:])@evecs.conj().swapaxes(-1,-2)
    u0=exponential(.12*algebra());direction=algebra();fields=[exponential(i*.001*direction)@u0 for i in range(5)]
    odd=(np.indices((4,4,4,4)).sum(axis=0).ravel()%2)==1
    spin=(256,4,3);phis={k:(rng.standard_normal(spin)+1j*rng.standard_normal(spin))/math.sqrt(2.) for k in MASSES}
    for v in phis.values():v[~odd]=0
    result={}
    for tol in (1e-12,1e-14):
        group={k:[] for k in MASSES};result[str(tol)]=group
        with tempfile.TemporaryDirectory() as folder,closing(Native("light",dims,Path(folder))) as w:
            for i,u in enumerate(fields):
                for tag,(mass,den) in MASSES.items():
                    phis[tag].tofile(w.folder/"phi.bin")
                    opt=dict(mass=mass,den=9. if den is None else den,tol=tol,limit=SOLVERS[tag][1])
                    def solve(mode):
                        rows=w.call(u,mode=mode,**opt);r=next(x for x in rows if x[0]=="RATIO")
                        return rows,float(r[2]),float(r[4]),int(r[5]),w.gauge("force.bin")
                    _,a,ra,ia,fa=solve(0);rows,b,rb,ib,fb=solve(13)
                    c=next(x for x in rows if x[0]=="CHRONO")
                    absolute=float(np.linalg.norm(fa-fb));relative=absolute/float(np.linalg.norm(fa))
                    if max(ra,rb)>1.5*tol or relative>3e-8 or abs(a-b)/abs(a)>1e-10:
                        raise RuntimeError(f"cold/chrono comparison failed for {tag} at gauge {i}")
                    group[tag].append(dict(field=i,cold_iterations=ia,chrono_iterations=ib,
                          cold_true_residual=ra,chrono_true_residual=rb,used=int(c[2]),starting_residual=float(c[3]),
                          action_abs_difference=abs(a-b),force_absolute_L2=absolute,force_relative_L2=relative))
    if args.report:args.report.write_text(json.dumps(result,indent=2)+"\n")
    for t,group in result.items():
        allrows=[r for rows in group.values() for r in rows]
        print(t,"cold iterations",sum(r["cold_iterations"] for r in allrows),
              "chrono iterations",sum(r["chrono_iterations"] for r in allrows),
              "max true residual",max(max(r["cold_true_residual"],r["chrono_true_residual"]) for r in allrows),
              "max force relative L2",max(r["force_relative_L2"] for r in allrows))
if __name__=="__main__":main()
