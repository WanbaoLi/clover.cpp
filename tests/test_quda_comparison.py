"""Synthetic independent inverse transform with an asymmetric source/sink field."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


def tsv(path,header,data):
    with path.open("w",newline="") as stream:
        writer=csv.writer(stream,delimiter="\t")
        writer.writerow(header)
        writer.writerows(data)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp); folder=root/"unrelated-output-name";folder.mkdir()
        dims=(4,3,2,4); grid=(2,1,1,2); lx,ly,lz,lt=(d//g for d,g in zip(dims,grid))
        spin_pairs=(((1,1),(3,-1)),((0,-1),(2,1)),((1,1),(3,1)),((0,-1),(2,-1)))
        rng=np.random.default_rng(73610)
        # [source_spin,source_color, t,z,y,x, sink_spin,sink_color]
        shape=(4,3,*dims[::-1],4,3)
        known=rng.normal(size=shape)+1j*rng.normal(size=shape)
        dgr=np.zeros_like(known)
        # Explicit index sums for inverse, independent of NumPy matrix
        # construction and both matrix multiplications in the production script.
        for dst in range(4):
            for src in range(4):
                for feed_dst,sd in spin_pairs[dst]:
                    for feed_src,ss in spin_pairs[src]:
                        dgr[src,...,dst,:] += .5*sd*ss*known[feed_src,...,feed_dst,:]
        known=known.reshape(12,*dims[::-1],4,3)
        ref=root/"gold.npy"
        np.save(ref,dgr.reshape(12,*dims[::-1],4,3))
        digest=hashlib.sha256(ref.read_bytes()).hexdigest()
        order=list(range(11,-1,-1))
        tsv(folder/"run.tsv",("key","value"),
            {"fields":"full","source":"point","rhs":",".join(map(str,order)),
             "lattice":"4,3,2,4","grid":"2,1,1,2","mass":"-0.279","position":"1,2,0,1",
             "basis":"Bridge Dirac","normalization":"2*kappa","independent_acceptance_strict":"1e-12"}.items())
        tsv(folder/"columns.tsv",("rhs","spin","color","accepted","independent_relative"),
            [(r,r//3,r%3,1,1e-13) for r in order])
        def write_fixtures():
            for px in range(2):
                for pt in range(2):
                    roi=known[order,pt*lt:(pt+1)*lt,:,:,px*lx:(px+1)*lx]
                    roi.astype(">c16").tofile(folder/f"field_{px}_0_0_{pt}.f64")
            density=np.sum(abs(known)**2,axis=(2,3,4,5,6))
            tsv(folder/"timeslices.tsv",("rhs","t","density"),[(i,t,density[i,t]) for i in range(12) for t in range(4)])
        def run(expected):
            result=subprocess.run([sys.executable,str(Path(__file__).with_name("compare_quda_field.py")),str(folder),
                                   str(ref),"--sha256",digest,"--endian","big","--rtol","1e-12"],capture_output=True,text=True)
            assert result.returncode==expected,(result.returncode,result.stderr)
            return json.loads(result.stdout)
        write_fixtures()
        exact=run(0)
        assert max(c["field_relative_l2"] for c in exact["columns"])<8e-16
        assert exact["summed_timeslice_density_relative_l2"]<8e-16 and exact["field_max_abs"]<3e-15
        # Touch one off-diagonal sink/source element on a global shard boundary.
        rhs=7; offset=.125+.5j
        old=known[rhs,2,0,1,2,1,2]
        known[rhs,2,0,1,2,1,2]=old+offset
        write_fixtures()
        one=run(1); column=one["columns"][rhs]
        assert abs(column["field_max_abs"]-abs(offset))<2e-15
        change=abs(abs(old+offset)**2-abs(old)**2)
        assert math.isclose(column["timeslice_density_max_abs"],change,abs_tol=2e-12)
        assert all(s["field_relative_l2"]<8e-16 for s in one["columns"] if s["rhs"]!=rhs)
        info={"lattice":dims,"grid":grid,"byte_order":"big","source_order":"reversed",
              "exact_max_field_relative_l2":max(s["field_relative_l2"] for s in exact["columns"]),
              "exact_max_field_absolute":exact["field_max_abs"],"whole_timeslice_relative_l2":exact["summed_timeslice_density_relative_l2"],
              "perturbed_rhs":rhs,"injected_abs":abs(offset),"observed_field_max_abs":column["field_max_abs"],
              "expected_timeslice_density_change":change,"measured_timeslice_density_change":column["timeslice_density_max_abs"],
              "control_exit":0,"perturbed_exit":1}
        if args.output: args.output.write_text(json.dumps(info,indent=2)+"\n")
        print(json.dumps(info))

if __name__=="__main__":main()
