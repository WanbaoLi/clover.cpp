"""Check true residuals, complete binary fields, timeslices and agreement of two methods."""
import argparse
from array import array
import csv
import json
import math
from pathlib import Path

cli=argparse.ArgumentParser()
cli.add_argument("left",type=Path)
cli.add_argument("right",type=Path,nargs="?")
args=cli.parse_args()
def read(path):
    assert (path/"COMPLETE").is_file(),path
    def rows(name):
        with (path/name).open() as f: return list(csv.DictReader(f,delimiter="\t"))
    meta={r["key"]:r["value"] for r in rows("run.tsv")}
    cols=rows("columns.tsv")
    rhs=[int(i) for i in meta["rhs"].split(",")]
    lattice=[int(i) for i in meta["lattice"].split(",")]
    grid=[int(i) for i in meta["grid"].split(",")]
    assert [int(r["rhs"]) for r in cols]==rhs
    residuals=[float(r["independent_relative"]) for r in cols]
    assert all(r<float(meta["independent_acceptance_strict"]) for r in residuals)
    ts=rows("timeslices.tsv")
    sums=[float(r["density"]) for r in ts if r["rhs"]=="sum"]
    assert len(sums)==lattice[3] and all(math.isfinite(v) and v>=0 for v in sums)
    for t,total in enumerate(sums):
        per=[float(r["density"]) for r in ts if r["rhs"]!="sum" and int(r["t"])==t]
        assert len(per)==len(rhs) and math.isclose(math.fsum(per),total,rel_tol=2e-13,abs_tol=1e-300)
    shards=sorted(path.glob("field_*.f64"))
    if meta["fields"]=="full":
        total_bytes=math.prod(lattice)*24*8*len(rhs)
        assert len(shards)==math.prod(grid) and all(f.stat().st_size==total_bytes//len(shards) for f in shards)
        reduced={int(r["rhs"])*lattice[3]+int(r["t"]):float(r["density"]) for r in ts if r["rhs"]!="sum"}
        rebuilt={k:0.0 for k in reduced}
        nx,ny,nz,nt=[l//g for l,g in zip(lattice,grid)]
        for shard in shards:
            coords=[int(v) for v in shard.stem.split("_")[1:]]
            assert len(coords)==4 and all(0<=c<g for c,g in zip(coords,grid))
            with shard.open("rb") as stream:
                for col in rhs:
                    for t in range(nt):
                        n=nx*ny*nz*24
                        chunks=[]
                        while n:
                            take=min(n,8192)
                            values=array("d");values.fromfile(stream,take)
                            chunks.append(math.fsum(v*v for v in values))
                            n-=take
                        rebuilt[col*lattice[3]+nt*coords[3]+t]+=math.fsum(chunks)
                assert not stream.read(1)
        error=max(abs(rebuilt[k]-v)/v if v else abs(rebuilt[k]) for k,v in reduced.items())
        assert error<2e-12,error
    else:
        total_bytes=0
        error=None
        assert not shards
    return meta,shards,{"method":meta["solver"],"columns":len(cols),"max_independent_relative":max(residuals),
                        "raw_plaquette":float(meta["raw_plaquette"]),"stout_plaquette":float(meta["stout_plaquette"]),
                        "operator_relative":meta["scalar_qxs_relative_l2"],"field_bytes":total_bytes,"norm2":math.fsum(sums),
                        "field_readback_max_timeslice_relative":error}
left,files,result=read(args.left)
if args.right:
    right,others,reference=read(args.right)
    for key in ("lattice","grid","position","source","mass","csw","rho","rhs"):
        assert left[key]==right[key],key
    assert files and [f.name for f in files]==[f.name for f in others]
    norms=[];errors=[]
    for a,b in zip(files,others):
        with a.open("rb") as first,b.open("rb") as second:
            while True:
                aa=array("d");bb=array("d")
                aa.frombytes(first.read(1024*64));bb.frombytes(second.read(1024*64))
                if not aa: break
                norms.append(math.fsum(x*x for x in bb))
                errors.append(math.fsum((x-y)**2 for x,y in zip(aa,bb)))
    delta=math.sqrt(math.fsum(errors)/math.fsum(norms))
    assert delta<1e-10,delta
    assert math.isclose(math.fsum(norms),reference["norm2"],rel_tol=2e-12)
    result={"first":result,"second":reference,"all_field_components_relative_l2":delta}
print(json.dumps(result,indent=2))
