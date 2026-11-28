"""Compare two chain-to-propagator results on global lattice coordinates.

The right-hand result defines the denominator for relative L2 errors.
Fields contain 12 mass-normalized, complex128 propagator columns in Bridge
Dirac basis. The gauge reader verifies ILDG/SciDAC extents and its data
checksum, when present. Gauge fields fit in memory; propagator shards are
memory mapped and compared one intersection and time slice at a time.
"""
import argparse
import csv
import hashlib
import importlib.util
from itertools import product
import json
import math
from pathlib import Path
import re

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def gauge_reader(runtime):
    source = runtime / "gauge_io.py"
    require(source.is_file(), "missing gauge_io.py; specify --hmc-runtime")
    spec = importlib.util.spec_from_file_location("clover_gauge_io", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.read_lime


def rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


class Difference:
    def __init__(self):
        self.d2 = 0.
        self.b2 = 0.
        self.maximum = 0.
        self.count = 0

    def add(self, left, right):
        require(left.shape == right.shape, "overlap shape mismatch")
        diff = left - right
        d2 = float(np.vdot(diff, diff).real)
        b2 = float(np.vdot(right, right).real)
        require(math.isfinite(d2) and math.isfinite(b2), "nonfinite field or norm")
        self.count += right.size
        self.d2 += d2
        self.b2 += b2
        self.maximum = max(self.maximum, float(np.max(np.abs(diff))))

    def result(self):
        require(self.count > 0 and self.b2 > 0, "empty or zero-norm reference field")
        return {"complex_values": self.count, "relative_l2": math.sqrt(self.d2 / self.b2),
                "max_abs": self.maximum}


def vector(text):
    return tuple(int(x) for x in text.split(","))


class Shards:
    def __init__(self, directory, source, endian):
        if source is not None:
            directory = directory / source
        metadata = rows(directory / "run.tsv")
        self.meta = {r["key"]: r["value"] for r in metadata}
        require(len(metadata) == len(self.meta), "duplicate key in run.tsv")
        require(self.meta["fields"] == "full" and self.meta["source"] in ("point","wall") and
                (source is None or self.meta["source"] == source), "missing full source field")
        require(self.meta["basis"] == "Bridge Dirac" and self.meta["normalization"] == "2*kappa",
                "unsupported field basis or normalization")
        self.dims = vector(self.meta["lattice"])
        self.grid = vector(self.meta["grid"])
        self.rhs = vector(self.meta["rhs"])
        require(len(self.rhs) == 12 and sorted(self.rhs) == list(range(12)), "expected all 12 unique source columns")
        require(len(self.grid) == len(self.dims) == 4 and
                all(d > 0 and g > 0 and d % g == 0 for d, g in zip(self.dims, self.grid)), "invalid MPI grid")
        cols = rows(directory / "columns.tsv")
        require(tuple(int(r["rhs"]) for r in cols) == self.rhs, "column table differs from on-disk column order")
        tol = float(self.meta["independent_acceptance_strict"])
        require(all(int(r["spin"])*3+int(r["color"]) == int(r["rhs"]) and int(r["accepted"]) == 1 and
                    0 <= float(r["independent_relative"]) < tol for r in cols), "source column failed independent acceptance")
        self.max_residual = max(float(r["independent_relative"]) for r in cols)
        self.local = tuple(d//g for d, g in zip(self.dims, self.grid))
        shape = (12, *self.local[::-1], 4, 3)
        dtype = "<c16" if endian == "little" else ">c16"
        self.files = {}
        for path in directory.glob("field_*.f64"):
            match = re.fullmatch(r"field_(\d+)_(\d+)_(\d+)_(\d+)\.f64", path.name)
            require(match, "invalid propagation shard name")
            coords = tuple(map(int, match.groups()))
            require(coords not in self.files and all(c < g for c, g in zip(coords, self.grid)), "invalid shard coordinate")
            require(path.stat().st_size == math.prod(shape)*16, "incomplete propagator shard")
            self.files[coords] = np.memmap(path, dtype=dtype, mode="r", shape=shape)
        require(len(self.files) == math.prod(self.grid), "incomplete MPI shard set")


def compare_propagators(left, right, dims, source, endians):
    a, b = (Shards(path, source, endian) for path, endian in zip((left,right), endians))
    require(a.dims == b.dims == dims, "gauge and propagator extents differ")
    require(vector(a.meta["position"]) == vector(b.meta["position"]), "physical source positions differ")
    for key in ("mass", "kappa", "csw", "rho"):
        require(float(a.meta[key]) == float(b.meta[key]), "physics parameter differs: "+key)
    norms = [Difference() for _ in range(12)]
    for coords, field in a.files.items():
        origin = tuple(x*n for x,n in zip(coords,a.local))
        end = tuple(x+n for x,n in zip(origin,a.local))
        overlaps = [range(lo//n, (hi-1)//n+1) for lo,hi,n in zip(origin,end,b.local)]
        for target in product(*overlaps):
            other = b.files[target]
            start = tuple(x*n for x,n in zip(target,b.local))
            lo = tuple(max(i,j) for i,j in zip(origin,start))
            hi = tuple(min(i,j+n) for i,j,n in zip(end,start,b.local))
            asi = tuple(slice(lo[i]-origin[i],hi[i]-origin[i]) for i in (2,1,0))
            bsi = tuple(slice(lo[i]-start[i],hi[i]-start[i]) for i in (2,1,0))
            for rhs,norm in enumerate(norms):
                for t in range(lo[3],hi[3]):
                    norm.add(field[(a.rhs.index(rhs),t-origin[3],*asi)],
                             other[(b.rhs.index(rhs),t-start[3],*bsi)])
    volume = math.prod(dims)
    require(all(n.count == volume*12 for n in norms), "incomplete global coordinate coverage")
    return {"grid_left": a.grid, "grid_right": b.grid, "source_position": vector(a.meta["position"]),
            "mass": float(a.meta["mass"]), "basis": a.meta["basis"], "normalization": a.meta["normalization"],
            "max_solver_residual_left": a.max_residual, "max_solver_residual_right": b.max_residual,
            "relative_l2": math.sqrt(math.fsum(n.d2 for n in norms)/math.fsum(n.b2 for n in norms)),
            "max_abs": max(n.maximum for n in norms), "complex_values": volume*144,
            "columns": [dict(rhs=i, **n.result()) for i,n in enumerate(norms)]}


def read_selected(directory, reader):
    state = json.loads((directory/"trajectory.json").read_text())
    dims = tuple(state["dims"])
    require(len(dims) == 4 and min(dims) > 0, "invalid checkpoint dimensions")
    field = reader(directory/("cfg-%d.lime" % state["number"]), dims)
    selected = hashlib.sha256(np.ascontiguousarray(field,dtype="<c16").view(np.uint8)).hexdigest()
    require(selected == state["gauge_sha256"], "exported gauge differs from checkpoint")
    return state,field


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("left",type=Path)
    p.add_argument("right",type=Path,help="reference for all relative L2 norms")
    p.add_argument("--labels",nargs=2,metavar=("LEFT","RIGHT"),default=("left","right"),
                   help="display names of the two results")
    p.add_argument("--hmc-runtime",type=Path,default=Path(__file__).resolve().parents[1]/"hmc"/"runtime")
    p.add_argument("--left-endian",choices=("little","big"),default="little")
    p.add_argument("--right-endian",choices=("little","big"),default="little")
    p.add_argument("--output",type=Path,help="write JSON here; default is stdout")
    args = p.parse_args()
    read = gauge_reader(args.hmc_runtime)
    (a,u),(b,v) = (read_selected(path,read) for path in (args.left,args.right))
    for key in ("dims","profile","config","seed","number"):
        require(a[key] == b[key], "chain setting differs: "+key)
    require(a["step"]["initial_sha256"] == b["step"]["initial_sha256"], "starting gauge fields differ")
    gauge = Difference()
    for offset in range(0,len(u),4096):
        gauge.add(u[offset:offset+4096],v[offset:offset+4096])
    summary = {"left_label":args.labels[0],"right_label":args.labels[1],"relative_l2_reference":"right",
               "lattice":a["dims"], "seed":a["seed"],"trajectory":a["number"],
               "profile_sha256":a["profile"],"config":a["config"],"initial_gauge_sha256":a["step"]["initial_sha256"],
               "rng_state_equal_after_update":a["rng"] == b["rng"],
               "updates":{side:{"accepted":s["step"]["accepted"],"delta_h":s["step"]["delta_h"],"uniform":s["step"]["uniform"],
                                "selected_gauge_sha256":s["gauge_sha256"]} for side,s in (("left",a),("right",b))},
               "gauge":gauge.result()}
    del u,v
    for source in ("point","wall"):
        summary[source] = compare_propagators(args.left,args.right,tuple(a["dims"]),source,(args.left_endian,args.right_endian))
    result=json.dumps(summary,indent=2,allow_nan=False)+"\n"
    if args.output: args.output.write_text(result,encoding="utf-8")
    else: print(result,end="")


if __name__=="__main__":
    main()
