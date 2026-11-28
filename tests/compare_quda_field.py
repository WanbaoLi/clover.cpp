"""Compare a complete CPU propagator with a mass-normalized QUDA NPY.

The NPY must contain FP64 complex values in DeGrand–Rossi basis with axes
[source spin*3+color, t, z, y, x, sink spin, sink color]. The reference and
CPU solve must use the same gauge, source position, mass and boundary
conditions. Pass the digest of the complete NPY, including its header.
Only NumPy and the sibling compare_chains.py are required.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

import numpy as np

from compare_chains import Shards, require, rows


def file_digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(8*1024*1024), b""):
            digest.update(data)
    return digest.hexdigest()


def basis_change():
    # The first factor rotates the gamma convention. gamma_2 is the
    # DeGrand-Rossi matrix used for the Bridge Dirac spin transformation.
    identity = np.eye(2)
    rotation = np.block([[identity, -identity], [identity, identity]])/math.sqrt(2)
    gamma_2 = np.array([[0,0,0,-1], [0,0,1,0], [0,1,0,0], [-1,0,0,0]])
    return np.kron(rotation @ gamma_2, np.eye(3))


def norm_by_source(matrix):
    return np.sum(np.abs(matrix)**2, axis=(0,1))


def compare(args):
    expected = args.sha256.lower()
    require(re.fullmatch("[a-f0-9]{64}", expected), "expected a 64-digit NPY SHA256")
    require(file_digest(args.reference) == expected, "reference NPY SHA256 mismatch")
    field = Shards(args.solve, None, args.endian)
    reference = np.load(args.reference, mmap_mode="r", allow_pickle=False)
    require(isinstance(reference,np.memmap) and reference.dtype.kind == "c" and reference.dtype.itemsize == 16
            and reference.shape == (12,*field.dims[::-1],4,3), "incompatible reference NPY axes, shape or precision")
    u = basis_change()
    nx,ny,nz,nt = field.local
    z_stride = max(1, 8192//(nx*ny))
    order = [field.rhs.index(rhs) for rhs in range(12)]
    error2 = np.zeros(12)
    max_abs = np.zeros(12)
    cpu_times = np.zeros((12,field.dims[3]))
    ref_times = np.zeros_like(cpu_times)
    sites = 0
    for (px,py,pz,pt),shard in sorted(field.files.items()):
        x0,y0,z0,t0 = px*nx,py*ny,pz*nz,pt*nt
        for t in range(nt):
            for z in range(0,nz,z_stride):
                z_end = min(z+z_stride,nz)
                gold = reference[:,t0+t,z0+z:z0+z_end,y0:y0+ny,x0:x0+nx]
                # A site matrix has [sink,source] spin-color indices.
                gold = np.ascontiguousarray(gold).reshape(12,-1,12).transpose(1,2,0)
                rotated = u @ gold @ u.conj().T
                cpu = shard[order,t,z:z_end].reshape(12,-1,12).transpose(1,2,0)
                delta = cpu-rotated
                error2 += norm_by_source(delta)
                max_abs = np.maximum(max_abs, np.max(np.abs(delta),axis=(0,1)))
                cpu_times[:,t0+t] += norm_by_source(cpu)
                ref_times[:,t0+t] += norm_by_source(rotated)
                sites += rotated.shape[0]
    require(sites == math.prod(field.dims), "incomplete global site coverage")
    require(all(np.all(np.isfinite(a)) for a in (error2,max_abs,cpu_times,ref_times)), "nonfinite field")
    reference2 = ref_times.sum(axis=1)
    require(np.all(reference2 > 0), "reference has a zero-norm source column")
    field_rel = np.sqrt(error2/reference2)
    time_delta = cpu_times-ref_times
    time_rel = np.linalg.norm(time_delta,axis=1)/np.linalg.norm(ref_times,axis=1)
    time_abs = np.max(np.abs(time_delta),axis=1)
    density_cpu = cpu_times.sum(axis=0)
    density_ref = ref_times.sum(axis=0)
    result = {
        "reference_npy_sha256":expected, "lattice_xyzt":field.dims, "grid_xyzt":field.grid,
        "source":field.meta["source"], "position_xyzt":field.meta["position"], "mass":float(field.meta["mass"]),
        "basis_change":"U @ S_DeGrandRossi @ U.conj().T", "normalization":field.meta["normalization"],
        "field_relative_l2":float(np.sqrt(error2.sum()/reference2.sum())), "field_max_abs":float(max_abs.max()),
        "complex_values":sites*144, "max_cpu_independent_relative":field.max_residual,
        "summed_timeslice_density_relative_l2":float(np.linalg.norm(density_cpu-density_ref)/np.linalg.norm(density_ref)),
        "summed_timeslice_density_max_abs":float(np.max(np.abs(density_cpu-density_ref))),
        "columns":[dict(rhs=i,field_relative_l2=float(field_rel[i]),field_max_abs=float(max_abs[i]),
                        timeslice_density_relative_l2=float(time_rel[i]),timeslice_density_max_abs=float(time_abs[i]))
                   for i in range(12)]}
    times_file = args.solve/"timeslices.tsv"
    if times_file.exists():
        table = np.full_like(cpu_times,np.nan)
        for row in rows(times_file):
            if row["rhs"] == "sum": continue
            i,t = int(row["rhs"]),int(row["t"])
            require(0 <= i < 12 and 0 <= t < field.dims[3] and np.isnan(table[i,t]), "invalid timeslice key")
            table[i,t] = float(row["density"])
        require(np.all(np.isfinite(table)), "incomplete timeslice table")
        result["cpu_timeslice_table_readback_relative_l2"] = float(np.linalg.norm(table-cpu_times)/np.linalg.norm(cpu_times))
        require(result["cpu_timeslice_table_readback_relative_l2"] < 1e-10, "recorded densities disagree with the CPU fields")
    result["per_column_relative_limit"] = args.rtol
    result["passed"] = bool(np.max(field_rel) < args.rtol and np.max(time_rel) < args.rtol)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solve", type=Path, help="complete MPI solver output directory")
    parser.add_argument("reference", type=Path)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--endian", choices=("little","big"), default="little", help="CPU field byte order")
    parser.add_argument("--rtol", type=float, default=1e-8, help="strict relative L2 limit per field column and its timeslice density")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(math.isfinite(args.rtol) and args.rtol > 0, "invalid relative tolerance")
    result = compare(args)
    data = json.dumps(result,indent=2,allow_nan=False)+"\n"
    if args.output: args.output.write_text(data,encoding="utf-8")
    else: print(data,end="")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
