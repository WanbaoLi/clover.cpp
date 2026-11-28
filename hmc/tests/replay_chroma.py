"""Replay one native Chroma draw through the installed clover-hmc trajectory step."""
import argparse
from contextlib import closing
import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import time

import numpy as np


# U gamma_Chroma U^dagger = gamma_Bridge for all four axes.
SPIN = np.array([[0,1,0,1],[-1,0,-1,0],[0,-1,0,1],[1,0,-1,0]],complex)/math.sqrt(2)


def raw(path, shape):
    if path.stat().st_size != math.prod(shape)*16:
        raise ValueError(f"wrong field length or geometry: {path.name}")
    return np.memmap(path, dtype="<c16", mode="r", shape=shape)


def difference(candidate, reference):
    squared = []
    denom = []
    maxima = []
    if candidate.shape != reference.shape:
        raise ValueError("fields have different geometry")
    for i in range(0, len(candidate), 8192):
        y = reference[i:i+8192]
        diff = candidate[i:i+8192]-y
        squared.append(float(np.vdot(diff,diff).real))
        denom.append(float(np.vdot(y,y).real))
        maxima.append(float(np.abs(diff).max()))
    denominator = math.fsum(denom)
    if denominator <= 0:
        raise ValueError("zero reference norm")
    return dict(relative_l2=math.sqrt(math.fsum(squared)/denominator), max_abs=max(maxima))


class ChromaNoise:
    """Supply the exact array interface used by Dynamics.refresh and generator.step."""
    def __init__(self, draws, dims, tags, uniform):
        self.draws = draws
        self.volume = math.prod(dims)
        self.tags = ["momentum"]+list(tags)
        self.index = 0
        self.component = 0
        self.pending = None
        self.uniform = uniform
        self.used_uniform = False
        self.odd = np.zeros(tuple(reversed(dims)), dtype=np.uint8)
        for axis, size in enumerate(reversed(dims)):
            shape = [1]*4
            shape[axis] = size
            self.odd ^= (np.arange(size, dtype=np.uint8)%2).reshape(shape)
        self.odd = self.odd.ravel().astype(bool)

    def standard_normal(self, shape):
        tag = self.tags[self.index]
        if self.component == 0:
            if tag == "momentum":
                self.pending = raw(self.draws/"momentum_raw.bin",(self.volume,4,3,3))
            else:
                eta = raw(self.draws/(tag+"_eta.bin"),(self.volume,4,3))
                if np.any(eta[~self.odd]):
                    raise ValueError(f"nonzero even Schur noise for {tag}")
                # The QDP export already includes the Gaussian normalization used
                # by refresh(). Undo it once before passing the two real draws.
                self.pending = np.einsum("ab,vbc->vac",SPIN,eta[self.odd])/math.sqrt(.5)
            if self.pending.shape != tuple(shape):
                raise ValueError(f"unexpected generator draw at {tag}: {shape}")
            self.component = 1
            return self.pending.real
        self.component = 0
        self.index += 1
        result = self.pending.imag
        self.pending = None
        return result

    def random(self):
        if self.index != len(self.tags) or self.component or self.used_uniform:
            raise ValueError("unexpected Metropolis random draw")
        self.used_uniform = True
        return self.uniform


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("runtime","gauge","draws","output","scratch"):
        parser.add_argument(name,type=Path)
    parser.add_argument("--dims",nargs=4,required=True,type=int)
    parser.add_argument("--uniform",type=float,required=True)
    parser.add_argument("--light-solver",choices=("cg","mg"),default="cg")
    parser.add_argument("--save-heatbath",action="store_true",help="save each refreshed Bridge-basis spinor for independent comparison")
    parser.add_argument("--reference",type=Path,help="Chroma selected gauge in LIME or raw x-fast order")
    parser.add_argument("--proposal",type=Path,help="optional independent pre-Metropolis gauge")
    parser.add_argument("--momentum",type=Path,help="optional independent final momentum")
    args = parser.parse_args()
    if not 0 <= args.uniform < 1:
        raise ValueError("invalid uniform")
    sys.path.insert(0,str(args.runtime.resolve()))
    from engine import Dynamics, MASSES, PROFILE, SEQUENCE
    from generator import load_input, step
    args.output.mkdir(exist_ok=False)
    dims = tuple(args.dims)
    u = load_input(args.gauge,dims)
    tags = [t for t in SEQUENCE if t in MASSES or t=="rat_3strange"]
    rng = ChromaNoise(args.draws,dims,tags,args.uniform)
    first_p = raw(args.draws/"momentum_chroma.bin",(math.prod(dims),4,3,3))
    details = {}
    start = time.monotonic()

    class ObservedDynamics(Dynamics):
        def residuals(self,tag,rows,refresh):
            value=super().residuals(tag,rows,refresh)
            if refresh:
                if args.save_heatbath:
                    (args.output/(tag+"_phi.bin")).write_bytes(self.phi[tag])
                receipt=dict(tag=tag,max_true_residual=value,
                             phi_sha256=hashlib.sha256(self.phi[tag]).hexdigest(),
                             seconds=time.monotonic()-start)
                with (args.output/"heatbath.jsonl").open("a") as log:
                    log.write(json.dumps(receipt,allow_nan=False)+"\n")
                print("Heatbath",tag,value,flush=True)
            return value

        def refresh(self,u,rng):
            p,record = super().refresh(u,rng)
            details["initial_momentum"] = difference(p,first_p)
            if details["initial_momentum"]["relative_l2"] >= 1e-13:
                raise ValueError("Chroma/Bridge initial momentum mismatch")
            print("Completed native refresh",time.monotonic()-start,flush=True)
            return p,record

        def evolve(self,u,p,config):
            v,q = super().evolve(u,p,config)
            v.astype("<c16",copy=False).tofile(args.output/"proposal.bin")
            q.astype("<c16",copy=False).tofile(args.output/"momentum.bin")
            if args.proposal:
                details["proposal"] = difference(v,load_input(args.proposal,dims))
            if args.momentum:
                details["final_momentum"] = difference(q,raw(args.momentum,q.shape))
            return v,q

    cfg = PROFILE["default_integration"]
    with tempfile.TemporaryDirectory(dir=args.scratch) as scratch:
        with closing(ObservedDynamics(dims,Path(scratch),light_solver=args.light_solver)) as dyn:
            selected,record = step(dyn,u,dict(dims=dims,config=cfg),rng)
    if rng.index != len(rng.tags) or rng.component or rng.used_uniform != (record["uniform"] is not None):
        raise ValueError("wrong number of Chroma draws consumed")
    selected.astype("<c16",copy=False).tofile(args.output/"selected.bin")
    if args.reference:
        details["selected_gauge"] = difference(selected,load_input(args.reference,dims))
    result = dict(dims=dims,integration=cfg,light_solver=args.light_solver,monomial_draw_order=tags,
                  conditional_uniform=args.uniform,metropolis_draw_consumed=rng.used_uniform,
                  step=record,comparison=details,seconds=time.monotonic()-start)
    (args.output/"result.json").write_text(json.dumps(result,allow_nan=False,indent=2)+"\n")
    print(json.dumps(details), flush=True)


if __name__ == "__main__":
    main()
