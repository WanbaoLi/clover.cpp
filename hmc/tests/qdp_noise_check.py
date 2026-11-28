"""Compare the installed QDP47 stream with raw fields exported from native QDP++."""
import argparse
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

from replay_chroma import ChromaNoise, difference


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("runtime","draws"):
        p.add_argument(name,type=Path)
    p.add_argument("--dims",nargs=4,type=int,required=True)
    p.add_argument("--seed",nargs=4,type=int,required=True)
    p.add_argument("--uniform",type=float,required=True)
    args=p.parse_args()
    sys.path.insert(0,str(args.runtime.resolve()))
    from engine import MASSES,SEQUENCE
    from qdp_stream import QDPNoise

    tags=[name for name in SEQUENCE if name in MASSES or name=="rat_3strange"]
    actual=QDPNoise(args.dims,args.seed)
    reference=ChromaNoise(args.draws,args.dims,tags,args.uniform)
    volume=math.prod(args.dims)
    comparisons={}
    for tag in ["momentum",*tags]:
        shape=(volume,4,3,3) if tag=="momentum" else (volume//2,4,3)
        x=actual.standard_normal(shape)+1j*actual.standard_normal(shape)
        y=reference.standard_normal(shape)+1j*reference.standard_normal(shape)
        comparisons[tag]=difference(x,y)
        del x,y
    before=[int(node.text) for node in ET.parse(args.draws/"accept_rng.xml").findall("record/Seed/elem")]
    observed=actual.words
    scalar=actual.random()
    passed=(before==observed and args.uniform==scalar and
            all(c["relative_l2"]<1e-14 and c["max_abs"]<1e-13 for c in comparisons.values()))
    print(json.dumps(dict(dims=args.dims,initial_seed=args.seed,fields=comparisons,
                          pre_metropolis=observed,uniform=scalar,post_metropolis=actual.words,
                          passed=passed),indent=2,allow_nan=False))
    return 0 if passed else 1


if __name__=="__main__":
    sys.exit(main())
