"""Compare a completed clover-hmc checkpoint with the same saved Chroma update."""
import argparse
import json
from pathlib import Path
import sys

from replay_chroma import difference


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("runtime","chain","archive"):
        parser.add_argument(name,type=Path)
    parser.add_argument("--relative-limit",type=float,required=True)
    args=parser.parse_args()
    if not 0<args.relative_limit<float("inf"):
        parser.error("--relative-limit must be finite and positive")
    sys.path.insert(0,str(args.runtime.resolve()))
    from gauge_io import chroma_restart
    from generator import load_input,read_state,statepath
    dims,number,words=chroma_restart(args.archive)
    path=statepath(args.chain,number)
    if not path.is_file():
        parser.error(f"the chain has no completed checkpoint for archive update {number}")
    gauge,state=read_state(path)
    if tuple(state["dims"])!=dims:
        parser.error("the two checkpoint geometries differ")
    selected=load_input(args.archive,dims)
    comparison=difference(gauge,selected)
    seed_match=state["rng"]=={"bit_generator":"QDP47","words":words}
    measured=state["step"] is not None
    passed=bool(measured and seed_match and comparison["relative_l2"]<args.relative_limit)
    print(json.dumps(dict(update=number,dims=dims,comparison=comparison,chain_executed_update=measured,
                          archive_QDP_seed=words,chain_rng=state["rng"],seed_match=seed_match,
                          relative_limit=args.relative_limit,passed=passed),indent=2,allow_nan=False))
    return 0 if passed else 1


if __name__=="__main__":
    sys.exit(main())
