#!/usr/bin/env python3
"""Place one MPI rank on physical cores within one Linux NUMA domain."""
import argparse
import json
import os
from pathlib import Path
import sys

from bind_lingsheng import ids, rank_slots


def physical_domains(root, allowed):
    domains = {}
    for path in (root / "node").glob("node[0-9]*/cpulist"):
        cores = []
        seen = set()
        for cpu in sorted(allowed.intersection(ids(path.read_text()))):
            if cpu in seen:
                continue
            siblings = root / "cpu" / f"cpu{cpu}" / "topology/thread_siblings_list"
            seen.update(ids(siblings.read_text()))
            cores.append(cpu)
        if cores:
            domains[int(path.parent.name[4:])] = cores
    return domains


def local_rank(explicit_count):
    env = os.environ
    for rank_name, size_name in (("OMPI_COMM_WORLD_LOCAL_RANK", "OMPI_COMM_WORLD_LOCAL_SIZE"),
                                 ("MPI_LOCALRANKID", "MPI_LOCALNRANKS")):
        if rank_name in env:
            count = int(env[size_name]) if explicit_count is None else explicit_count
            return int(env[rank_name]), count
    # Some Hydra versions only expose the world rank. It identifies a local
    # slot exactly when the caller confirms the entire world fits on one node.
    if explicit_count and int(env.get("PMI_SIZE", "0")) == explicit_count:
        return int(env["PMI_RANK"]), explicit_count
    # A nested Hydra launcher can inherit its Slurm parent's local rank.
    if "SLURM_LOCALID" in env and ("PMI_RANK" not in env or env.get("SLURM_PROCID") == env["PMI_RANK"]):
        count = int(env["SLURM_NTASKS_PER_NODE"]) if explicit_count is None else explicit_count
        return int(env["SLURM_LOCALID"]), count
    raise ValueError("local MPI rank unavailable; with single-node Hydra supply --ranks-per-node")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rank-cpus", type=int, help="reserve this many physical cores for each rank")
    ap.add_argument("--ranks-per-node", type=int, help="local MPI rank count when the launcher omits it")
    ap.add_argument("--show", action="store_true", help="print this rank's planned CPU mask and exit")
    ap.add_argument("command", nargs=argparse.REMAINDER)
    args = ap.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    try:
        threads = int(os.environ.get("OMP_NUM_THREADS", "1"))
        width = threads if args.rank_cpus is None else args.rank_cpus
        if not 1 <= threads <= width:
            raise ValueError("each OpenMP thread requires a reserved physical core")
        rank, count = local_rank(args.ranks_per_node)
        if not 0 <= rank < count:
            raise ValueError("local rank is outside the requested range")
        domains = physical_domains(Path("/sys/devices/system"), os.sched_getaffinity(0))
        node, cpus = rank_slots(domains, count, width)[rank]
        if args.show:
            print(json.dumps(dict(local_rank=rank, numa=node, cpus=cpus)))
            return
        if not command:
            ap.error("supply a command after -- or use --show")
        os.sched_setaffinity(0, cpus)
        os.environ["OMP_PLACES"] = ",".join("{"+str(c)+"}" for c in cpus)
        os.environ["OMP_PROC_BIND"] = "close"
        os.environ["OMP_NUM_THREADS"] = str(threads)
        os.execvp(command[0], command)
    except (KeyError, OSError, ValueError) as error:
        sys.exit(f"x86 rank binding: {error}")


if __name__ == "__main__":
    main()
