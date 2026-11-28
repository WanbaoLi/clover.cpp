#!/usr/bin/env python3
"""Select cores and a preferred high-bandwidth NUMA node on the measured LingSheng topology."""
import argparse
import os
from pathlib import Path
import shutil
import sys


def ids(cpulist):
    values = []
    for item in cpulist.strip().split(","):
        if item:
            low, _, high = item.partition("-")
            values.extend(range(int(low), int(high or low) + 1))
    return values


def rank_slots(domains, count, width):
    """Assign whole teams within domains, spreading ranks while capacity permits."""
    if width <= 0 or count <= 0:
        raise ValueError("rank and core counts must be positive")
    slots = {node: 0 for node in sorted(domains)}
    for _ in range(count):
        available = [n for n in slots if (slots[n] + 1) * width <= len(domains[n])]
        if not available:
            raise ValueError("insufficient domain-local physical cores in the inherited CPU mask")
        node = min(available, key=lambda n: (slots[n], n))
        slots[node] += 1
    return [(n, domains[n][i * width:(i + 1) * width])
            for n, ranks in slots.items() for i in range(ranks)]


def rank_binding(rank, count, threads, node_root, allowed_cpus=None, allowed_memory=None):
    if not (0 <= rank < count and threads > 0):
        raise ValueError("invalid local rank, rank count, or thread count")
    nodes = sorted(node_root.glob("node[0-9]*"), key=lambda n: int(n.name[4:]))
    if [n.name for n in nodes] != [f"node{i}" for i in range(32)]:
        raise ValueError("expected the LingSheng topology with NUMA nodes 0 through 31")
    cpus = [ids((n / "cpulist").read_text()) for n in nodes]
    if not all(cpus[:16]) or any(cpus[16:]):
        raise ValueError("expected CPU domains 0–15 and memory-only domains 16–31")
    domains = {node: [cpu for cpu in local if allowed_cpus is None or cpu in allowed_cpus]
               for node, local in enumerate(cpus[:16])
               if allowed_memory is None or node + 16 in allowed_memory}
    node, team = rank_slots(domains, count, threads)[rank]
    return team, node + 16


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rank-cpus",type=int,help="reserve this many cores per rank, including for single-threaded workers")
    ap.add_argument("--show", action="store_true", help="print this rank's binding and exit")
    ap.add_argument("command", nargs=argparse.REMAINDER)
    args = ap.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    try:
        rank, count, threads = [int(os.environ[key]) for key in
                                ("OMPI_COMM_WORLD_LOCAL_RANK", "OMPI_COMM_WORLD_LOCAL_SIZE", "OMP_NUM_THREADS")]
        slot_size = threads if args.rank_cpus is None else args.rank_cpus
        if threads < 1 or threads > slot_size:
            raise ValueError("each rank needs at least one reserved core per OpenMP thread")
        allowed_memory = None
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("Mems_allowed_list:"):
                allowed_memory = set(ids(line.split(":", 1)[1]))
                break
        cpus, memory = rank_binding(rank, count, slot_size, Path("/sys/devices/system/node"),
                                    os.sched_getaffinity(0), allowed_memory)
        cpulist = ",".join(map(str, cpus))
        if args.show:
            print(f"local_rank={rank} cpus={cpulist} preferred_memory={memory}")
            return
        if not command:
            ap.error("supply a command after -- or use --show")
        numactl = shutil.which("numactl")
        if not numactl:
            raise ValueError("numactl is required")
        os.environ.setdefault("OMP_PLACES", "cores")
        os.environ.setdefault("OMP_PROC_BIND", "close")
        os.execv(numactl, [numactl, "--physcpubind=" + cpulist, "--preferred=" + str(memory), *command])
    except (OSError, ValueError, KeyError) as error:
        sys.exit(f"LingSheng rank binding: {error}")


if __name__ == "__main__":
    main()
