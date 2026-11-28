"""Check placement on full, cropped and uneven CPU and memory allocations."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "solver/examples"))
from bind_lingsheng import rank_binding
from bind_x86 import local_rank


class NumaBindings(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        for node in range(32):
            dest = self.root / f"node{node}"
            dest.mkdir()
            cpus = f"{38*node}-{38*node+37}" if node < 16 else ""
            (dest / "cpulist").write_text(cpus)

    def test_whole_application_node(self):
        allowed = {cpu for cpu in range(608) if cpu % 38 != 37}
        teams = [rank_binding(i, 16, 37, self.root, allowed) for i in range(16)]
        self.assertEqual(set().union(*(set(t) for t, _ in teams)), allowed)
        self.assertEqual([n for _, n in teams], list(range(16, 32)))
        self.assertTrue(all(len(t) == 37 for t, _ in teams))

    def test_uneven_cpus_and_allowed_hbm(self):
        allowed = {1, 3, 5, 7, 40, 41, 79, 84}
        teams = [rank_binding(i, 3, 2, self.root, allowed, {16, 18}) for i in range(3)]
        self.assertEqual(teams, [([1, 3], 16), ([5, 7], 16), ([79, 84], 18)])
        with self.assertRaises(ValueError):
            rank_binding(0, 4, 2, self.root, allowed, {16, 18})

    def test_nested_and_direct_mpi(self):
        step = {"SLURM_LOCALID": "0", "SLURM_PROCID": "0", "SLURM_NTASKS_PER_NODE": "1"}
        with patch.dict(os.environ, {**step, "PMI_SIZE": "4", "PMI_RANK": "3"}, clear=True):
            self.assertEqual(local_rank(4), (3, 4))
            with self.assertRaises(ValueError):
                local_rank(2)
        direct = {**step, "SLURM_LOCALID": "1", "SLURM_PROCID": "5", "SLURM_NTASKS_PER_NODE": "4",
                  "PMI_SIZE": "8", "PMI_RANK": "5"}
        with patch.dict(os.environ, direct, clear=True):
            self.assertEqual(local_rank(None), (1, 4))


if __name__ == "__main__":
    unittest.main()
