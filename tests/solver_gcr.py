"""Exercise restarted MG-GCR in both outer precisions."""
import csv
from pathlib import Path
import subprocess
import sys
import tempfile

launcher, count_flag, *solver = sys.argv[1:]
with tempfile.TemporaryDirectory(prefix="clover-gcr-") as workspace:
    for precision, tolerance in (("single", "1e-6"), ("double", "1e-9")):
        out = Path(workspace) / precision
        subprocess.run(
            [launcher, count_flag, "1", *solver, "--method", "mg", "--mg-outer", "gcr",
             "--gcr-restart", "1", "--mg-max-iters", "100", "--lattice", "8,8,4,4",
             "--grid", "1,1,1,1", "--gauge", "unit", "--mass", "1", "--rhs", "0,1",
             "--threads", "1", "--precision", precision, "--nvec", "4", "--setup-steps", "1",
             "--tol", tolerance, "--accept", tolerance, "--fields", "none", "--output", str(out)],
            check=True, timeout=60,
        )
        assert (out / "COMPLETE").is_file()
        with (out / "columns.tsv").open() as table:
            rows = list(csv.DictReader(table, delimiter="\t"))
        assert [int(row["rhs"]) for row in rows] == [0, 1]
        for row in rows:
            iterations = int(row["outer_loops"])
            # Restarting each step refreshes the residual after each operator image.
            assert iterations > 0 and int(row["nconv"]) == 1 + 2*iterations
            assert row["accepted"] == "1" and float(row["independent_relative"]) < float(tolerance)
