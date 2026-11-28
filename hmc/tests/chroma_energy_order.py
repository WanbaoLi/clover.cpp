"""Reproduce the original 24³×64 Chroma energy scalars from all twenty monomial values."""
import json
from pathlib import Path
import sys

try:
    import numpy
except ImportError:
    sys.exit(77)

sys.path.insert(0,sys.argv[1])
from generator import PROFILE,chroma_energy_difference

v=json.loads(Path(sys.argv[2]).read_text())
reference=v["original_chroma_fields"]
names=PROFILE["hamiltonian"]
assert names==reference["monomial_order"]
start={"action":dict(zip(names,reference["initial_monomials"])),"kinetic":reference["initial_kinetic"]}
end={"action":dict(zip(names,reference["final_monomials"])),"kinetic":reference["final_kinetic"]}
diff=chroma_energy_difference(start,end,{})
for key in ("initial_potential","final_potential","delta_h"):
    if diff[key]!=reference[key]:
        raise AssertionError(f"{key}: expected {reference[key].hex()}, got {diff[key].hex()}")
