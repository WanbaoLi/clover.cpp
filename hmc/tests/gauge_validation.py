"""Exercise both SU(3) acceptance boundaries and a final partial batch."""
from pathlib import Path
import sys
try:
    import numpy as np
except ModuleNotFoundError:
    sys.exit(77)

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"runtime"))
from generator import checked_gauge

dims=[1,1,1,2050]
gauge=np.broadcast_to(np.eye(3,dtype=np.complex128),(2050,4,3,3)).copy()
for phase,scale in ((0,1),(1e-11,1),(0,1+5e-12),(1e-4,1)):
    link=gauge[-1,-1]
    link[0,0]=scale*np.exp(1j*phase)
    err=max(np.abs(gauge.conj().swapaxes(-1,-2)@gauge-np.eye(3)).max(),
            np.abs(np.linalg.det(gauge)-1).max())
    for tol in (1e-11,np.nextafter(err,np.inf),err):
        try:
            checked_gauge(gauge,dims,tol)
            accepted=True
        except ValueError:
            accepted=False
        assert accepted==(err<tol), (phase,scale,tol,err)
gauge[-1,-1,0,0]=np.nan
try:
    checked_gauge(gauge,dims,1e-11)
except ValueError:
    pass
else:
    raise AssertionError("nonfinite final link accepted")
