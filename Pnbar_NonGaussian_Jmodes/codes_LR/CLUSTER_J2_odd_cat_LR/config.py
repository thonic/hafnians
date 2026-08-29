"""J=2 odd cat — edit CP, J, cutoff before qsub."""

from __future__ import annotations

import numpy as np
from scipy.special import factorial

CAT_ALPHA = 0.7
ODD_CAT_NMAX = 6  # even length for Fiurasek (|0>..|5>)

STATE = "odd_cat"


def _cp_odd_cat(alpha: complex = CAT_ALPHA, nmax: int = ODD_CAT_NMAX) -> np.ndarray:
    alpha = complex(alpha)
    cp = np.zeros(nmax, dtype=np.complex128)
    damp = np.exp(-0.5 * abs(alpha) ** 2)
    for n in range(nmax):
        cn = damp * (alpha ** n) / np.sqrt(float(factorial(n)))
        cn_neg = damp * ((-alpha) ** n) / np.sqrt(float(factorial(n)))
        cp[n] = cn - cn_neg
    return cp


CP = _cp_odd_cat()
J = 2
CUTOFF = 6
APPLY_HADAMARD = True
N_WORKERS: int | None = None
ZERO_TOL = 1e-7  # export floor: LR (0,0) renorm P ~5e-8 vs Walrus ~1e-31 (see investigate_j2_odd_cat.txt)
RAW_SKIP_NBARS: tuple[tuple[int, ...], ...] = ((0, 0),)

OUTPUT_DIRNAME = "output"
RESULT_JSON = "Pnbar_J2_odd_cat_LR.json"
RESULT_TXT = "Pnbar_J2_odd_cat_LR.txt"
