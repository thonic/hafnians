"""J=2 even cat — edit CP, J, cutoff before qsub."""

from __future__ import annotations

import numpy as np
from scipy.special import factorial

CAT_ALPHA = 0.7
INPUT_NMAX = 5

STATE = "even_cat"


def _cp_even_cat(alpha: complex = CAT_ALPHA, nmax: int = INPUT_NMAX) -> np.ndarray:
    alpha = complex(alpha)
    cp = np.zeros(nmax, dtype=np.complex128)
    damp = np.exp(-0.5 * abs(alpha) ** 2)
    for n in range(nmax):
        cn = damp * (alpha ** n) / np.sqrt(float(factorial(n)))
        cn_neg = damp * ((-alpha) ** n) / np.sqrt(float(factorial(n)))
        cp[n] = cn + cn_neg
    return cp


CP = _cp_even_cat()
J = 2
CUTOFF = 6
INTERFEROMETER = "haar"
HAAR_RANDOM_SEED = 20250810
N_WORKERS: int | None = None
ZERO_TOL = 1e-8

OUTPUT_DIRNAME = "output"
RESULT_JSON = "Pnbar_J2_even_cat_LR.json"
RESULT_TXT = "Pnbar_J2_even_cat_LR.txt"
