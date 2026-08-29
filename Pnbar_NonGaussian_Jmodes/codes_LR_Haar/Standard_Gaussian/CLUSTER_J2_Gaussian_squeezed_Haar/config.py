"""J=2 standard Gaussian squeezed vacuum — edit state parameters before qsub."""

from __future__ import annotations

import numpy as np
from scipy.special import factorial

SQUEEZE_R = float(np.arcsinh(1.0))  # sinh(r)=1, r≈0.8814
INPUT_NMAX = 5

STATE = "gaussian_squeezed"


def _cp_gaussian_squeezed(r: float = SQUEEZE_R, nmax: int = INPUT_NMAX) -> np.ndarray:
    """Truncated Fock coefficients of |S(r)> (even photons only)."""
    cp = np.zeros(nmax, dtype=np.complex128)
    sech_r = 1.0 / np.cosh(r)
    for n in range(0, nmax, 2):
        amplitude = (
            np.sqrt(float(factorial(n)) * np.tanh(r) ** n)
            / float(factorial(n // 2))
            / 2.0 ** (n / 2)
        )
        cp[n] = amplitude
    return cp * np.sqrt(sech_r)


CP = _cp_gaussian_squeezed()
J = 2
CUTOFF = 6
INTERFEROMETER = "haar"
HAAR_BASE_SEED = 20250810
HAAR_RANDOM_SEED = HAAR_BASE_SEED  # single-shot alias
N_ENSEMBLE = 1000          # production ensemble size
ENSEMBLE_BATCH_SIZE = 10  # realizations per PBS array task
ZERO_SAVE_TOL = 1e-10  # sparse P(n̄) storage threshold
N_WORKERS: int | None = None
ZERO_TOL = 1e-8

OUTPUT_DIRNAME = "output"
RESULT_JSON = "Pnbar_J2_gaussian_squeezed_LR.json"
RESULT_TXT = "Pnbar_J2_gaussian_squeezed_LR.txt"
