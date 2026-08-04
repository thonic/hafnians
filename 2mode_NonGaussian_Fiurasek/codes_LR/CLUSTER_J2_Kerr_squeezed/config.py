"""J=2 Kerr-squeezed state — edit state parameters before qsub."""

from __future__ import annotations

import numpy as np
from scipy.special import factorial

KERR_SQUEEZED_KAPPA = 0.1
SQUEEZE_R = float(np.arcsinh(1.0))  # sinh(r)=1, r≈0.8814
INPUT_NMAX = 5

STATE = "kerr_squeezed"


def _cp_kerr_squeezed(
    kappa: float = KERR_SQUEEZED_KAPPA,
    r: float = SQUEEZE_R,
    nmax: int = INPUT_NMAX,
) -> np.ndarray:
    """Truncated coefficients of exp(-i*kappa*n²)|S(r)>."""
    cp = np.zeros(nmax, dtype=np.complex128)
    sech_r = 1.0 / np.cosh(r)
    for n in range(0, nmax, 2):
        phase = np.exp(-1j * kappa * n * n)
        amplitude = (
            np.sqrt(float(factorial(n)) * np.tanh(r) ** n)
            / float(factorial(n // 2))
            / 2.0 ** (n / 2)
        )
        cp[n] = phase * amplitude
    return cp * np.sqrt(sech_r)


CP = _cp_kerr_squeezed()
J = 2
CUTOFF = 6
APPLY_HADAMARD = True
N_WORKERS: int | None = None
ZERO_TOL = 1e-8

OUTPUT_DIRNAME = "output"
RESULT_JSON = "Pnbar_J2_kerr_squeezed_LR.json"
RESULT_TXT = "Pnbar_J2_kerr_squeezed_LR.txt"
