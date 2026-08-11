"""J=4 Kerr-squeezed state (kappa=1) — two active copies plus vacuum."""

from __future__ import annotations

import numpy as np
from scipy.special import factorial

KERR_SQUEEZED_KAPPA = 1.0
SQUEEZE_R = float(np.arcsinh(1.0))  # sinh(r)=1, r≈0.8814
INPUT_NMAX = 5
J = 4
CUTOFF = 6
INTERFEROMETER = "haar"
HAAR_BASE_SEED = 20250810
HAAR_RANDOM_SEED = HAAR_BASE_SEED  # single-shot alias
N_ENSEMBLE = 1000          # production ensemble size
ENSEMBLE_BATCH_SIZE = 10  # realizations per PBS array task
ZERO_SAVE_TOL = 1e-10  # sparse P(n̄) storage threshold
N_WORKERS: int | None = None
ZERO_TOL = 1e-8

STATE = "kerr_squeezed_kappa_1"
CP_VACUUM = np.array([1.0], dtype=np.complex128)

# Compatibility with the existing J=4 result schema, whose cat jobs expose
# an "alpha" field. Kerr-squeezed jobs write null for that in JSON.
CAT_ALPHA = None

GEOMETRIES: dict[str, tuple[int, int]] = {
    "A": (1, 2),  # (state, state, vacuum, vacuum)
    "B": (1, 3),  # (state, vacuum, state, vacuum)
    "C": (1, 4),  # (state, vacuum, vacuum, state)
}
GEOMETRY_ORDER: tuple[str, ...] = ("A", "B", "C")

OUTPUT_DIRNAME = "output"
SUMMARY_JSON = "Pnbar_J4_kerr_squeezed_kappa_1_LR_all_geometries.json"
SUMMARY_TXT = "Pnbar_J4_kerr_squeezed_kappa_1_LR_all_geometries.txt"


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


def mode_position_states(active_modes: tuple[int, int]) -> tuple[np.ndarray, ...]:
    active = set(active_modes)
    return tuple(
        np.asarray(CP, dtype=np.complex128)
        if (position + 1) in active
        else np.asarray(CP_VACUUM, dtype=np.complex128)
        for position in range(J)
    )


def geometry_layout(active_modes: tuple[int, int]) -> str:
    active = set(active_modes)
    slots = ["s" if (position + 1) in active else "vac" for position in range(J)]
    return f"({', '.join(slots)})"


def result_json_name(geometry: str) -> str:
    return f"Pnbar_J4_kerr_squeezed_kappa_1_LR_geom{geometry}.json"


def result_txt_name(geometry: str) -> str:
    return f"Pnbar_J4_kerr_squeezed_kappa_1_LR_geom{geometry}.txt"
