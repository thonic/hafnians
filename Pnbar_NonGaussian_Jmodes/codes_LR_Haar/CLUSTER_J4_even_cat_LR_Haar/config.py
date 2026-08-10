"""J=4 even cat — mode-position geometries (two active copies + vacuum)."""

from __future__ import annotations

import numpy as np
from scipy.special import factorial

# Match probability_nbar_J2/study_params.py
CAT_ALPHA = 0.7
INPUT_NMAX = 5
J = 4
CUTOFF = 6
INTERFEROMETER = "haar"
HAAR_BASE_SEED = 20250810
HAAR_RANDOM_SEED = HAAR_BASE_SEED  # single-shot alias
N_ENSEMBLE = 10          # <-- Change ONLY this to 1000 for production
ENSEMBLE_BATCH_SIZE = 10  # realizations per PBS array task
ZERO_SAVE_TOL = 1e-10  # sparse P(n̄) storage threshold
N_WORKERS: int | None = None
ZERO_TOL = 1e-8

STATE = "even_cat"
CP_VACUUM = np.array([1.0], dtype=np.complex128)

# Mode-position: active copies (1-based), vacuum elsewhere. Order A → B → C.
GEOMETRIES: dict[str, tuple[int, int]] = {
    "A": (1, 2),  # (s, s, vac, vac)
    "B": (1, 3),  # (s, vac, s, vac)
    "C": (1, 4),  # (s, vac, vac, s)
}
GEOMETRY_ORDER: tuple[str, ...] = ("A", "B", "C")

OUTPUT_DIRNAME = "output"
SUMMARY_JSON = "Pnbar_J4_even_cat_LR_all_geometries.json"
SUMMARY_TXT = "Pnbar_J4_even_cat_LR_all_geometries.txt"


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


def mode_position_states(active_modes: tuple[int, int]) -> tuple[np.ndarray, ...]:
    active = set(active_modes)
    return tuple(
        np.asarray(CP, dtype=np.complex128)
        if (j + 1) in active
        else np.asarray(CP_VACUUM, dtype=np.complex128)
        for j in range(J)
    )


def geometry_layout(active_modes: tuple[int, int]) -> str:
    active = set(active_modes)
    slots = ["s" if (j + 1) in active else "vac" for j in range(J)]
    return f"({', '.join(slots)})"


def result_json_name(geometry: str) -> str:
    return f"Pnbar_J4_even_cat_LR_geom{geometry}.json"


def result_txt_name(geometry: str) -> str:
    return f"Pnbar_J4_even_cat_LR_geom{geometry}.txt"
