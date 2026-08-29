"""J=4 superposition — mode-position geometries (two active copies + vacuum)."""

from __future__ import annotations

import numpy as np

STATE = "superposition"
CP = np.array([1.0, 1.0], dtype=np.complex128)
CP_VACUUM = np.array([1.0], dtype=np.complex128)
INPUT_NMAX = 2
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

# Compatibility with the existing J=4 result schema (cat jobs expose "alpha").
CAT_ALPHA = None

# Mode-position: active copies (1-based), vacuum elsewhere. Order A → B → C.
GEOMETRIES: dict[str, tuple[int, int]] = {
    "A": (1, 2),  # (s, s, vac, vac)
    "B": (1, 3),  # (s, vac, s, vac)
    "C": (1, 4),  # (s, vac, vac, s)
}
GEOMETRY_ORDER: tuple[str, ...] = ("A", "B", "C")

OUTPUT_DIRNAME = "output"
SUMMARY_JSON = "Pnbar_J4_superposition_LR_all_geometries.json"
SUMMARY_TXT = "Pnbar_J4_superposition_LR_all_geometries.txt"


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
    return f"Pnbar_J4_superposition_LR_geom{geometry}.json"


def result_txt_name(geometry: str) -> str:
    return f"Pnbar_J4_superposition_LR_geom{geometry}.txt"
