"""J=2 superposition — edit CP, J, cutoff before qsub."""

from __future__ import annotations

import numpy as np

STATE = "superposition"
CP = np.array([1.0, 1.0], dtype=np.complex128)
J = 2
CUTOFF = 6
INTERFEROMETER = "haar"
HAAR_BASE_SEED = 20250810
HAAR_RANDOM_SEED = HAAR_BASE_SEED  # single-shot alias
N_ENSEMBLE = 10          # <-- Change ONLY this to 1000 for production
ENSEMBLE_BATCH_SIZE = 10  # realizations per PBS array task
ZERO_SAVE_TOL = 1e-10  # sparse P(n̄) storage threshold
N_WORKERS: int | None = None
ZERO_TOL = 1e-8

OUTPUT_DIRNAME = "output"
RESULT_JSON = "Pnbar_superposition_J2_parallel_LR.json"
RESULT_TXT = "Pnbar_superposition_J2_parallel_LR.txt"
