"""Row expansion for Walrus-compatible ``reps`` (cluster-local copy)."""

from __future__ import annotations

import numpy as np


def reps_to_row_indices(reps: np.ndarray) -> np.ndarray:
    idx: list[int] = []
    for i, r in enumerate(np.asarray(reps, dtype=int).tolist()):
        idx.extend([i] * int(r))
    return np.asarray(idx, dtype=np.int64)


def expand_rows(
    G: np.ndarray, mu: np.ndarray, reps: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    idx = reps_to_row_indices(reps)
    if idx.size == 0:
        return np.zeros((0, G.shape[1]), dtype=complex), np.zeros(0, dtype=complex)
    return G[idx, :].copy(), mu[idx].copy()
