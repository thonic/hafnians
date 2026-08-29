"""Low-rank loop hafnian for Cerit (self-contained under LR_Cluster_Codes)."""

from __future__ import annotations

import numpy as np

from .cluster_reps import expand_rows
from .fast_loop_hafnian import low_rank_loop_hafnian_fast, takagi_factor


def takagi_factor_b_mat(b_mat: np.ndarray, *, tol: float = 1e-11) -> np.ndarray:
    return takagi_factor(np.asarray(b_mat, dtype=complex), tol=tol)


def loop_hafnian_lr_from_G(
    G: np.ndarray,
    mu: np.ndarray,
    reps: np.ndarray | None = None,
    *,
    use_numba: bool = True,
) -> complex:
    G = np.asarray(G, dtype=complex)
    mu = np.asarray(mu, dtype=complex).reshape(-1)
    if reps is not None:
        G, mu = expand_rows(G, mu, reps)
    if G.shape[0] == 0:
        return 1.0 + 0.0j
    return low_rank_loop_hafnian_fast(G, mu, use_numba=use_numba)


def loop_hafnian_lr(
    b_mat: np.ndarray,
    mu: np.ndarray,
    reps: np.ndarray | None = None,
    *,
    G: np.ndarray | None = None,
    use_numba: bool = True,
    tol: float = 1e-11,
) -> complex:
    if G is None:
        G = takagi_factor_b_mat(b_mat, tol=tol)
    return loop_hafnian_lr_from_G(G, mu, reps, use_numba=use_numba)
