"""
Photon-number amplitudes from (Σ, Δ) via Walrus hafnians (Hamilton–Jex).

Two engines for the P(n̄) numerator:
  • dense       — hafnian_batched (one tensor, all patterns)
  • per_pattern — loop_hafnian once per output pattern (same math, no big tensor)

P(n̄) uses |amplitude|² with herald modes fixed at photon number 1 in the
Fiurášek mode layout.  This is NOT conditional probability (no / P_herald).
"""

from __future__ import annotations

import cmath
from math import factorial, prod

import numpy as np
from thewalrus import hafnian_batched, loop_hafnian
from thewalrus.quantum import Xmat


def slice_probabilities(hafnian, n):
    """Index hafnian tensor for photon pattern n = (n1, n2, …)."""
    idx = []
    for i in range(2 * len(n)):
        idx += [n[i % len(n)]]
    return hafnian[tuple(idx)]


def post_select_on_herald_modes(hafnian, n, MK_dict):
    """
    |amplitude|² for system pattern n with one photon in each herald mode.

    MK_dict[k] = number of modes in copy k (signal first, then herald(s)).
    """
    idx = []
    k_len = len(MK_dict)
    for k in MK_dict.keys():
        m = MK_dict[k]
        idx += [n[k % k_len]] + [1] * (m - 1)
    return hafnian[tuple(idx)] * np.conj(hafnian[tuple(idx)])


def reps_for_nbar(nbar: tuple[int, ...], MK_dict: dict) -> np.ndarray:
    """Repetition vector matching post_select_on_herald_modes indexing."""
    reps: list[int] = []
    j = len(nbar)
    for k in MK_dict.keys():
        m = MK_dict[k]
        reps.append(int(nbar[k % j]))
        reps.extend([1] * (m - 1))
    return np.asarray(reps, dtype=np.int64)


def probability(covariance_matrix, displacement_vector, cutoff):
    """Gaussian photon-number probabilities (optional herald path; not used for P(n̄))."""
    n_modes = len(covariance_matrix) // 2
    q = covariance_matrix + np.identity(2 * n_modes) / 2
    inv_q = np.linalg.inv(q)
    a_mat = Xmat(n_modes) @ (np.identity(2 * n_modes) - inv_q)
    replacement_vector = inv_q @ displacement_vector

    detq = np.linalg.det(q)
    exp_factor = (
        -0.5
        * np.transpose(np.conj(displacement_vector))
        @ inv_q
        @ displacement_vector
    )
    exp_norm = cmath.exp(exp_factor)

    hafnians = hafnian_batched(a_mat, cutoff, mu=replacement_vector, renorm=False)
    hafnians = hafnians / cmath.sqrt(detq)
    hafnians = hafnians * exp_norm
    return hafnians


def _hafnian_objects(covariance_matrix, displacement_vector):
    """Build B, μ, and amplitude scale for the non-Gaussian P(n̄) path."""
    n_modes = len(covariance_matrix) // 2
    q = covariance_matrix + np.identity(2 * n_modes) / 2
    inv_q = np.linalg.inv(q)
    a_mat = Xmat(n_modes) @ (np.identity(2 * n_modes) - inv_q)
    b_mat = a_mat[0:n_modes, 0:n_modes]
    mu = (inv_q @ displacement_vector)[:n_modes]
    detq = np.linalg.det(q)
    exp_factor = (
        -0.5
        * np.transpose(np.conj(displacement_vector))
        @ inv_q
        @ displacement_vector
    )
    scale = cmath.sqrt(cmath.exp(exp_factor)) / cmath.sqrt(cmath.sqrt(detq))
    return b_mat, mu, scale


def loop_hafnian_renorm(b_mat, mu, reps) -> complex:
    """loop_hafnian / √(∏ kᵢ!) — matches hafnian_batched(..., renorm=True)."""
    raw = loop_hafnian(b_mat, D=mu, reps=reps, glynn=True)
    denom = cmath.sqrt(prod(factorial(int(r)) for r in np.asarray(reps).tolist()))
    return raw / denom


def non_gaussian_probability(
    covariance_matrix,
    displacement_vector,
    cutoff,
    *,
    engine: str = "dense",
):
    """
    Non-Gaussian photon-number amplitudes (P(n̄) path).

    engine="dense" returns a cutoff^n_modes tensor (hafnian_batched).
    For per-pattern evaluation use amplitude_per_pattern / driver_code.
    """
    if engine != "dense":
        raise ValueError(
            f"engine={engine!r}: use driver_code.find_probabilities(..., engine='per_pattern') "
            "for per-pattern evaluation; this function only builds the dense tensor."
        )
    b_mat, mu, scale = _hafnian_objects(covariance_matrix, displacement_vector)
    hafnians = hafnian_batched(b_mat, cutoff, mu=mu, renorm=True)
    return hafnians * scale


def amplitude_per_pattern(covariance_matrix, displacement_vector, nbar, MK_dict) -> complex:
    """One-pattern amplitude (Walrus loop_hafnian). Does not use Craig low-rank."""
    b_mat, mu, scale = _hafnian_objects(covariance_matrix, displacement_vector)
    reps = reps_for_nbar(nbar, MK_dict)
    return loop_hafnian_renorm(b_mat, mu, reps) * scale
