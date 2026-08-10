"""
Build Walrus hafnian geometry (B, μ, scale) from Fiurášek cv/dv + interferometer.

Used once per cluster job before the parallel pattern loop in run_parallel_J2_even_cat.py.

Pipeline:
  cp → generate_alpha_corrected → generate_cv_and_dv → rearrange → Haar interferometer → _hafnian_objects
"""

from __future__ import annotations

import numpy as np
from scipy.stats import unitary_group

import config

from alpha_fiurasek import generate_alpha_corrected
from generate_cv_and_dv import (
    generate_cv_and_dv,
    generate_u_cv_and_dv_udag,
    rearrange_cv_and_dv,
)
from hafnian_batched_statistics import _hafnian_objects

def _haar_unitary(J: int, seed: int) -> np.ndarray:
    """Haar-random J×J unitary with reproducible seed."""
    rng = np.random.default_rng(seed)
    U = np.asarray(unitary_group.rvs(J, random_state=rng), dtype=np.complex128)
    identity = np.eye(J, dtype=np.complex128)
    if not np.allclose(U.conj().T @ U, identity, atol=1e-10, rtol=1e-10):
        err = float(np.max(np.abs(U.conj().T @ U - identity)))
        raise ValueError(
            f"Haar unitary failed unitarity check: max|U†U - I| = {err:.3e} "
            f"for J={J}, seed={seed}"
        )
    return U



def build_B_mu_scale(
    cp: np.ndarray,
    *,
    J: int,
    interferometer: str | None = None,
):
    """
    J identical Fiurášek copies of ``cp``, Haar-random unitary on signal modes (or identity).

    Returns
    -------
    b_mat, mu, scale, MK_dict
    """
    cp = np.asarray(cp, dtype=np.complex128).reshape(-1)
    states = tuple([cp] * J)
    mode = interferometer if interferometer is not None else config.INTERFEROMETER
    if mode == "haar" and J > 1:
        uint = _haar_unitary(J, config.HAAR_RANDOM_SEED)
    elif mode == "identity" or J <= 1:
        uint = np.eye(J, dtype=np.complex128)
    else:
        raise ValueError(f"unknown interferometer {mode!r}; use 'haar' or 'identity'")

    MK_dict: dict[int, int] = {}
    cvs, dvs, count, N = {}, {}, 0, 0
    for state_cp in states:
        state_cp = state_cp / np.sqrt(np.conj(state_cp) @ state_cp)
        m = len(state_cp)
        MK_dict[count] = m
        N += m
        alpha = generate_alpha_corrected(state_cp)
        cv, dv = generate_cv_and_dv(alpha, 1, m, m, single_mode=True)
        cvs[count], dvs[count] = cv, dv
        count += 1

    cov, disp = rearrange_cv_and_dv(cvs, dvs, count, size=N)
    cov, disp = generate_u_cv_and_dv_udag(cov, disp, MK_dict, uint)
    b_mat, mu, scale = _hafnian_objects(cov, disp)
    return b_mat, mu, scale, MK_dict
