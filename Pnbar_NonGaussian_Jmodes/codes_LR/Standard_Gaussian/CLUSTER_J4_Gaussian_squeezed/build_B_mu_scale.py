"""
Build Walrus hafnian geometry (B, μ, scale) from Fiurášek cv/dv + interferometer.

Used once per geometry before the parallel pattern loop in run_parallel_J4_Gaussian_squeezed_LR.py.

Supports mixed mode-position layouts: active state copies + vacuum copies.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import hadamard

from alpha_fiurasek import generate_alpha_corrected
from generate_cv_and_dv import (
    generate_cv_and_dv,
    generate_u_cv_and_dv_udag,
    rearrange_cv_and_dv,
)
from hafnian_batched_statistics import _hafnian_objects


def build_B_mu_scale(
    states: tuple[np.ndarray, ...],
    *,
    J: int | None = None,
    apply_hadamard: bool = True,
):
    """
    Fiurášek copies ``states`` (one cp array per copy), optional Hadamard on signals.

    Returns
    -------
    b_mat, mu, scale, MK_dict
    """
    states = tuple(np.asarray(s, dtype=np.complex128).reshape(-1) for s in states)
    J = int(J if J is not None else len(states))
    if len(states) != J:
        raise ValueError(f"len(states)={len(states)} != J={J}")

    uint = (
        hadamard(J).astype(np.complex128) / np.sqrt(J)
        if apply_hadamard and J > 1
        else np.eye(J, dtype=np.complex128)
    )

    MK_dict: dict[int, int] = {}
    cvs: dict[int, np.ndarray] = {}
    dvs: dict[int, np.ndarray] = {}
    N = 0
    for count, state_cp in enumerate(states):
        state_cp = state_cp / np.sqrt(np.conj(state_cp) @ state_cp)
        m = len(state_cp)
        MK_dict[count] = m
        N += m
        alpha = generate_alpha_corrected(state_cp)
        cv, dv = generate_cv_and_dv(alpha, 1, m, m, single_mode=True)
        cvs[count], dvs[count] = cv, dv

    cov, disp = rearrange_cv_and_dv(cvs, dvs, J, size=N)
    cov, disp = generate_u_cv_and_dv_udag(cov, disp, MK_dict, uint)
    b_mat, mu, scale = _hafnian_objects(cov, disp)
    return b_mat, mu, scale, MK_dict
