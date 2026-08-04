"""
Build Walrus hafnian geometry (B, μ, scale) from Fiurášek cv/dv + interferometer.

Used once per cluster job before the parallel pattern loop in run_parallel_J2_Gaussian_squeezed_LR.py.

Pipeline:
  cp → generate_alpha_corrected → generate_cv_and_dv → rearrange → Hadamard → _hafnian_objects
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
    cp: np.ndarray,
    *,
    J: int,
    apply_hadamard: bool = True,
):
    """
    J identical Fiurášek copies of ``cp``, optional Hadamard on signal modes.

    Returns
    -------
    b_mat, mu, scale, MK_dict
    """
    cp = np.asarray(cp, dtype=np.complex128).reshape(-1)
    states = tuple([cp] * J)
    uint = (
        hadamard(J).astype(np.complex128) / np.sqrt(J)
        if apply_hadamard and J > 1
        else np.eye(J, dtype=np.complex128)
    )

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
