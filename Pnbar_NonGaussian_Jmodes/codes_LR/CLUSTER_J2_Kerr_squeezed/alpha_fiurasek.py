"""Minimal Fiurášek α from Fock coefficients (corrected). No sympy."""

from __future__ import annotations

import math

import numpy as np


def _setup_ladder(n: int):
    nm = n - 1
    maxn = 2 * n
    am = np.zeros((maxn + 2, maxn + 2))
    ad = np.zeros((maxn + 2, maxn + 2))
    for j in range(maxn + 1):
        am[j, j + 1] = np.sqrt(j + 1)
        ad[j + 1, j] = np.sqrt(j + 1)
    sq = np.arcsinh(1.0)
    a_op = np.cosh(sq) * am + ad
    am_pow = {j: np.linalg.matrix_power(a_op, j)[:n, :n] for j in range(n + 1)}
    cvac = np.zeros(n)
    cvac[0] = 1.0
    return am_pow, cvac, nm


def generate_alpha_corrected(cp: np.ndarray, t: float = 0.99999999) -> np.ndarray:
    """α from Fock coefficients (corrected; matches production pipeline)."""
    n = len(cp)
    am_pow, cvac, nm = _setup_ladder(n)
    h = np.zeros(n, dtype=np.complex128)
    psi = np.copy(cp)
    for j in range(n):
        h[j] = psi[n - j - 1] / np.sqrt(math.factorial(nm - j))
        psi = psi - h[j] * (am_pow[n - j - 1] @ cvac)
    beta = np.roots(h)
    m = np.zeros((n - 1, n - 1))
    for j in range(n - 1, 0, -1):
        for k in range(n - 1, j - 1, -1):
            m[j - 1, k - 1] = t ** ((n - 1) - k)
    alpha = np.zeros(n, dtype=np.complex128)
    if n > 1:
        alpha[1:] = np.linalg.inv(m) @ beta
    sq = np.arcsinh(1.0)
    s1 = s2 = 0.0
    for j in range(1, n):
        s1 += alpha[j] * t ** (n - j)
        s2 += np.conj(alpha[j]) * t ** (j - n)
    s1 *= np.cosh(sq)
    sdiff = (s2 - s1) / np.cosh(sq)
    alpha[0] = (
        np.real(sdiff) / (t ** n - t ** (-n) / np.cosh(sq))
        + 1j * np.imag(sdiff) / (t ** n + t ** (-n) / np.cosh(sq))
    )
    return alpha
