"""
Ordinary low-rank hafnian of A = G @ G.T (The Walrus semantics).

This is the quantity computed by ``thewalrus.low_rank_hafnian``
(Björklund Appendix C) — **not** the loop hafnian.

Production path
---------------
``thewalrus.low_rank_hafnian`` uses a SymPy monomial expansion that is
correct (for thewalrus >= 0.22) but far too slow for Fiurášek pattern
sweeps.  We therefore evaluate the same ordinary hafnian with the
project's Numba kernel ``low_rank_hafnian_fast`` (loop hafnian with
μ ≡ 0), which matches ``thewalrus.hafnian(G @ G.T)`` and the corrected
Algorithm-C implementation.

thewalrus 0.21.0 bug
--------------------
``factorial2(2*pi - 1)`` is 0 when ``pi == 0``, so stock
``low_rank_hafnian`` returns 0 for most rank > 1 matrices.  Fixed in
0.22 via ``extend="complex"``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_LR_ROOT = Path(__file__).resolve().parent.parent
if str(_LR_ROOT) not in sys.path:
    sys.path.insert(0, str(_LR_ROOT))

from lowrank.fast_loop_hafnian import low_rank_hafnian_fast  # noqa: E402


def low_rank_hafnian(G: np.ndarray) -> complex:
    """
    Hafnian of the low-rank matrix A = G @ G.T.

    Parameters
    ----------
    G :
        Factor matrix of shape (n, r). Odd n → 0; empty n → 1.
    """
    G = np.asarray(G, dtype=complex)
    if G.ndim != 2:
        raise ValueError(f"G must be 2-D, got shape {G.shape}")
    n = G.shape[0]
    if n == 0:
        return 1.0 + 0.0j
    if n % 2 != 0:
        return 0.0 + 0.0j
    return complex(low_rank_hafnian_fast(G))


def verify_against_walrus_hafnian(*, seed: int = 0) -> None:
    """Raise if the production kernel disagrees with ``thewalrus.hafnian``."""
    from thewalrus import hafnian

    rng = np.random.default_rng(seed)
    G = rng.normal(size=(6, 3)) + 1j * rng.normal(size=(6, 3))
    A = G @ G.T
    if not np.allclose(low_rank_hafnian(G), hafnian(A), rtol=1e-8, atol=1e-10):
        raise RuntimeError("low_rank_hafnian disagrees with thewalrus.hafnian(G@G.T)")
