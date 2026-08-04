"""
Fast low-rank (loop) hafnian with Numba kernels for ranks r = 1, 2, 3, 4.

A = B @ B.T  with  B.shape = (N, r).
The loop hafnian with loop weights mu is

    lhaf(A, mu) = E_{x ~ CN(0,I_r)}[ prod_{i=1}^N (mu_i + sum_k B_ik x_k) ]

which expands to

    lhaf = sum_{c in N^r, all c_k even}  coeff_Q(c) * prod_k (c_k - 1)!!,

    Q(x) = prod_i (mu_i + sum_k B_ik x_k).

For mu = 0 the ordinary hafnian is recovered.  Odd N is allowed for the loop
hafnian (it gives 0 only for the plain hafnian).
"""

from __future__ import annotations

import numpy as np
from numpy.linalg import eigh

try:
    from numba import njit

    _HAS_NUMBA = True
except Exception:  # pragma: no cover
    _HAS_NUMBA = False

    def njit(*args, **kwargs):
        def _decorator(f):
            return f

        return _decorator


def _double_factorial_arr(m: int) -> np.ndarray:
    """(e-1)!! for even e in [0, m]; odd entries keep value 1 (unused)."""
    df = np.ones(m + 1, dtype=np.float64)
    for e in range(2, m + 1, 2):
        df[e] = df[e - 2] * (e - 1)
    return df


def takagi_factor(
    A: np.ndarray, rank: int | None = None, tol: float = 1e-11
) -> np.ndarray:
    """
    Autonne–Takagi factorization of a complex symmetric matrix:
        A = U diag(d) U.T  ->  return B = U sqrt(d) so that A = B @ B.T.
    """
    A = np.asarray(A, dtype=complex)
    A = (A + A.T) / 2.0
    N = A.shape[0]
    if N == 0:
        return np.zeros((0, 0), dtype=complex)

    P, Q = A.real, A.imag
    M = np.block([[P, Q], [Q, -P]])
    w, X = eigh(M)

    order = np.argsort(w)[::-1]
    w, X = w[order], X[:, order]
    scale = max(1.0, np.max(np.abs(w)))
    keep = w > tol * scale

    d = w[keep]
    a = X[:N, keep]
    b = X[N:, keep]
    B = (a + 1j * b) * np.sqrt(d)[None, :]

    if rank is not None:
        rank = min(rank, B.shape[1])
        B = B[:, :rank]
    return B


def _loop_hafnian_py(B: np.ndarray, mu: np.ndarray, df: np.ndarray) -> complex:
    N, r = B.shape
    base = N + 1
    strides = [base**k for k in range(r)]

    coeffs: dict[int, complex] = {0: 1.0 + 0j}
    use_mu = mu is not None

    for i in range(N):
        row = B[i]
        mi = mu[i] if use_mu else 0.0 + 0j
        new: dict[int, complex] = {}
        for key, val in coeffs.items():
            if use_mu:
                new[key] = new.get(key, 0.0 + 0j) + val * mi
            for k in range(r):
                b = row[k]
                if b != 0.0 + 0j:
                    nkey = key + strides[k]
                    new[nkey] = new.get(nkey, 0.0 + 0j) + val * b
        coeffs = new

    total = 0.0 + 0j
    for key, val in coeffs.items():
        tmp = key
        weight = 1.0
        ok = True
        for _ in range(r):
            e = tmp % base
            tmp //= base
            if e & 1:
                ok = False
                break
            weight *= df[e]
        if ok:
            total += val * weight
    return total


@njit(cache=True)
def _pairwise_sum_terms(terms: np.ndarray, n: int) -> complex:
    """Pairwise sum of terms[0:n] (n <= 5 in our kernels)."""
    while n > 1:
        j = 0
        i = 0
        while i < n:
            if i + 1 < n:
                terms[j] = terms[i] + terms[i + 1]
                i += 2
            else:
                terms[j] = terms[i]
                i += 1
            j += 1
        n = j
    return terms[0]


@njit(cache=True)
def _neumaier_sum_terms(terms: np.ndarray, n: int) -> complex:
    """Neumaier compensated sum of terms[0:n] (n <= 5)."""
    s = 0.0 + 0.0j
    c = 0.0 + 0.0j
    for i in range(n):
        t = terms[i]
        z = s + t
        if abs(s) >= abs(t):
            c += (s - z) + t
        else:
            c += (t - z) + s
        s = z
    return s + c


@njit(cache=True)
def _cell_sum_terms(terms: np.ndarray, n: int) -> complex:
    """Pairwise default; Neumaier for two-term (cancellation) and five-term cells."""
    if n == 2 or n == 5:
        return _neumaier_sum_terms(terms, n)
    if n <= 0:
        return 0.0 + 0.0j
    return _pairwise_sum_terms(terms, n)


@njit(cache=True)
def _lr_loop_haf_r1(B: np.ndarray, mu: np.ndarray, df: np.ndarray) -> complex:
    N = B.shape[0]
    C = np.zeros(N + 1, dtype=np.complex128)
    C[0] = 1.0 + 0.0j
    terms = np.zeros(5, dtype=np.complex128)

    for i in range(N):
        b = B[i, 0]
        m = mu[i]
        for e in range(i + 1, -1, -1):
            nt = 0
            if m != 0.0 + 0.0j:
                terms[nt] = m * C[e]
                nt += 1
            if e > 0:
                terms[nt] = b * C[e - 1]
                nt += 1
            C[e] = _cell_sum_terms(terms, nt) if nt > 0 else 0.0 + 0.0j

    total = 0.0 + 0.0j
    for e in range(N + 1):
        if (e & 1) == 0:
            total += C[e] * df[e]
    return total


@njit(cache=True)
def _lr_loop_haf_r2(B: np.ndarray, mu: np.ndarray, df: np.ndarray) -> complex:
    N = B.shape[0]
    C = np.zeros((N + 1, N + 1), dtype=np.complex128)
    C[0, 0] = 1.0 + 0.0j
    terms = np.zeros(5, dtype=np.complex128)

    for i in range(N):
        a = B[i, 0]
        b = B[i, 1]
        m = mu[i]
        for e0 in range(i + 1, -1, -1):
            for e1 in range(i + 1 - e0, -1, -1):
                nt = 0
                if m != 0.0 + 0.0j:
                    terms[nt] = m * C[e0, e1]
                    nt += 1
                if e0 > 0:
                    terms[nt] = a * C[e0 - 1, e1]
                    nt += 1
                if e1 > 0:
                    terms[nt] = b * C[e0, e1 - 1]
                    nt += 1
                C[e0, e1] = _cell_sum_terms(terms, nt) if nt > 0 else 0.0 + 0.0j

    total = 0.0 + 0.0j
    for e0 in range(0, N + 1, 2):
        for e1 in range(0, N + 1 - e0, 2):
            total += C[e0, e1] * df[e0] * df[e1]
    return total


@njit(cache=True)
def _lr_loop_haf_r3(B: np.ndarray, mu: np.ndarray, df: np.ndarray) -> complex:
    N = B.shape[0]
    C = np.zeros((N + 1, N + 1, N + 1), dtype=np.complex128)
    C[0, 0, 0] = 1.0 + 0.0j
    terms = np.zeros(5, dtype=np.complex128)

    for i in range(N):
        a0 = B[i, 0]
        a1 = B[i, 1]
        a2 = B[i, 2]
        m = mu[i]
        for e0 in range(i + 1, -1, -1):
            for e1 in range(i + 1 - e0, -1, -1):
                for e2 in range(i + 1 - e0 - e1, -1, -1):
                    nt = 0
                    if m != 0.0 + 0.0j:
                        terms[nt] = m * C[e0, e1, e2]
                        nt += 1
                    if e0 > 0:
                        terms[nt] = a0 * C[e0 - 1, e1, e2]
                        nt += 1
                    if e1 > 0:
                        terms[nt] = a1 * C[e0, e1 - 1, e2]
                        nt += 1
                    if e2 > 0:
                        terms[nt] = a2 * C[e0, e1, e2 - 1]
                        nt += 1
                    C[e0, e1, e2] = _cell_sum_terms(terms, nt) if nt > 0 else 0.0 + 0.0j

    total = 0.0 + 0.0j
    for e0 in range(0, N + 1, 2):
        for e1 in range(0, N + 1 - e0, 2):
            for e2 in range(0, N + 1 - e0 - e1, 2):
                total += C[e0, e1, e2] * df[e0] * df[e1] * df[e2]
    return total


@njit(cache=True)
def _lr_loop_haf_r4(B: np.ndarray, mu: np.ndarray, df: np.ndarray) -> complex:
    """Rank-4 kernel. Table size (N+1)^4 complex128 — avoid very large N."""
    N = B.shape[0]
    C = np.zeros((N + 1, N + 1, N + 1, N + 1), dtype=np.complex128)
    C[0, 0, 0, 0] = 1.0 + 0.0j
    terms = np.zeros(5, dtype=np.complex128)

    for i in range(N):
        a0 = B[i, 0]
        a1 = B[i, 1]
        a2 = B[i, 2]
        a3 = B[i, 3]
        m = mu[i]
        for e0 in range(i + 1, -1, -1):
            for e1 in range(i + 1 - e0, -1, -1):
                for e2 in range(i + 1 - e0 - e1, -1, -1):
                    for e3 in range(i + 1 - e0 - e1 - e2, -1, -1):
                        nt = 0
                        if m != 0.0 + 0.0j:
                            terms[nt] = m * C[e0, e1, e2, e3]
                            nt += 1
                        if e0 > 0:
                            terms[nt] = a0 * C[e0 - 1, e1, e2, e3]
                            nt += 1
                        if e1 > 0:
                            terms[nt] = a1 * C[e0, e1 - 1, e2, e3]
                            nt += 1
                        if e2 > 0:
                            terms[nt] = a2 * C[e0, e1, e2 - 1, e3]
                            nt += 1
                        if e3 > 0:
                            terms[nt] = a3 * C[e0, e1, e2, e3 - 1]
                            nt += 1
                        C[e0, e1, e2, e3] = (
                            _cell_sum_terms(terms, nt) if nt > 0 else 0.0 + 0.0j
                        )

    total = 0.0 + 0.0j
    for e0 in range(0, N + 1, 2):
        for e1 in range(0, N + 1 - e0, 2):
            for e2 in range(0, N + 1 - e0 - e1, 2):
                for e3 in range(0, N + 1 - e0 - e1 - e2, 2):
                    total += C[e0, e1, e2, e3] * df[e0] * df[e1] * df[e2] * df[e3]
    return total


def low_rank_loop_hafnian_fast(
    B: np.ndarray,
    mu: np.ndarray | None = None,
    *,
    use_numba: bool = True,
) -> complex:
    """Loop hafnian of A = B @ B.T."""
    B = np.asarray(B, dtype=complex)
    N, r = B.shape
    if N == 0:
        return 1.0 + 0j

    if mu is None:
        mu = np.einsum("ik,ik->i", B, B)
    else:
        mu = np.broadcast_to(np.asarray(mu, dtype=complex), (N,)).copy()

    df = _double_factorial_arr(N)

    if use_numba and _HAS_NUMBA:
        if r == 1:
            return _lr_loop_haf_r1(B, mu, df)
        if r == 2:
            return _lr_loop_haf_r2(B, mu, df)
        if r == 3:
            return _lr_loop_haf_r3(B, mu, df)
        if r == 4:
            return _lr_loop_haf_r4(B, mu, df)

    return _loop_hafnian_py(B, mu, df)


def low_rank_hafnian_fast(B: np.ndarray, *, use_numba: bool = True) -> complex:
    """Ordinary hafnian of A = B @ B.T. Returns 0 for odd N."""
    B = np.asarray(B, dtype=complex)
    N = B.shape[0]
    if N % 2 != 0:
        return 0.0 + 0j
    mu = np.zeros(N, dtype=complex)
    return low_rank_loop_hafnian_fast(B, mu=mu, use_numba=use_numba)
