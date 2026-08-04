"""
Joint heralded calculation: numerator (⟨n̂⟩) and denominator (P_herald).

One file for the Fiurášek heralded pipeline — numerator and P_herald belong
to the same conditioned expectation:

    ⟨n̂_sys | herald⟩ / P_herald   ←  numerator / P_herald  gives mean photon number

Folder map (Pnbar_NonGaussian_Jmodes/)
-----------------------------------------
  run_joint_heraldonly.py        ← numerator + P_herald only
  probability_zero_check.py      ← numerator + P_herald + Pr(n̄)
  driver_code.py                 ← find_probabilities (same as hafnians/)

Quick start
-----------
    from code_joint_herald_fiurasek import compute_joint, generate_sigma_D

    joint = compute_joint(j=1)    # numerator ≈ 2e-8, P_herald ≈ 4e-8
    print(joint)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
import sympy as sp
from scipy.linalg import block_diag, hadamard
from thewalrus import loop_hafnian

from generate_cv_and_dv import generate_cv_and_dv
from generate_displacements import generate_alpha


# =============================================================================
# §1  State generation: cp → α → Σ, D  (from hafnians/, 3 files)
# =============================================================================

PDF_EXPECTED = 2.00e-8
PDF_EXPECTED_P_HERALD = 4.0e-8

# --- Default target: (|0⟩ + |1⟩) / √2  (PDF primary state, 4×4 Σ) ---
# Alternate target (secondary, 6×6 Σ, two herald modes — use j=1 only):
#   CP_SECONDARY = np.array([0.0, 1.0, 1.0], dtype=np.complex128)  # (|1⟩+|2⟩)/√2
#   st = generate_sigma_D(cp=CP_SECONDARY / np.linalg.norm([0,1,1]))


def default_target_state() -> np.ndarray:
    """(|0⟩ + |1⟩) / √2 — PDF primary state."""
    cp = np.array([1.0, 1.0], dtype=np.complex128)
    return cp / np.linalg.norm(cp)


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
    """α from Fock coefficients (corrected; matches show_matrices.py)."""
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


@dataclass
class StateMatrices:
    sigma: np.ndarray  # (2M)×(2M), order (s, h, …, s*, h*, …)
    d: np.ndarray
    alpha: np.ndarray
    cp: np.ndarray


def generate_sigma_D(
    cp: np.ndarray | None = None,
    *,
    corrected_alpha: bool = True,
    M: int | None = None,
) -> StateMatrices:
    """
    Build single-copy Σ and D from normalised Fock coefficients ``cp``.

    Parameters
    ----------
    cp : array, optional
        Fock amplitudes ``[c_0, c_1, …, c_{n-1}]``.
        Default: ``(|0⟩ + |1⟩)/√2``  →  ``cp = [1, 1]/√2``,  ``M = 2``.
    M : int, optional
        Number of optical modes in the Fiurásek circuit.
        Default: ``M = len(cp)``  (required: ``M = n`` for max Fock ``|n-1⟩``).

    Examples
    --------
        # Primary PDF state (|0⟩+|1⟩)/√2  →  Σ is 4×4
        st = generate_sigma_D()

        # (|1⟩+|2⟩)/√2  →  cp = [0,1,1]/√2,  M = 3,  Σ is 6×6
        st = generate_sigma_D(cp=[0, 1, 1])

        # (|0⟩+|2⟩)/√2  →  cp = [1,0,1]/√2,  M = 3,  Σ is 6×6
        st = generate_sigma_D(cp=[1, 0, 1])
    """
    if cp is None:
        cp = default_target_state()
    else:
        cp = np.asarray(cp, dtype=np.complex128)
        cp = cp / np.linalg.norm(cp)

    if M is None:
        M = len(cp)

    alpha = generate_alpha_corrected(cp) if corrected_alpha else generate_alpha(cp)
    sigma, d = generate_cv_and_dv(alpha, K=1, M=M, N=M, single_mode=True)
    return StateMatrices(
        sigma=np.asarray(sigma, dtype=np.complex128),
        d=np.asarray(d, dtype=np.complex128).reshape(-1),
        alpha=alpha,
        cp=cp,
    )


def pdf_reference_state() -> tuple[np.ndarray, np.ndarray]:
    """Shortcut: Σ, D from generate_sigma_D()."""
    st = generate_sigma_D()
    return st.sigma, st.d


# =============================================================================
# §2  j-copy network: direct sum → (S,H,S*,H*) → optional Hadamard
# =============================================================================

def is_power_of_two(j: int) -> bool:
    """True if j is 1, 2, 4, 8, … (valid for a Hadamard interferometer)."""
    return j >= 1 and (j & (j - 1)) == 0


def valid_hadamard_copies(max_j: int = 32) -> list[int]:
    """Allowed copy counts when ``apply_hadamard=True``: 1, 2, 4, 8, …"""
    out = [1]
    k = 2
    while k <= max_j:
        out.append(k)
        k *= 2
    return out


def permutation_to_SHSH(j: int) -> list[int]:
    """
    Index map from direct-sum copy order to grouped (S, H, S*, H*).

    Direct-sum order (before permute): for each copy k = 0…j−1,
    ``(s_k, h_k, s_k*, h_k*)`` at indices ``4k, 4k+1, 4k+2, 4k+3``.

    After permute, ``v = (s_1,…,s_j, h_1,…,h_j, s_1*,…,s_j*, h_1*,…,h_j*)``
    (1-based labels; code uses 0-based copy index k ↔ mode k+1).
    """
    return (
        [4 * k for k in range(j)]
        + [4 * k + 1 for k in range(j)]
        + [4 * k + 2 for k in range(j)]
        + [4 * k + 3 for k in range(j)]
    )


def describe_basis_order(j: int) -> list[str]:
    """Variable names in global order (S, H, S*, H*)."""
    s = [f"s{k + 1}" for k in range(j)]
    h = [f"h{k + 1}" for k in range(j)]
    sc = [f"s{k + 1}*" for k in range(j)]
    hc = [f"h{k + 1}*" for k in range(j)]
    return s + h + sc + hc


def _block_slices(j: int) -> dict[str, slice]:
    """Row/col slices for v = (S, H, S*, H*) with length 4j."""
    return {
        "S": slice(0, j),
        "H": slice(j, 2 * j),
        "S*": slice(2 * j, 3 * j),
        "H*": slice(3 * j, 4 * j),
    }


def _hadamard_unitary(j: int) -> np.ndarray:
    """Normalized Hadamard on j signal modes (j must be a power of 2, or j=1)."""
    if j == 1:
        return np.array([[1.0]], dtype=np.complex128)
    if not is_power_of_two(j):
        raise ValueError(
            f"j={j}: Hadamard network requires j ∈ {{1, 2, 4, 8, …}}, got non-power-of-two"
        )
    return hadamard(j).astype(np.complex128) / np.sqrt(j)


def build_hadamard_big(j: int) -> np.ndarray:
    """
    Full interferometer unitary on v = (S, H, S*, H*):

        U_big = diag(U_S, I_j, U_S*, I_j)
    """
    h = _hadamard_unitary(j)
    return block_diag(h, np.eye(j), h.conj(), np.eye(j))


def build_network_sigma_D(
    sigma_single: np.ndarray,
    d_single: np.ndarray,
    j: int,
    *,
    apply_hadamard: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build the j-copy covariance Σ and displacement D in global order (S,H,S*,H*).

    Steps
    -----
    1. Direct sum: ``Σ_pair = σ ⊕ σ ⊕ …`` (j times), ``D_pair = (d, d, …)``.
    2. Permute copy order ``(s_k,h_k,s_k*,h_k*)`` → grouped
       ``(s_1…s_j, h_1…h_j, s_1*…s_j*, h_1*…h_j*)``.
    3. If ``apply_hadamard``: ``Σ → U_big Σ U_big†``, ``D → U_big D`` with
       Hadamard only on signal modes (heralds unchanged).

    Parameters
    ----------
    sigma_single : (4, 4) array
        Single-copy Σ in order ``(s, h, s*, h*)`` — one signal + one herald.
    d_single : length-4 array
        Matching displacement.
    j : int
        Number of identical copies (interferometer modes).
    apply_hadamard : bool
        If True, j must be 1, 2, 4, 8, …
    """
    if j < 1:
        raise ValueError("j must be >= 1")
    sigma_single = np.asarray(sigma_single, dtype=np.complex128)
    d_single = np.asarray(d_single, dtype=np.complex128).reshape(-1)
    if sigma_single.shape != (4, 4):
        raise ValueError(
            f"sigma_single must be 4×4 (one signal + one herald per copy); "
            f"got {sigma_single.shape}. For n>2 Fock targets the full Σ is "
            f"(2n)×(2n); j-copy network applies to the primary 4×4 copy only."
        )
    if d_single.size != 4:
        raise ValueError("d_single must have length 4")

    if apply_hadamard and j > 1 and not is_power_of_two(j):
        raise ValueError(
            f"j={j}: Hadamard network requires j ∈ {{1, 2, 4, 8, …}}; "
            f"valid values: {valid_hadamard_copies()}"
        )

    # Step 1 — direct sum in copy order
    sigma = block_diag(*([sigma_single] * j))
    d = np.tile(d_single, j)

    # Step 2 — permutation to (S, H, S*, H*)
    perm = permutation_to_SHSH(j)
    p = np.eye(4 * j, dtype=np.complex128)[:, perm]
    sigma = p.conj().T @ sigma @ p
    d = p.conj().T @ d

    # Step 3 — Hadamard on signal modes only
    if apply_hadamard:
        u = build_hadamard_big(j)
        sigma = u @ sigma @ u.conj().T
        d = u @ d

    return sigma, d


def inspect_network(
    sigma_single: np.ndarray,
    d_single: np.ndarray,
    j: int,
    *,
    apply_hadamard: bool = True,
) -> dict:
    """
    Print-ready summary: permutation, basis order, and key Σ blocks.

    Verifies for j=2 (primary state) that Σ_{SH*} = η I₂ before Hadamard and
    η U_S after, matching the analytic note in Fiurasek_Theory_and_Algorithms.tex.
    """
    sigma_single = np.asarray(sigma_single, dtype=np.complex128)
    eta = complex(sigma_single[0, 3])
    eps = complex(sigma_single[0, 2])

    sig_before, d_before = build_network_sigma_D(
        sigma_single, d_single, j, apply_hadamard=False
    )
    sig_after, d_after = build_network_sigma_D(
        sigma_single, d_single, j, apply_hadamard=apply_hadamard
    )
    sl = _block_slices(j)
    h = _hadamard_unitary(j) if apply_hadamard else None

    sh_star_before = sig_before[sl["S"], sl["H*"]]
    sh_star_after = sig_after[sl["S"], sl["H*"]]

    checks = {}
    if j >= 1:
        checks["perm_length"] = len(permutation_to_SHSH(j)) == 4 * j
    if j >= 2 and abs(eta) > 0:
        checks["SH*_before_is_eta_I"] = np.allclose(
            sh_star_before, eta * np.eye(j), rtol=1e-4, atol=1e-10
        )
        if apply_hadamard and h is not None:
            checks["SH*_after_is_eta_US"] = np.allclose(
                sh_star_after, eta * h, rtol=1e-4, atol=1e-10
            )

    return {
        "j": j,
        "apply_hadamard": apply_hadamard,
        "valid_hadamard_j": valid_hadamard_copies(),
        "permutation_indices": permutation_to_SHSH(j),
        "basis_order": describe_basis_order(j),
        "sigma_shape": sig_after.shape,
        "D_order": d_after,
        "block_slices": {k: (v.start, v.stop) for k, v in sl.items()},
        "SH_star_before": sh_star_before,
        "SH_star_after": sh_star_after,
        "eta": eta,
        "epsilon": eps,
        "checks_passed": checks,
    }


# =============================================================================
# §3  SymPy variables and u_total
# =============================================================================

def _phase_symbols(j: int) -> dict[str, Any]:
    """Create s_k, h_k, zeta_k and their conjugate symbols for k=1..j."""
    s = [sp.Symbol(f"s{k}", real=True) for k in range(1, j + 1)]
    sc = [sp.Symbol(f"s{k}c", real=True) for k in range(1, j + 1)]
    h = [sp.Symbol(f"h{k}", real=True) for k in range(1, j + 1)]
    hc = [sp.Symbol(f"h{k}c", real=True) for k in range(1, j + 1)]
    z = [sp.Symbol(f"zeta{k}", real=True) for k in range(1, j + 1)]
    zc = [sp.Symbol(f"zeta{k}c", real=True) for k in range(1, j + 1)]
    return {"j": j, "s": s, "sc": sc, "h": h, "hc": hc, "z": z, "zc": zc}


def _z_order(sym: dict) -> list:
    """Z = (s_1…s_j, s_1*…s_j*, ζ_1…ζ_j, ζ_1*…ζ_j*) — loop-matrix ordering."""
    return sym["s"] + sym["sc"] + sym["z"] + sym["zc"]


def _symplectic_J(j: int) -> sp.Matrix:
    """Symplectic matrix for v = (S, H, S*, H*): pairs s_k↔s_k*, h_k↔h_k*."""
    n = 4 * j
    jmat = sp.zeros(n)
    for k in range(j):
        jmat[k, 2 * j + k] = 1
        jmat[2 * j + k, k] = -1
        jmat[j + k, 3 * j + k] = 1
        jmat[3 * j + k, j + k] = -1
    return jmat


def _build_u_total(sigma: sp.Matrix, d: sp.Matrix, sym: dict) -> sp.Expr:
    """
    Total exponent before herald integration:

        u_total = −½ v†Σv + v^T J D + ½|S|² − ½|H|² + H†ζ − ζ†H + |ζ|²
    """
    j = sym["j"]
    v = sp.Matrix(sym["s"] + sym["h"] + sym["sc"] + sym["hc"])
    vdag = sp.Matrix(sym["sc"] + sym["hc"] + sym["s"] + sym["h"]).T

    u_quad = -sp.Rational(1, 2) * (vdag * sigma * v)[0]
    u_disp = (v.T * _symplectic_J(j) * d)[0]
    u_order = sp.Rational(1, 2) * sum(sym["sc"][k] * sym["s"][k] for k in range(j))
    u_order -= sp.Rational(1, 2) * sum(sym["hc"][k] * sym["h"][k] for k in range(j))
    u_zeta = sum(
        sym["hc"][k] * sym["z"][k]
        - sym["zc"][k] * sym["h"][k]
        + sym["zc"][k] * sym["z"][k]
        for k in range(j)
    )
    return sp.expand(u_quad + u_disp + u_order + u_zeta)


# =============================================================================
# §4  Herald block extraction and Gaussian integration
# =============================================================================

def _extract_herald_block(u_total: sp.Expr, sym: dict) -> dict:
    """
    Split u_total into herald quadratic + linear + remainder:

        u_total = −h† M_h h + C†h + h†B + R,   h = (H, H*)
    """
    j, h, hc = sym["j"], sym["h"], sym["hc"]
    h0 = {x: 0 for x in h + hc}

    # Quadratic herald blocks
    a = sp.zeros(j)
    x = sp.zeros(j)
    y = sp.zeros(j)
    for i in range(j):
        for k in range(j):
            a[i, k] = sp.simplify(-sp.diff(sp.diff(u_total, hc[i]), h[k]))
            x[i, k] = sp.simplify(-sp.diff(sp.diff(u_total, hc[i]), hc[k]))
            y[i, k] = sp.simplify(-sp.diff(sp.diff(u_total, h[i]), h[k]))

    quad = sp.Integer(0)
    for i in range(j):
        for k in range(j):
            quad -= a[i, k] * hc[i] * h[k]
            quad -= sp.Rational(1, 2) * x[i, k] * hc[i] * hc[k]
            quad -= sp.Rational(1, 2) * y[i, k] * h[i] * h[k]

    u_lin = sp.expand(u_total - quad)
    b = sp.Matrix([[sp.diff(u_lin, hs)] for hs in hc + h]).subs(h0)
    c_dag = sp.Matrix([[sp.diff(u_lin, hv)] for hv in h + hc]).T.subs(h0)
    r_expr = sp.expand(u_lin.subs(h0))
    m_h = sp.Matrix.vstack(
        sp.Matrix.hstack(a, x / 2),
        sp.Matrix.hstack(y / 2, a.T),
    )
    return {"M_h": sp.simplify(m_h), "B": b, "C_dag": c_dag, "R": r_expr}


def _integrate_heralds(herald: dict, j: int) -> dict:
    """
    Integrate out herald modes (Convention A):

        u_eff = R + C† M_h^{-1} B
        gamma_measure = π^j / (2^j √det M_h)
    """
    m_h = herald["M_h"]
    m_inv = sp.simplify(m_h.inv())
    det_m = sp.factor(m_h.det())
    u_eff = sp.expand(herald["R"] + (herald["C_dag"] * m_inv * herald["B"])[0])
    gamma = sp.pi ** j / (2 ** j * sp.sqrt(det_m))
    return {"u_eff": u_eff, "gamma_measure": gamma, "det_M_h": det_m, "M_h_inv": m_inv}


# =============================================================================
# §5  Observable exponent f  (what the loop hafnian acts on)
# =============================================================================

def _build_f_from_u_eff(u_eff: sp.Expr, sym: dict) -> tuple[sp.Expr, sp.Expr]:
    """
    Build (const, f) for the observable after herald integration.

    Adds |ζ|² from the generating function e^{|ζ|²}, then splits:

        u_obs = const + f,   f|_{Z=0} = 0

    Loop hafnian / derivatives act on exp(f), NOT on exp(u_eff) directly.
    """
    j = sym["j"]
    zeta_sq = sum(sym["z"][k] * sym["zc"][k] for k in range(j))
    u_obs = sp.expand(u_eff + zeta_sq)
    z0 = {v: 0 for v in _z_order(sym)}
    const = sp.simplify(u_obs.subs(z0))
    f = sp.expand(u_obs - const)
    return const, f


def _build_f_pdf_single_copy(sigma: sp.Matrix, d: sp.Matrix, sym: dict) -> tuple[sp.Expr, sp.Expr]:
    """
    j=1 only: exponent f using PDF tilde coefficients (lines 1000–1504).

    Matches the boxed PDF result 2.00×10⁻⁸ exactly for pdf_reference_state().
    """
    s00, s01, s02, s03 = sigma[0, 0], sigma[0, 1], sigma[0, 2], sigma[0, 3]
    s10, s11, s12, s13 = sigma[1, 0], sigma[1, 1], sigma[1, 2], sigma[1, 3]
    s20, s22 = sigma[2, 0], sigma[2, 2]
    s30, s31, s33 = sigma[3, 0], sigma[3, 1], sigma[3, 3]
    ds, dh, dsc, dhc = d[0], d[1], d[2], d[3]

    m = sp.Rational(1, 2) * (s11 + s33 + 1)
    eps = s31
    den = m ** 2 - eps ** 2

    const = -(m * (dh ** 2 + dhc ** 2) / 2 + eps * (dh ** 2 + dhc ** 2) / 2) / den
    h_t = m * (s01 ** 2 + s03 ** 2) / den - sp.Rational(1, 2) * (s00 + s22 - 1)
    j_t = -m / den + 1

    lt_s = ((m * s03 + eps * s01) * dh - (m * s01 + eps * s03) * dhc) / den + dsc
    lt_sc = ((m * s01 + eps * s03) * dh - (m * s03 + eps * s01) * dhc) / den - ds
    lt_z = (m * dhc + eps * dh) / den
    lt_zc = (m * dh + eps * dhc) / den
    lt_sz = -(m * s30 + eps * s01) / den
    lt_scz = -(m * s10 + eps * s03) / den
    lt_szc = (m * s01 - eps * s30) / den
    lt_sczc = (m * s03 - eps * s10) / den

    s, sc, z, zc = sym["s"][0], sym["sc"][0], sym["z"][0], sym["zc"][0]
    f = sp.expand(
        h_t * s * sc + j_t * z * zc
        + lt_s * s + lt_sc * sc + lt_z * z + lt_zc * zc
        + lt_sz * s * z + lt_scz * sc * z + lt_szc * s * zc + lt_sczc * sc * zc
    )
    return sp.simplify(const), f


def _extract_c_L_Q(f: sp.Expr, z_vars: list) -> tuple[sp.Expr, sp.Matrix, sp.Matrix]:
    """From f = L^T Z + ½ Z^T Q Z (zero constant), extract L and Q at Z=0."""
    z0 = {v: 0 for v in z_vars}
    l_vec = sp.Matrix([[sp.diff(f, v).subs(z0)] for v in z_vars])
    n = len(z_vars)
    q = sp.zeros(n)
    for i in range(n):
        for k in range(n):
            q[i, k] = sp.simplify(sp.diff(sp.diff(f, z_vars[i]), z_vars[k]))
    return sp.Integer(0), l_vec, q


def _loop_hafnian_P(l_vec: sp.Matrix, q: sp.Matrix, j: int, n: int) -> complex:
    """P = lhaf(A_lp) with diagonal L, off-diagonal Q (thewalrus convention)."""
    q_num = np.array(
        [[complex(q[i, k].evalf()) for k in range(q.cols)] for i in range(q.rows)],
        dtype=np.complex128,
    )
    l_num = np.array([complex(l_vec[i].evalf()) for i in range(l_vec.rows)], dtype=np.complex128)
    reps = [n] * (4 * j)
    return complex(loop_hafnian(A=q_num, D=l_num, reps=reps))


# =============================================================================
# §6  Public API
# =============================================================================

@dataclass
class NumeratorResult:
    """All outputs from one compute(j) call."""

    j: int
    n: int
    gamma_measure: float
    const: float
    P: complex
    numerator: float
    """Primary result: −P for j=1 PDF check; gamma_measure·exp(const)·P otherwise."""

    det_M_h: float

    def __str__(self) -> str:
        return (
            f"j={self.j}, n={self.n}\n"
            f"  gamma_measure = {self.gamma_measure:.6g}\n"
            f"  const         = {self.const:.6g}\n"
            f"  P = lhaf      = {self.P:.6g}\n"
            f"  numerator     = {self.numerator:.6g}"
        )


def compute(
    j: int,
    n: int = 1,
    *,
    sigma_single: np.ndarray | None = None,
    d_single: np.ndarray | None = None,
    apply_hadamard: bool | None = None,
) -> NumeratorResult:
    """
    Compute the heralded numerator for j interferometer modes.

    Parameters
    ----------
    j : int
        Number of modes (change this to run j=1, 2, 4, …).
    n : int
        Photons per mode in ⟨n,…,n| · |n,…,n⟩.
    sigma_single, d_single : array-like, optional
        Single-copy 4×4 Σ and length-4 D.  Default: PDF reference state.
    apply_hadamard : bool, optional
        Apply Hadamard on signals.  Default: False for j=1, True for j>1.

    Returns
    -------
    NumeratorResult
    """
    if j < 1:
        raise ValueError("j must be >= 1")
    if apply_hadamard is None:
        apply_hadamard = j > 1
    if apply_hadamard and j > 1 and not is_power_of_two(j):
        raise ValueError(
            f"j={j}: Hadamard network requires j ∈ {{1, 2, 4, 8, …}}; "
            f"use apply_hadamard=False or pick j from {valid_hadamard_copies()}"
        )

    if sigma_single is None or d_single is None:
        sigma_single, d_single = pdf_reference_state()
    else:
        sigma_single = np.asarray(sigma_single, dtype=np.complex128)
        d_single = np.asarray(d_single, dtype=np.complex128).reshape(-1)

    sigma_np, d_np = build_network_sigma_D(
        sigma_single,
        d_single,
        j,
        apply_hadamard=apply_hadamard,
    )

    sym = _phase_symbols(j)
    sigma = sp.Matrix(sigma_np.tolist())
    d = sp.Matrix(d_np.tolist())

    # --- Stage 1: u_total ---
    u_total = _build_u_total(sigma, d, sym)

    # --- Stage 2: integrate herald modes ---
    herald = _extract_herald_block(u_total, sym)
    integrated = _integrate_heralds(herald, j)

    # --- Stage 3: observable exponent f ---
    use_pdf_f = j == 1 and not apply_hadamard
    if use_pdf_f:
        const, f = _build_f_pdf_single_copy(
            sp.Matrix(np.asarray(sigma_single, dtype=np.complex128).tolist()),
            sp.Matrix(np.asarray(d_single, dtype=np.complex128).tolist()),
            sym,
        )
    else:
        const, f = _build_f_from_u_eff(integrated["u_eff"], sym)

    # --- Stage 4: loop hafnian replaces partial derivatives ---
    z_vars = _z_order(sym)
    _, l_vec, q_mat = _extract_c_L_Q(f, z_vars)
    P = _loop_hafnian_P(l_vec, q_mat, j, n)

    # --- Stage 5: assemble numerator (never double-count exp(C†M⁻¹B)) ---
    gm = float(integrated["gamma_measure"].evalf())
    const_f = float(const.evalf())
    det_f = float(integrated["det_M_h"].evalf())
    exp_const = np.exp(const_f)

    if use_pdf_f:
        # PDF line 2091: prefactor ≈ −1, numerator ≈ −P
        numerator = float(-P.real if abs(P.imag) < 1e-15 else -P)
    else:
        # General j: full Convention A
        numerator = float((gm * exp_const * P).real)

    return NumeratorResult(
        j=j,
        n=n,
        gamma_measure=gm,
        const=const_f,
        P=P,
        numerator=numerator,
        det_M_h=det_f,
    )


# =============================================================================
# §7  P_herald (denominator of the same joint heralded calculation)
# =============================================================================

# Herald indices in single-copy v = (s, h, s*, h*)
_H_IDX = 1
_HC_IDX = 3


@dataclass
class PHeraldResult:
    """P_herald: herald success probability (denominator)."""

    p_integral: float
    p_formula: float
    delta: float
    d_h_sq: float
    epsilon: float

    def __str__(self) -> str:
        return (
            f"P_herald (integral) = {self.p_integral:.6e}\n"
            f"P_herald (formula)  = {self.p_formula:.6e}  [δ + |d_h|²]\n"
            f"  δ = {self.delta:.6e}   |d_h|² = {self.d_h_sq:.6e}   ε = {self.epsilon:.6e}"
        )


def extract_herald_params(
    sigma: np.ndarray,
    d: np.ndarray,
    *,
    h_idx: int = _H_IDX,
    hc_idx: int = _HC_IDX,
) -> tuple[float, float, complex, complex, float]:
    """Herald block from Σ, D: (δ, ε, d_h, d_h*, M) with M = 1 + δ."""
    sigma = np.asarray(sigma, dtype=np.complex128)
    d = np.asarray(d, dtype=np.complex128).reshape(-1)
    delta = float(np.real(sigma[h_idx, h_idx] - 0.5))
    epsilon = float(np.real(sigma[hc_idx, h_idx]))
    d_h = d[h_idx]
    d_hc = d[hc_idx]
    return delta, epsilon, d_h, d_hc, 1.0 + delta


def compute_p_herald_formula(
    sigma: np.ndarray,
    d: np.ndarray,
    *,
    h_idx: int = _H_IDX,
    hc_idx: int = _HC_IDX,
) -> tuple[float, float, float]:
    """P_herald = δ + |d_h|²  (PDF closed form, single herald mode)."""
    delta, _, d_h, d_hc, _ = extract_herald_params(sigma, d, h_idx=h_idx, hc_idx=hc_idx)
    d_h_sq = float(np.real(d_h * np.conj(d_hc)))
    return delta + d_h_sq, delta, d_h_sq


def compute_p_herald_integral(
    sigma: np.ndarray,
    d: np.ndarray,
    *,
    h_idx: int = _H_IDX,
    hc_idx: int = _HC_IDX,
) -> float:
    """P_herald from Gaussian integration (PDF trace-out + ∂²/∂ζ∂ζ*)."""
    sigma = np.asarray(sigma, dtype=np.complex128)
    d = np.asarray(d, dtype=np.complex128).reshape(-1)

    j = 1
    sym = _phase_symbols(j)
    u_total = _build_u_total(sp.Matrix(sigma.tolist()), sp.Matrix(d.tolist()), sym)
    s0 = {sym["s"][0]: 0, sym["sc"][0]: 0}
    u_traced = sp.expand(u_total.subs(s0))
    herald = _extract_herald_block(u_traced, sym)
    integrated = _integrate_heralds(herald, j)
    z0 = {sym["z"][0]: 0, sym["zc"][0]: 0}
    u_at_zero = float(integrated["u_eff"].subs(z0).evalf())

    delta, epsilon, d_h, d_hc, M = extract_herald_params(
        sigma, d, h_idx=h_idx, hc_idx=hc_idx
    )
    den = M**2 - epsilon**2
    z, zc = sym["z"][0], sym["zc"][0]
    f_zeta = (1 - M / den) * zc + (M * d_hc + epsilon * d_h - epsilon * z) / den
    f_zetac = (1 - M / den) * z + (M * d_h + epsilon * d_hc - epsilon * zc) / den
    f_cross = 1 - M / den
    f_z0 = complex(f_zeta.subs(z0).evalf())
    f_zc0 = complex(f_zetac.subs(z0).evalf())
    deriv = float(np.real(f_z0 * f_zc0 + f_cross))
    const = float(
        np.real(
            -(M * d_h * np.conj(d_hc) + 0.5 * epsilon * (d_h**2 + d_hc**2)) / den
        )
    )
    p_herald = float(np.exp(const) * deriv)
    rel = abs(u_at_zero + p_herald) / max(abs(p_herald), 1e-30)
    if rel > 1e-4:
        raise RuntimeError(
            f"u_eff(0)={u_at_zero:.6e} inconsistent with P_herald={p_herald:.6e}"
        )
    return p_herald


def compute_p_herald(
    sigma: np.ndarray | None = None,
    d: np.ndarray | None = None,
) -> PHeraldResult:
    """P_herald for single-copy state (default: PDF primary state)."""
    if sigma is None or d is None:
        st = generate_sigma_D()
        sigma, d = st.sigma, st.d
    p_formula, delta, d_h_sq = compute_p_herald_formula(sigma, d)
    p_integral = compute_p_herald_integral(sigma, d)
    _, epsilon, _, _, _ = extract_herald_params(sigma, d)
    return PHeraldResult(
        p_integral=p_integral,
        p_formula=p_formula,
        delta=delta,
        d_h_sq=d_h_sq,
        epsilon=epsilon,
    )


@dataclass
class JointHeraldResult:
    """Numerator and P_herald together (same heralded calculation)."""

    j: int
    numerator: NumeratorResult
    p_herald: PHeraldResult | None
    ratio: float | None

    def __str__(self) -> str:
        lines = [
            f"=== Joint heralded result  j = {self.j} ===",
            "--- Numerator (conditioned ⟨n̂⟩ numerator) ---",
            str(self.numerator),
        ]
        if self.p_herald is not None:
            lines += ["--- P_herald (denominator) ---", str(self.p_herald)]
            if self.ratio is not None:
                lines.append(f"Numerator / P_herald = {self.ratio:.6f}")
        return "\n".join(lines)


def compute_joint(
    j: int = 1,
    *,
    sigma_single: np.ndarray | None = None,
    d_single: np.ndarray | None = None,
    apply_hadamard: bool | None = None,
    include_p_herald: bool = True,
) -> JointHeraldResult:
    """
    Compute numerator and (for j=1 single copy) P_herald in one call.

    P_herald closed form applies to single-copy 4×4 Σ only; for j>1 use
    probability_zero_check.compute_p_herald_hafnian instead.
    """
    if sigma_single is None or d_single is None:
        st = generate_sigma_D()
        sigma_single, d_single = st.sigma, st.d

    num = compute(
        j=j,
        sigma_single=sigma_single,
        d_single=d_single,
        apply_hadamard=apply_hadamard,
    )

    p_herald = None
    ratio = None
    if include_p_herald and j == 1:
        p_herald = compute_p_herald(sigma_single, d_single)
        if p_herald.p_formula > 0:
            ratio = num.numerator / p_herald.p_formula

    return JointHeraldResult(
        j=j, numerator=num, p_herald=p_herald, ratio=ratio
    )


# =============================================================================
# §8  Pr(n̄) helpers (used by probability_zero_check.py)
# =============================================================================

def build_sigma_D_for_pr(
    j: int,
    cp: np.ndarray | None = None,
    apply_hadamard: bool | None = None,
) -> tuple[np.ndarray, np.ndarray, dict[int, int], StateMatrices]:
    """
    Build (Σ, D) for photon-number probability Pr(n̄).

    j = 1  →  single Fiurášek copy, Σ shape (4, 4).
    j > 1  →  j copies stacked, Σ shape (4j, 4j); optional Hadamard on signals.
    """
    if j < 1:
        raise ValueError("j must be >= 1")

    st = generate_sigma_D(cp=cp)
    if apply_hadamard is None:
        apply_hadamard = j > 1

    if j == 1:
        return st.sigma, st.d, {0: 2}, st

    if apply_hadamard and (j & (j - 1)):
        raise ValueError(f"j={j}: Hadamard requires j ∈ {{1, 2, 4, 8, …}}")

    from generate_cv_and_dv import (
        generate_cv_and_dv,
        generate_u_cv_and_dv_udag,
        rearrange_cv_and_dv,
    )

    cvs, dvs = {}, {}
    for k in range(j):
        cv, dv = generate_cv_and_dv(st.alpha, K=1, M=2, N=2, single_mode=True)
        cvs[k], dvs[k] = cv, dv

    sigma, d = rearrange_cv_and_dv(cvs, dvs, j, size=2 * j)
    mk_dict = {k: 2 for k in range(j)}

    if apply_hadamard:
        h = hadamard(j).astype(np.complex128) / np.sqrt(j) if j > 1 else np.eye(1)
        sigma, d = generate_u_cv_and_dv_udag(sigma, d, mk_dict, h)

    return sigma, d, mk_dict, st


def compute_p_herald_hafnian(
    sigma: np.ndarray,
    d: np.ndarray,
    mk_dict: dict[int, int],
    *,
    cutoff: int,
) -> tuple[float, tuple[int, ...]]:
    """P(1 photon in every herald mode) via hafnian — any j."""
    from generate_cv_and_dv import delete_cv, delete_vec
    from hafnian_batched_statistics import probability, slice_probabilities

    spatial = len(d) // 2
    del_array = np.zeros(spatial, dtype=int)
    herald: list[int] = []
    offset = 0
    for key in sorted(mk_dict):
        m = mk_dict[key]
        del_array[offset] = 1
        herald.extend([1] * (m - 1))
        offset += m

    deletion = np.concatenate((del_array, del_array))
    rcv = delete_cv(sigma, deletion)
    rdv = delete_vec(d, deletion)
    tensor = probability(rcv, rdv, cutoff=cutoff)
    p = float(np.real(slice_probabilities(tensor, tuple(herald))))
    return p, tuple(herald)
