"""
Fock-basis states for driver_code / probability_zero_check.

Pick a cp array below (same as writing cp = np.array([...]) in driver_code).
Each cp is UN-normalised; driver_code normalises it internally.

Usage in probability_zero_check.py
------------------------------------
    from state_definitions import CP_SUPERPOSITION
    cp = CP_SUPERPOSITION
"""

from __future__ import annotations

import numpy as np
from scipy.special import factorial

# =============================================================================
# Ready-to-use cp arrays  (pick ONE in probability_zero_check.py)
# =============================================================================

# Primary PDF state: (|0⟩ + |1⟩) / √2  — works with any K
CP_SUPERPOSITION = np.array([1.0, 1.0], dtype=np.complex128)

# Secondary: (|1⟩ + |2⟩) / √2  — use K = 1 only (6×6 Σ)
CP_SECONDARY = np.array([0.0, 1.0, 1.0], dtype=np.complex128)

# (|0⟩ + |2⟩) / √2  — use K = 1 only
CP_ZERO_TWO = np.array([1.0, 0.0, 1.0], dtype=np.complex128)


def cp_kerr_squeezed(kappa: float = 0.1, r: float | None = None, nmax: int = 10) -> np.ndarray:
    """Kerr squeezed: e^{-i κ n²}|S(r)⟩ (standard evolution convention)."""
    if r is None:
        r = float(np.arcsinh(0.5))
    cp = np.zeros(nmax, dtype=np.complex128)
    sech_r = 1.0 / np.cosh(r)
    for n in range(0, nmax, 2):
        phase = np.exp(-1j * kappa * n * n)
        amp = (
            np.sqrt(float(factorial(n)) * (np.tanh(r) ** n))
            / float(factorial(n // 2))
            / (2.0 ** (n / 2))
        )
        cp[n] = phase * amp
    return cp * np.sqrt(sech_r)


def cp_even_cat(alpha: complex = 0.7, nmax: int = 20) -> np.ndarray:
    """(|α⟩ + |-α⟩) / N"""
    alpha = complex(alpha)
    cp = np.zeros(nmax, dtype=np.complex128)
    damp = np.exp(-0.5 * abs(alpha) ** 2)
    for n in range(nmax):
        cn = damp * (alpha ** n) / np.sqrt(float(factorial(n)))
        cn_neg = damp * ((-alpha) ** n) / np.sqrt(float(factorial(n)))
        cp[n] = cn + cn_neg
    return cp


def cp_odd_cat(alpha: complex = 0.7, nmax: int = 4) -> np.ndarray:
    """(|α⟩ − |-α⟩) / N  (use even nmax; default 4)."""
    alpha = complex(alpha)
    cp = np.zeros(nmax, dtype=np.complex128)
    damp = np.exp(-0.5 * abs(alpha) ** 2)
    for n in range(nmax):
        cn = damp * (alpha ** n) / np.sqrt(float(factorial(n)))
        cn_neg = damp * ((-alpha) ** n) / np.sqrt(float(factorial(n)))
        cp[n] = cn - cn_neg
    return cp


def cp_kerr_coherent(alpha: complex = 0.7, kappa: float = 0.5, nmax: int = 7) -> np.ndarray:
    """Kerr coherent: e^{-i κ n²}|α⟩ truncated |0⟩…|nmax-1⟩ (standard evolution convention)."""
    alpha = complex(alpha)
    cp = np.zeros(nmax, dtype=np.complex128)
    damp = np.exp(-0.5 * abs(alpha) ** 2)
    for n in range(nmax):
        cp[n] = damp * (alpha ** n) / np.sqrt(float(factorial(n))) * np.exp(-1j * kappa * n * n)
    return cp


# Normalised versions (for code_joint_herald_fiurasek / generate_sigma_D)
def psi_superposition() -> np.ndarray:
    return CP_SUPERPOSITION / np.linalg.norm(CP_SUPERPOSITION)


def psi_secondary() -> np.ndarray:
    return CP_SECONDARY / np.linalg.norm(CP_SECONDARY)


def psi_kerr_squeezed(**kw) -> np.ndarray:
    cp = cp_kerr_squeezed(**kw)
    return cp / np.linalg.norm(cp)


def psi_even_cat(**kw) -> np.ndarray:
    cp = cp_even_cat(**kw)
    return cp / np.linalg.norm(cp)


def psi_odd_cat(**kw) -> np.ndarray:
    cp = cp_odd_cat(**kw)
    return cp / np.linalg.norm(cp)


def psi_kerr_coherent(**kw) -> np.ndarray:
    cp = cp_kerr_coherent(**kw)
    return cp / np.linalg.norm(cp)


def print_coefficients(cp: np.ndarray, *, label: str = "state") -> None:
    """Print n and c_n."""
    cp = np.asarray(cp, dtype=np.complex128).reshape(-1)
    cp_n = cp / np.linalg.norm(cp)
    print(f"--- {label}  (len={len(cp)}) ---")
    print(f"{'n':>4}  {'Re(c_n)':>12}  {'Im(c_n)':>12}  {'|c_n|²':>10}")
    for n in range(len(cp)):
        c = cp_n[n]
        print(f"{n:4d}  {np.real(c):12.6e}  {np.imag(c):12.6e}  {abs(c)**2:10.6e}")


if __name__ == "__main__":
    for label, cp in [
        ("CP_SUPERPOSITION", CP_SUPERPOSITION),
        ("CP_SECONDARY", CP_SECONDARY),
        ("cp_kerr_squeezed", cp_kerr_squeezed()),
        ("cp_even_cat", cp_even_cat()),
        ("cp_odd_cat", cp_odd_cat()),
        ("cp_kerr_coherent", cp_kerr_coherent()),
    ]:
        print_coefficients(cp, label=label)
        print()
