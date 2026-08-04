"""
Secondary-state verification: (|1⟩ + |2⟩) / √2  (n = 3 modes)

Checks that P_herald = δ + |d_h|² holds **per herald mode** and compares
marginal and joint herald probabilities against hafnian_batched.

Not part of the main API — run explicitly:
    python verify_secondary_state.py
"""

from __future__ import annotations

import sys

import numpy as np

from generate_cv_and_dv import delete_cv, delete_vec, generate_cv_and_dv
from generate_displacements import generate_alpha
from hafnian_batched_statistics import probability

from code_joint_herald_fiurasek import compute_p_herald_formula, extract_herald_params

# Target: (|1⟩ + |2⟩) / √2  →  cp = [0, 1, 1] / √2
SECONDARY_CP = np.array([0.0, 1.0, 1.0], dtype=np.complex128)
SECONDARY_CP /= np.linalg.norm(SECONDARY_CP)

# Mode layout for M = 3: v = (m0, m1, m2, m0*, m1*, m2*)
HERALD_MODES = [
    {"label": "Herald mode 1", "mode": 1, "h_idx": 1, "hc_idx": 4, "delete": [1, 0, 1]},
    {"label": "Herald mode 2", "mode": 2, "h_idx": 2, "hc_idx": 5, "delete": [1, 1, 0]},
]

# Reference values (hafnian_batched, cutoff = 2)
REF_MARGINAL = {
    1: 1.074170e-07,
    2: 2.915154e-08,
}
REF_JOINT_11 = 1.600000e-15


def build_secondary_sigma_D() -> tuple[np.ndarray, np.ndarray]:
    alpha = generate_alpha(SECONDARY_CP)
    sigma, d = generate_cv_and_dv(alpha, K=1, M=3, N=3, single_mode=True)
    return np.asarray(sigma, dtype=np.complex128), np.asarray(d, dtype=np.complex128).reshape(-1)


def marginal_herald_probability(
    sigma: np.ndarray,
    d: np.ndarray,
    delete_modes: list[int],
    *,
    cutoff: int = 2,
) -> float:
    """P( exactly 1 photon in the kept herald mode )."""
    del_array = np.asarray(delete_modes, dtype=int)
    deletion = np.concatenate((del_array, del_array))
    rcv = delete_cv(sigma, deletion)
    rdv = delete_vec(d, deletion)
    return float(np.real(probability(rcv, rdv, cutoff=cutoff)[1, 1]))


def joint_herald_probability(
    sigma: np.ndarray,
    d: np.ndarray,
    *,
    signal_mode: int = 0,
    cutoff: int = 2,
) -> float:
    """P( 1 photon in herald 1 AND 1 photon in herald 2 ) — delete signal only."""
    M = len(d) // 2
    del_array = np.zeros(M, dtype=int)
    del_array[signal_mode] = 1  # 1 = delete in delete_cv
    deletion = np.concatenate((del_array, del_array))
    rcv = delete_cv(sigma, deletion)
    rdv = delete_vec(d, deletion)
    prob = probability(rcv, rdv, cutoff=cutoff)
    arr = np.asarray(prob)
    if arr.ndim == 4:
        return float(np.real(arr[1, 1, 1, 1]))
    if arr.ndim == 2:
        return float(np.real(arr[1, 1]))
    raise ValueError(f"Unexpected probability shape {arr.shape}")


def verify(*, rtol: float = 1e-4) -> bool:
    sigma, d = build_secondary_sigma_D()
    ok = True

    print("=" * 72)
    print("Secondary state: (|1⟩ + |2⟩) / √2   (Fiurásek n = 3 modes)")
    print("=" * 72)
    print(f"cp = {np.round(SECONDARY_CP, 4)}")
    print(f"Σ shape: {sigma.shape}   D length: {d.size}")
    print()

    print("--- Per-herald blocks (δ + |d_h|² vs hafnian marginal) ---")
    marginals = []
    for hm in HERALD_MODES:
        p_form, delta, d_h_sq = compute_p_herald_formula(
            sigma, d, h_idx=hm["h_idx"], hc_idx=hm["hc_idx"]
        )
        p_haf = marginal_herald_probability(sigma, d, hm["delete"])
        ref = REF_MARGINAL[hm["mode"]]
        err = abs(p_form - p_haf) / max(abs(p_haf), 1e-30)
        pass_m = err < rtol
        ok &= pass_m
        marginals.append(p_haf)
        print(f"{hm['label']}  (indices {hm['h_idx']}, {hm['hc_idx']}):")
        print(f"  δ           = {delta:.6e}")
        print(f"  |d_h|²      = {d_h_sq:.6e}")
        print(f"  formula     = {p_form:.6e}")
        print(f"  hafnian     = {p_haf:.6e}  (ref {ref:.6e})")
        print(f"  rel. error  = {err:.2%}  {'PASS' if pass_m else 'FAIL'}")
        print()

    p_joint = joint_herald_probability(sigma, d)
    p_prod = marginals[0] * marginals[1]
    print("--- Joint herald event |1,1⟩ (both heralds click) ---")
    print(f"  hafnian P(|1,1⟩)_herald = {p_joint:.6e}  (ref {REF_JOINT_11:.6e})")
    print(f"  product of marginals    = {p_prod:.6e}")
    print("  Note: joint ≠ δ₁+δ₂ and joint ≈ P₁·P₂ only approximately")
    print("        (correlations modify the joint by O(10⁻¹⁵)).")
    err_j = abs(p_joint - REF_JOINT_11) / REF_JOINT_11
    pass_j = err_j < 0.05  # 5% on tiny joint probability
    ok &= pass_j
    print(f"  joint vs ref: {err_j:.2%}  {'PASS' if pass_j else 'FAIL'}")
    print()
    print("=" * 72)
    print("OVERALL:", "ALL CHECKS PASSED" if ok else "SOME CHECKS FAILED")
    print("=" * 72)
    return ok


if __name__ == "__main__":
    success = verify()
    sys.exit(0 if success else 1)
