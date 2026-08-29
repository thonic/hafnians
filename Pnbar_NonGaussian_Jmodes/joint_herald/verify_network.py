"""
Verify j-copy network: permutation (S,H,S*,H*) and Hadamard unitary.

Run:
    python verify_network.py
    python verify_network.py --j 4
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

from code_joint_herald_fiurasek import (
    build_network_sigma_D,
    describe_basis_order,
    generate_sigma_D,
    inspect_network,
    is_power_of_two,
    permutation_to_SHSH,
    valid_hadamard_copies,
)


def _print_report(info: dict) -> bool:
    j = info["j"]
    ok = all(info["checks_passed"].values()) if info["checks_passed"] else True

    print("=" * 72)
    print(f"j = {j} copies   Hadamard = {info['apply_hadamard']}")
    print("=" * 72)
    print(f"Valid j with Hadamard: {info['valid_hadamard_j']}")
    print(f"Basis order v = ({', '.join(info['basis_order'])})")
    print(f"Permutation (old index per new slot): {info['permutation_indices']}")
    print(f"Block slices in Σ, D (length {4 * j}):")
    for name, (a, b) in info["block_slices"].items():
        print(f"  {name:3s} → indices [{a}:{b})")
    print(f"\nΣ shape after network: {info['sigma_shape']}")
    print(f"D after permute (+ Hadamard if on):\n  {np.round(info['D_order'], 8)}")
    print(f"\nη (single-copy σ[0,3]) = {info['eta']:.6e}")
    print(f"ε (single-copy σ[0,2]) = {info['epsilon']:.6e}")
    print("\nΣ_{S,H*} block (signal–herald* correlations):")
    print("  BEFORE Hadamard (expect η·I_j):")
    print(info["SH_star_before"])
    if info["apply_hadamard"]:
        print("  AFTER Hadamard (expect η·U_S):")
        print(info["SH_star_after"])

    print("\nChecks:")
    for name, passed in info["checks_passed"].items():
        print(f"  [{('PASS' if passed else 'FAIL')}] {name}")
    print("=" * 72)
    print("OVERALL:", "OK" if ok else "FAILED")
    return ok


def verify_all(*, j_list: list[int] | None = None) -> bool:
    st = generate_sigma_D()
    if j_list is None:
        j_list = [1, 2, 4]

    all_ok = True
    for j in j_list:
        if j == 1:
            info = inspect_network(st.sigma, st.d, j, apply_hadamard=False)
        else:
            if not is_power_of_two(j):
                print(f"\nSkipping j={j} (not a power of 2 for Hadamard)")
                continue
            info = inspect_network(st.sigma, st.d, j, apply_hadamard=True)
        all_ok &= _print_report(info)

    # invalid j should raise
    print("\nExpect error for j=3 with Hadamard:")
    try:
        build_network_sigma_D(st.sigma, st.d, 3, apply_hadamard=True)
        print("  FAIL: no error raised")
        all_ok = False
    except ValueError as e:
        print(f"  PASS: {e}")

    return all_ok


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Verify permutation + Hadamard network")
    parser.add_argument(
        "--j",
        type=int,
        nargs="*",
        default=None,
        help="copy counts to test (default: 1 2 4)",
    )
    args = parser.parse_args()
    ok = verify_all(j_list=args.j)
    sys.exit(0 if ok else 1)
