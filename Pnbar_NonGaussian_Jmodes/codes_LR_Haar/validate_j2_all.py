#!/usr/bin/env python3
"""J=2 LR vs Walrus for all CLUSTER_J2_*_LR folders."""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import sys
from math import factorial, prod
from pathlib import Path

import numpy as np
from thewalrus import loop_hafnian

ROOT = Path(__file__).resolve().parents[1]
PNBAR = ROOT / "probability_nbar_J2"
LR_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PNBAR))
sys.path.insert(0, str(LR_ROOT))

from lowrank.lr_production_kernel import loop_hafnian_lr, takagi_factor_b_mat  # noqa: E402

STATES = (
    (
        "CLUSTER_J2_superposition_LR_Haar",
        ROOT / "CLUSTER_J2_superposition" / "output" / "Pnbar_superposition_J2_parallel.json",
        "Pnbar_superposition_J2_parallel_LR.json",
    ),
    (
        "CLUSTER_J2_even_cat_LR_Haar",
        ROOT / "CLUSTER_J2_even_cat" / "output" / "Pnbar_J2_even_cat.json",
        "Pnbar_J2_even_cat_LR.json",
    ),
    (
        "CLUSTER_J2_odd_cat_LR_Haar",
        ROOT / "CLUSTER_J2_odd_cat" / "output" / "Pnbar_J2_odd_cat.json",
        "Pnbar_J2_odd_cat_LR.json",
    ),
)


def run_kernel_check(lr_folder: str) -> dict:
    here = LR_ROOT / lr_folder
    spec = importlib.util.spec_from_file_location("cfg", here / "config.py")
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    sys.path.insert(0, str(here))
    from build_B_mu_scale import build_B_mu_scale  # noqa: E402
    from hafnian_batched_statistics import reps_for_nbar  # noqa: E402

    patterns = list(itertools.product(*[range(cfg.CUTOFF) for _ in range(cfg.J)]))
    b_mat, mu, scale, MK = build_B_mu_scale(cfg.CP, J=cfg.J, interferometer=cfg.INTERFEROMETER)
    G = takagi_factor_b_mat(b_mat)
    max_da = max_dp = 0.0
    for nbar in patterns:
        reps = reps_for_nbar(nbar, MK)
        w = loop_hafnian(A=b_mat, D=mu, reps=reps, glynn=True)
        lr = loop_hafnian_lr(b_mat, mu, reps, G=G)
        max_da = max(max_da, abs(w - lr))
        denom = np.sqrt(prod(factorial(int(r)) for r in np.asarray(reps).tolist()))
        max_dp = max(max_dp, abs(abs((w / denom) * scale) ** 2 - abs((lr / denom) * scale) ** 2))
    return {
        "folder": lr_folder,
        "state": cfg.STATE,
        "patterns": len(patterns),
        "takagi_rank": int(G.shape[1]),
        "max_abs_loop_haf_diff": max_da,
        "max_abs_prob_diff": max_dp,
    }


def compare_json(lr_folder: str, ref_path: Path, lr_json_name: str) -> dict:
    lr_path = LR_ROOT / lr_folder / "output" / lr_json_name
    if not lr_path.exists():
        return {"folder": lr_folder, "error": f"missing {lr_path}"}
    if not ref_path.exists():
        return {"folder": lr_folder, "error": f"missing ref {ref_path}"}
    p_lr = json.loads(lr_path.read_text()).get("probabilities", {})
    p_ref = json.loads(ref_path.read_text()).get("probabilities", {})
    keys = set(p_lr) | set(p_ref)
    max_dp = max(abs(float(p_lr.get(k, 0)) - float(p_ref.get(k, 0))) for k in keys)
    return {"folder": lr_folder, "max_prob_diff": max_dp}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compare-json", action="store_true")
    args = parser.parse_args()
    print("=== Kernel: LR vs Walrus ===")
    for folder, _, _ in STATES:
        print(run_kernel_check(folder))
    if args.compare_json:
        print("\n=== JSON vs Walrus CLUSTER output ===")
        for folder, ref, lr_name in STATES:
            print(compare_json(folder, ref, lr_name))


if __name__ == "__main__":
    main()
