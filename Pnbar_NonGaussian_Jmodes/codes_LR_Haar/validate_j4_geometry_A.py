#!/usr/bin/env python3
"""LR vs Walrus on J=4 geometry A (1296 patterns) for even + odd cat LR folders."""

from __future__ import annotations

import importlib.util
import itertools
import sys
from math import factorial, prod
from pathlib import Path

import numpy as np
from thewalrus import loop_hafnian

ROOT = Path(__file__).resolve().parents[1]
LR_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "probability_nbar_J2"))
sys.path.insert(0, str(LR_ROOT))

from lowrank.lr_production_kernel import loop_hafnian_lr, takagi_factor_b_mat  # noqa: E402


def check(folder: str) -> dict:
    here = LR_ROOT / folder
    spec = importlib.util.spec_from_file_location("cfg", here / "config.py")
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    sys.path.insert(0, str(here))
    from build_B_mu_scale import build_B_mu_scale  # noqa: E402
    from hafnian_batched_statistics import reps_for_nbar  # noqa: E402

    g = "A"
    active = cfg.GEOMETRIES[g]
    states = cfg.mode_position_states(active)
    b_mat, mu, scale, MK = build_B_mu_scale(states, J=cfg.J, interferometer=cfg.INTERFEROMETER)
    G = takagi_factor_b_mat(b_mat)
    patterns = list(itertools.product(*[range(cfg.CUTOFF) for _ in range(cfg.J)]))
    max_da = max_dp = 0.0
    for nbar in patterns:
        reps = reps_for_nbar(nbar, MK)
        w = loop_hafnian(A=b_mat, D=mu, reps=reps, glynn=True)
        lr = loop_hafnian_lr(b_mat, mu, reps, G=G)
        max_da = max(max_da, abs(w - lr))
        denom = np.sqrt(prod(factorial(int(r)) for r in np.asarray(reps).tolist()))
        max_dp = max(max_dp, abs(abs((w / denom) * scale) ** 2 - abs((lr / denom) * scale) ** 2))
    return {
        "folder": folder,
        "patterns": len(patterns),
        "takagi_rank": int(G.shape[1]),
        "max_abs_loop_haf_diff": max_da,
        "max_abs_prob_diff": max_dp,
    }


def main() -> None:
    for folder in ("CLUSTER_J4_even_cat_LR_Haar", "CLUSTER_J4_odd_cat_LR_Haar"):
        print(check(folder))


if __name__ == "__main__":
    main()
