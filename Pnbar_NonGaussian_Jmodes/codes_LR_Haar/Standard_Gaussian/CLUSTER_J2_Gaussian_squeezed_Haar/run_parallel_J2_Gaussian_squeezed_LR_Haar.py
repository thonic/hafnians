"""
J=2 P(n̄) — parallel per-pattern The Walrus low-rank hafnian (ordinary, not loop).

Same Fiurášek layout as the even/odd-cat and Kerr LR jobs; only the hafnian
engine differs: Takagi factor G of B, expand by herald reps, then
``low_rank_hafnian(G_exp)`` for A = G G^T (displacements unused / μ = 0).

Usage
-----
    python run_parallel_J2_Gaussian_squeezed_LR_Haar.py
    python run_parallel_J2_Gaussian_squeezed_LR_Haar.py --workers 8
    qsub submit_J2_Gaussian_squeezed_LR_Haar.slurm
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from math import factorial, prod
from pathlib import Path

import numpy as np

import config
from build_B_mu_scale import build_B_mu_scale
from hafnian_batched_statistics import reps_for_nbar

HERE = Path(__file__).resolve().parent
SG_ROOT = HERE.parent
LR_ROOT = SG_ROOT.parent
for p in (str(SG_ROOT), str(LR_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

import ensemble_common as eco  # noqa: E402
from lowrank.cluster_reps import expand_rows  # noqa: E402
from lowrank.lr_production_kernel import takagi_factor_b_mat  # noqa: E402
from walrus_low_rank_hafnian import (  # noqa: E402
    low_rank_hafnian,
    verify_against_walrus_hafnian,
)

ENGINE = "ordinary_low_rank_hafnian_AeqGGT_per_pattern_parallel"

_G_G = None
_G_SCALE = None
_G_MK = None


def _init_worker(G_factor, scale, MK_dict) -> None:
    global _G_G, _G_SCALE, _G_MK
    _G_G = G_factor
    _G_SCALE = scale
    _G_MK = MK_dict


def _amplitude_prob(nbar: tuple[int, ...]) -> tuple[tuple[int, ...], float]:
    reps = reps_for_nbar(nbar, _G_MK)
    # Ordinary hafnian of A = G G^T: expand Takagi rows; do not use μ (loop).
    zeros = np.zeros(_G_G.shape[0], dtype=complex)
    G_exp, _ = expand_rows(_G_G, zeros, reps)
    raw = low_rank_hafnian(G_exp)
    denom = np.sqrt(prod(factorial(int(r)) for r in np.asarray(reps).tolist()))
    amp = (raw / denom) * _G_SCALE
    p = float((amp * np.conj(amp)).real)
    return nbar, p


def resolve_workers(requested: int | None) -> int:
    if requested is not None and requested > 0:
        return int(requested)
    if config.N_WORKERS is not None and config.N_WORKERS > 0:
        return int(config.N_WORKERS)
    pbs = os.environ.get("PBS_NCPUS") or os.environ.get("PBS_NP")
    if pbs:
        return max(1, int(pbs))
    slurm = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm:
        return max(1, int(slurm))
    return max(1, os.cpu_count() or 1)


def run(
    *,
    workers: int | None = None,
    cutoff: int | None = None,
    J: int | None = None,
    haar_seed: int | None = None,
    out_json: Path | None = None,
    skip_existing: bool = False,
) -> dict | None:
    J = int(J if J is not None else config.J)
    cutoff = int(cutoff if cutoff is not None else config.CUTOFF)
    n_workers = resolve_workers(workers)
    patterns = list(itertools.product(*[range(cutoff) for _ in range(J)]))
    n_patterns = len(patterns)

    seed = int(
        haar_seed if haar_seed is not None else eco.haar_base_seed(config)
    )
    save_tol = eco.zero_save_tol(config)

    if skip_existing and out_json is not None and eco.is_valid_realization(
        out_json, expected_seed=seed
    ):
        print(f"SKIP existing valid realization: {out_json}")
        return None

    print("=" * 60)
    print(f"Cluster LR: J={J}  state={config.STATE}")
    print("=" * 60)
    print(f"  cutoff    : {cutoff}  →  {n_patterns} configurations")
    print(f"  Interferometer : {config.INTERFEROMETER}")
    print(f"  Haar seed : {seed}")
    print(f"  workers   : {n_workers}")
    print(f"  engine    : {ENGINE}")
    print(f"  squeeze_r : {config.SQUEEZE_R:.6f}  (sinh(r)=1)")
    print()

    verify_against_walrus_hafnian()

    t0 = time.perf_counter()
    print("Building (B, μ, scale) + Takagi G once ...", flush=True)
    b_mat, mu, scale, MK_dict = build_B_mu_scale(
        config.CP, J=J, interferometer=config.INTERFEROMETER, haar_seed=seed
    )
    G = takagi_factor_b_mat(b_mat)
    t_build = time.perf_counter() - t0
    print(
        f"  B shape = {b_mat.shape}, Takagi rank = {G.shape[1]}, "
        f"||μ|| = {np.linalg.norm(mu):.3e} (unused by ordinary hafnian), "
        f"build time = {t_build:.3f} s",
        flush=True,
    )

    print(f"Evaluating {n_patterns} patterns with {n_workers} workers ...", flush=True)
    t1 = time.perf_counter()
    raw: dict[tuple[int, ...], float] = {}
    chunksize = max(1, n_patterns // (n_workers * 4))

    with ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=_init_worker,
        initargs=(G, scale, MK_dict),
    ) as pool:
        for nbar, p in pool.map(_amplitude_prob, patterns, chunksize=chunksize):
            if p > 0.0:
                raw[nbar] = p

    t_haf = time.perf_counter() - t1
    total_raw = float(sum(raw.values()))
    probs = {k: float(v / total_raw) for k, v in raw.items()} if total_raw > 0.0 else {}
    sparse = eco.sparse_probability_list(probs, tol=save_tol)
    sum_p = eco.sum_sparse(sparse)
    t_total = time.perf_counter() - t0

    print(f"  hafnian wall time = {t_haf:.3f} s", flush=True)
    print(f"  Σ raw             = {total_raw:.6e}")
    print(f"  sparse entries    = {len(sparse)} / {n_patterns}  (tol={save_tol:g})")
    print(f"  sum_P (sparse)    = {sum_p:.16f}")
    print(f"  total wall time   = {t_total:.3f} s")

    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "state": config.STATE,
        "J": J,
        "cutoff": cutoff,
        "squeeze_r": float(config.SQUEEZE_R),
        "input_nmax": int(config.INPUT_NMAX),
        "n_patterns": n_patterns,
        "n_workers": n_workers,
        "interferometer": config.INTERFEROMETER,
        "seed": seed,
        "haar_base_seed": eco.haar_base_seed(config),
        "engine": ENGINE,
        "B_shape": list(b_mat.shape),
        "takagi_rank": int(G.shape[1]),
        "mu_norm_unused": float(np.linalg.norm(mu)),
        "timing_s": {
            "build_B_mu_takagi": t_build,
            "parallel_hafnians": t_haf,
            "total": t_total,
        },
        "total_raw_before_renorm": total_raw,
        "sum_P": sum_p,
        "probabilities": sparse,
        "zero_save_tol": save_tol,
    }

    if out_json is not None:
        json_path = Path(out_json)
        eco.write_json_atomic(json_path, result)
        print(f"\nWrote {json_path}")
        return result

    out_dir = HERE / config.OUTPUT_DIRNAME
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / config.RESULT_JSON
    txt_path = out_dir / config.RESULT_TXT
    eco.write_json_atomic(json_path, result)

    lines = [
        f"J={J} {config.STATE} — Haar sparse P(n̄)",
        f"seed={seed}  workers={n_workers}  takagi_rank={G.shape[1]}",
        f"build={t_build:.3f}s  hafnian={t_haf:.3f}s  total={t_total:.3f}s",
        f"sum_P={sum_p:.8f}  sparse={len(sparse)}/{n_patterns}",
        "",
        "nbar  P",
    ]
    for nbar, p in sparse:
        lines.append(f"{tuple(nbar)}  {p:.8e}")
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"\nWrote {json_path}")
    print(f"Wrote {txt_path}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parallel J=2 standard-Gaussian low-rank P(n̄)"
    )
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--cutoff", type=int, default=None)
    parser.add_argument("--J", type=int, default=None)
    parser.add_argument("--haar-seed", type=int, default=None)
    args = parser.parse_args()
    run(workers=args.workers, cutoff=args.cutoff, J=args.J, haar_seed=args.haar_seed)


if __name__ == "__main__":
    main()
