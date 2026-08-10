"""
J=2 P(n̄) — parallel per-pattern low-rank loop hafnian.

Same layout as Walrus run_parallel_J2_*.py; only the hafnian engine differs.
Takagi factor G is built once; workers use loop_hafnian_lr_from_G.

Usage
-----
    python run_parallel_J2_even_cat_LR_Haar.py
    python run_parallel_J2_even_cat_LR_Haar.py --workers 8
    qsub submit_J2_even_cat_LR_Haar.slurm
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
LR_ROOT = HERE.parent
if str(LR_ROOT) not in sys.path:
    sys.path.insert(0, str(LR_ROOT))

from lowrank.lr_production_kernel import loop_hafnian_lr_from_G, takagi_factor_b_mat  # noqa: E402

ENGINE = "low_rank_loop_hafnian_per_pattern_parallel"

_G_G = None
_G_MU = None
_G_SCALE = None
_G_MK = None


def _init_worker(G_factor, mu, scale, MK_dict) -> None:
    global _G_G, _G_MU, _G_SCALE, _G_MK
    _G_G = G_factor
    _G_MU = mu
    _G_SCALE = scale
    _G_MK = MK_dict


def _amplitude_prob(nbar: tuple[int, ...]) -> tuple[tuple[int, ...], float]:
    reps = reps_for_nbar(nbar, _G_MK)
    raw = loop_hafnian_lr_from_G(_G_G, _G_MU, reps)
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
) -> dict:
    J = int(J if J is not None else config.J)
    cutoff = int(cutoff if cutoff is not None else config.CUTOFF)
    n_workers = resolve_workers(workers)
    patterns = list(itertools.product(*[range(cutoff) for _ in range(J)]))
    n_patterns = len(patterns)

    print("=" * 60)
    print(f"Cluster LR: J={J}  state={config.STATE}")
    print("=" * 60)
    print(f"  cutoff    : {cutoff}  →  {n_patterns} configurations")
    print(f"  Interferometer : {config.INTERFEROMETER}")
    print(f"  workers   : {n_workers}")
    print(f"  engine    : {ENGINE}")
    print()

    t0 = time.perf_counter()
    print("Building (B, μ, scale) + Takagi G once ...", flush=True)
    b_mat, mu, scale, MK_dict = build_B_mu_scale(
        config.CP, J=J, interferometer=config.INTERFEROMETER
    )
    G = takagi_factor_b_mat(b_mat)
    t_build = time.perf_counter() - t0
    print(
        f"  B shape = {b_mat.shape}, Takagi rank = {G.shape[1]}, "
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
        initargs=(G, mu, scale, MK_dict),
    ) as pool:
        for nbar, p in pool.map(_amplitude_prob, patterns, chunksize=chunksize):
            if p > 0.0:
                raw[nbar] = p

    t_haf = time.perf_counter() - t1
    total_raw = float(sum(raw.values()))
    probs = {k: float(v / total_raw) for k, v in raw.items()} if total_raw > 0.0 else {}
    active = {k: v for k, v in probs.items() if v >= config.ZERO_TOL}
    t_total = time.perf_counter() - t0

    print(f"  hafnian wall time = {t_haf:.3f} s", flush=True)
    print(f"  Σ raw             = {total_raw:.6e}")
    print(f"  nonzero (tol)     = {len(active)} / {n_patterns}")
    print(f"  total wall time   = {t_total:.3f} s")

    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "state": config.STATE,
        "J": J,
        "cutoff": cutoff,
        "n_patterns": n_patterns,
        "n_workers": n_workers,
        "interferometer": config.INTERFEROMETER,
        "haar_random_seed": config.HAAR_RANDOM_SEED,
        "engine": ENGINE,
        "B_shape": list(b_mat.shape),
        "takagi_rank": int(G.shape[1]),
        "timing_s": {
            "build_B_mu_takagi": t_build,
            "parallel_hafnians": t_haf,
            "total": t_total,
        },
        "total_raw_before_renorm": total_raw,
        "sum_P": float(sum(active.values())),
        "probabilities": {str(k): float(v) for k, v in sorted(active.items())},
        "zero_tol": config.ZERO_TOL,
    }

    out_dir = HERE / config.OUTPUT_DIRNAME
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / config.RESULT_JSON
    txt_path = out_dir / config.RESULT_TXT
    json_path.write_text(json.dumps(result, indent=2), encoding="utf-8")

    lines = [
        f"J={J} {config.STATE} — parallel per-pattern low-rank",
        f"workers={n_workers}  takagi_rank={G.shape[1]}",
        f"build={t_build:.3f}s  hafnian={t_haf:.3f}s  total={t_total:.3f}s",
        f"sum_P={result['sum_P']:.8f}",
        "",
        "nbar  P",
    ]
    for nbar, p in sorted(active.items(), key=lambda kv: -kv[1]):
        lines.append(f"{nbar}  {p:.8e}")
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"\nWrote {json_path}")
    print(f"Wrote {txt_path}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Parallel J=2 low-rank P(n̄)")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--cutoff", type=int, default=None)
    parser.add_argument("--J", type=int, default=None)
    args = parser.parse_args()
    run(workers=args.workers, cutoff=args.cutoff, J=args.J)


if __name__ == "__main__":
    main()
