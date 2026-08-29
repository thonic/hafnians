"""
J=4 P(n̄) — Kerr-squeezed κ=1.0; parallel low-rank loop hafnian, geometries A → B → C.

Mirrors Walrus run_parallel_J4_*.py; Takagi G once per geometry.
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

import ensemble_common as eco  # noqa: E402
from lowrank.lr_production_kernel import loop_hafnian_lr_from_G, takagi_factor_b_mat  # noqa: E402

ENGINE = "low_rank_loop_hafnian_per_pattern_parallel"

_G_G = None
_G_MU = None
_G_SCALE = None
_G_MK = None


def _input_nmax() -> int:
    if hasattr(config, "INPUT_NMAX"):
        return int(config.INPUT_NMAX)
    return int(config.ODD_CAT_NMAX)


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


def run_geometry(
    geometry: str,
    *,
    workers: int | None = None,
    cutoff: int | None = None,
    J: int | None = None,
    haar_seed: int | None = None,
    out_json: Path | None = None,
    skip_existing: bool = False,
) -> dict | None:
    if geometry not in config.GEOMETRIES:
        raise ValueError(f"unknown geometry {geometry!r}")

    J = int(J if J is not None else config.J)
    cutoff = int(cutoff if cutoff is not None else config.CUTOFF)
    n_workers = resolve_workers(workers)
    active = config.GEOMETRIES[geometry]
    layout = config.geometry_layout(active)
    states = config.mode_position_states(active)
    patterns = list(itertools.product(*[range(cutoff) for _ in range(J)]))
    n_patterns = len(patterns)

    seed = int(
        haar_seed if haar_seed is not None else eco.haar_base_seed(config)
    )
    save_tol = eco.zero_save_tol(config)

    if skip_existing and out_json is not None and eco.is_valid_realization(
        out_json, expected_seed=seed, expected_geometry=geometry
    ):
        print(f"SKIP existing valid realization: {out_json}")
        return None

    print("=" * 60)
    print(f"Cluster LR: J={J}  state={config.STATE}  geometry={geometry}  {layout}")
    print("=" * 60)
    print(f"  engine        : {ENGINE}")
    print(f"  cutoff        : {cutoff}  →  {n_patterns} configurations")
    print(f"  workers       : {n_workers}")
    print(f"  Haar seed     : {seed}")
    print()

    t0 = time.perf_counter()
    b_mat, mu, scale, MK_dict = build_B_mu_scale(
        states, J=J, interferometer=config.INTERFEROMETER, haar_seed=seed
    )
    G = takagi_factor_b_mat(b_mat)
    t_build = time.perf_counter() - t0
    print(
        f"  B shape = {b_mat.shape}, Takagi rank = {G.shape[1]}, "
        f"build = {t_build:.3f} s",
        flush=True,
    )

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
    sparse = eco.sparse_probability_list(probs, tol=save_tol)
    sum_p = eco.sum_sparse(sparse)
    t_total = time.perf_counter() - t0

    print(
        f"  hafnian wall time = {t_haf:.3f} s  "
        f"sum_P={sum_p:.6f}  sparse={len(sparse)}/{n_patterns}"
    )

    result = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "state": config.STATE,
        "kappa": float(config.KERR_SQUEEZED_KAPPA),
        "squeeze_r": float(config.SQUEEZE_R),
        "J": J,
        "geometry": geometry,
        "layout": layout,
        "active_copies": list(active),
        "cutoff": cutoff,
        "alpha": config.CAT_ALPHA,
        "input_nmax": _input_nmax(),
        "n_patterns": n_patterns,
        "n_workers": n_workers,
        "interferometer": config.INTERFEROMETER,
        "seed": seed,
        "haar_base_seed": eco.haar_base_seed(config),
        "engine": ENGINE,
        "B_shape": list(b_mat.shape),
        "takagi_rank": int(G.shape[1]),
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
        print(f"Wrote {json_path}")
        return result

    out_dir = HERE / config.OUTPUT_DIRNAME
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / config.result_json_name(geometry)
    txt_path = out_dir / config.result_txt_name(geometry)
    eco.write_json_atomic(json_path, result)

    lines = [
        f"J={J} {config.STATE} Haar geometry={geometry} {layout}",
        f"seed={seed}  takagi_rank={G.shape[1]}  workers={n_workers}",
        f"build={t_build:.3f}s  hafnian={t_haf:.3f}s  sum_P={sum_p:.8f}",
        "",
        "nbar  P",
    ]
    for nbar, p in sparse:
        lines.append(f"{tuple(nbar)}  {p:.8e}")
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {json_path}")
    return result


def run_all(
    *,
    workers: int | None = None,
    cutoff: int | None = None,
    geometries: tuple[str, ...] | None = None,
    haar_seed: int | None = None,
) -> dict[str, dict]:
    geom_list = geometries or config.GEOMETRY_ORDER
    results = {
        g: run_geometry(g, workers=workers, cutoff=cutoff, haar_seed=haar_seed)
        for g in geom_list
    }

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "state": config.STATE,
        "kappa": float(config.KERR_SQUEEZED_KAPPA),
        "squeeze_r": float(config.SQUEEZE_R),
        "J": config.J,
        "engine": ENGINE,
        "cutoff": int(cutoff if cutoff is not None else config.CUTOFF),
        "alpha": config.CAT_ALPHA,
        "input_nmax": _input_nmax(),
        "interferometer": config.INTERFEROMETER,
        "haar_base_seed": eco.haar_base_seed(config),
        "geometry_order": list(geom_list),
        "geometries": results,
    }
    out_dir = HERE / config.OUTPUT_DIRNAME
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_json = out_dir / config.SUMMARY_JSON
    summary_txt = out_dir / config.SUMMARY_TXT
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    txt_lines = [
        f"J={config.J} {config.STATE} LR — all geometries",
        f"order: {' → '.join(geom_list)}",
        "",
    ]
    for g in geom_list:
        r = results[g]
        txt_lines.append(f"=== geometry {g}  {r['layout']}  sum_P={r['sum_P']:.8f} ===")
        for nbar, p in r["probabilities"]:
            txt_lines.append(f"{tuple(nbar)}  {float(p):.8e}")
        txt_lines.append("")
    summary_txt.write_text("\n".join(txt_lines) + "\n", encoding="utf-8")
    print(f"Wrote summary {summary_json}")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description="Parallel J=4 low-rank P(n̄)")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--cutoff", type=int, default=None)
    parser.add_argument(
        "--geometry",
        type=str,
        default=None,
        choices=tuple(config.GEOMETRIES),
    )
    parser.add_argument("--haar-seed", type=int, default=None)
    args = parser.parse_args()
    if args.geometry:
        run_geometry(
            args.geometry,
            workers=args.workers,
            cutoff=args.cutoff,
            haar_seed=args.haar_seed,
        )
    else:
        run_all(
            workers=args.workers,
            cutoff=args.cutoff,
            haar_seed=args.haar_seed,
        )


if __name__ == "__main__":
    main()
