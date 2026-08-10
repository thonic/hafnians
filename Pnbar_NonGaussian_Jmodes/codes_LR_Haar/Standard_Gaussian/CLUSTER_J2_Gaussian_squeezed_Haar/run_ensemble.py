#!/usr/bin/env python3
"""
Haar ensemble batch driver (J=2).

Option B: one PBS array task runs a contiguous batch of realizations.
Each realization keeps the existing pattern-level ProcessPoolExecutor.
No nested multiprocessing. No correlation computation.

Usage
-----
    ./submit_ensemble.sh                          # PBS: array size from config.py
    python run_ensemble.py --batch-index 0        # one batch
    python run_ensemble.py --workers 8            # all N_ENSEMBLE locally
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
# Parent may be codes_LR_Haar or Standard_Gaussian/
for p in (HERE.parent, HERE.parent.parent):
    if (p / "ensemble_common.py").is_file() and str(p) not in sys.path:
        sys.path.insert(0, str(p))

import ensemble_common as eco  # noqa: E402
import config  # noqa: E402


def _load_runner():
    """Import run_parallel_*.py under its real module name (needed for ProcessPool pickling)."""
    matches = sorted(HERE.glob("run_parallel_*.py"))
    if len(matches) != 1:
        raise RuntimeError(f"expected one run_parallel_*.py in {HERE}, found {matches}")
    if str(HERE) not in sys.path:
        sys.path.insert(0, str(HERE))
    return importlib.import_module(matches[0].stem)


def main() -> None:
    parser = argparse.ArgumentParser(description="Haar ensemble batch (J=2)")
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--cutoff", type=int, default=None)
    parser.add_argument("--batch-index", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--start", type=int, default=None, help="inclusive realization index")
    parser.add_argument("--end", type=int, default=None, help="exclusive realization index")
    args = parser.parse_args()

    n_total = eco.n_ensemble(config)
    base = eco.haar_base_seed(config)
    batch_size = int(args.batch_size if args.batch_size is not None else eco.ensemble_batch_size(config))

    env_idx = os.environ.get("PBS_ARRAY_INDEX") or os.environ.get("SLURM_ARRAY_TASK_ID")
    if args.start is not None or args.end is not None:
        start = int(args.start if args.start is not None else 0)
        end = int(args.end if args.end is not None else n_total)
        indices = range(max(0, start), min(end, n_total))
        batch_index = None
    elif args.batch_index is not None or env_idx is not None:
        batch_index = int(args.batch_index if args.batch_index is not None else env_idx)
        indices = eco.batch_realization_range(
            batch_index=batch_index, batch_size=batch_size, n_total=n_total
        )
    else:
        # Local default: all realizations from config.N_ENSEMBLE
        indices = range(0, n_total)
        batch_index = None

    ens_dir = eco.ensemble_root(HERE, config)
    ens_dir.mkdir(parents=True, exist_ok=True)
    eco.write_manifest(
        ens_dir,
        {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "state": config.STATE,
            "J": config.J,
            "cutoff": config.CUTOFF,
            "interferometer": "haar",
            "haar_base_seed": base,
            "n_ensemble": n_total,
            "ensemble_batch_size": batch_size,
            "zero_save_tol": eco.zero_save_tol(config),
            "note": "Sparse P(n̄) only; correlations computed offline.",
        },
    )

    runner = _load_runner()
    print(
        f"Ensemble batch: realizations {list(indices)[:3]}..."
        f"{list(indices)[-1] if indices else 'empty'} "
        f"(batch_index={batch_index}, n={len(indices)})"
    )

    done = skipped = ran = 0
    for idx in indices:
        seed = eco.realization_seed(base, idx)
        out = eco.realization_json_path(ens_dir, idx)
        if eco.is_valid_realization(out, expected_seed=seed):
            print(f"[{idx:06d}] SKIP seed={seed}")
            skipped += 1
            done += 1
            continue
        print(f"[{idx:06d}] RUN  seed={seed} -> {out.name}")
        runner.run(
            workers=args.workers,
            cutoff=args.cutoff,
            haar_seed=seed,
            out_json=out,
            skip_existing=True,
        )
        ran += 1
        done += 1

    print(f"Batch finished: done={done} ran={ran} skipped={skipped}")


if __name__ == "__main__":
    main()
