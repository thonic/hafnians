"""
Shared Haar-ensemble helpers for codes_LR_Haar.

Cluster jobs write sparse P(n̄) only. Correlations / ensemble stats are offline.

Ensemble size, batch size, and base seed are read ONLY from each job's config.py:
    N_ENSEMBLE, ENSEMBLE_BATCH_SIZE, HAAR_BASE_SEED
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


def zero_save_tol(config) -> float:
    return float(config.ZERO_SAVE_TOL)


def haar_base_seed(config) -> int:
    """Single source of truth: config.HAAR_BASE_SEED."""
    return int(config.HAAR_BASE_SEED)


def n_ensemble(config) -> int:
    """Single source of truth: config.N_ENSEMBLE."""
    return int(config.N_ENSEMBLE)


def ensemble_batch_size(config) -> int:
    """Single source of truth: config.ENSEMBLE_BATCH_SIZE."""
    return int(config.ENSEMBLE_BATCH_SIZE)


def n_batches(config) -> int:
    """Number of PBS array tasks: ceil(N_ENSEMBLE / ENSEMBLE_BATCH_SIZE)."""
    n = n_ensemble(config)
    b = ensemble_batch_size(config)
    if b <= 0:
        raise ValueError(f"ENSEMBLE_BATCH_SIZE must be positive, got {b}")
    return max(1, int(math.ceil(n / b)))


def realization_seed(base_seed: int, realization_index: int) -> int:
    """Deterministic seed: HAAR_BASE_SEED + realization_index."""
    return int(base_seed) + int(realization_index)


def ensemble_root(here: Path, config) -> Path:
    """output/ensemble/R{N}_base{seed}/ — N and seed from config only."""
    base = haar_base_seed(config)
    n = n_ensemble(config)
    return here / config.OUTPUT_DIRNAME / "ensemble" / f"R{n}_base{base}"


def realization_json_path(
    ens_dir: Path,
    realization_index: int,
    *,
    geometry: str | None = None,
) -> Path:
    stem = f"seed_{int(realization_index):06d}"
    if geometry is not None:
        stem = f"{stem}_geom{geometry}"
    return ens_dir / f"{stem}.json"


def sparse_probability_list(
    probs: dict[tuple[int, ...], float],
    *,
    tol: float,
) -> list[list]:
    """
    Convert {nbar: P} → [[list(nbar), P], ...] for P >= tol, sorted by -P then nbar.
    """
    items = [(tuple(int(x) for x in k), float(v)) for k, v in probs.items() if float(v) >= tol]
    items.sort(key=lambda kv: (-kv[1], kv[0]))
    return [[list(nbar), p] for nbar, p in items]


def sum_sparse(probabilities: list) -> float:
    return float(sum(float(p) for _, p in probabilities))


def is_valid_realization(
    path: Path,
    *,
    expected_seed: int | None = None,
    expected_geometry: str | None = None,
    min_sum_p: float = 1.0 - 1e-6,
) -> bool:
    """True if JSON exists, parses, and looks like a completed sparse P(n̄) dump."""
    if not path.is_file():
        return False
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    if data.get("interferometer") != "haar":
        return False
    if "probabilities" not in data or not isinstance(data["probabilities"], list):
        return False
    if expected_seed is not None and int(data.get("seed", -1)) != int(expected_seed):
        return False
    if expected_geometry is not None and data.get("geometry") != expected_geometry:
        return False
    sum_p = data.get("sum_P")
    if sum_p is None:
        sum_p = sum_sparse(data["probabilities"])
    if float(sum_p) < min_sum_p:
        return False
    return True


def write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    """Write via temp file then rename for crash safety."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)


def batch_realization_range(
    *,
    batch_index: int,
    batch_size: int,
    n_total: int,
) -> range:
    """Realizations [batch_index * batch_size, min(...)] for PBS array index."""
    start = int(batch_index) * int(batch_size)
    end = min(start + int(batch_size), int(n_total))
    if start >= n_total:
        return range(0, 0)
    return range(start, end)


def write_manifest(ens_dir: Path, payload: dict[str, Any]) -> None:
    write_json_atomic(ens_dir / "manifest.json", payload)
