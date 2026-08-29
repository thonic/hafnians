#!/usr/bin/env python3
"""
Build codes_LR_Haar cluster folders from Walrus CLUSTER_* templates.

Folder names: CLUSTER_J2_even_cat_LR_Haar (Walrus + _LR + _Haar).
Scripts: run_parallel_J2_even_cat_LR_Haar.py, submit_J2_even_cat_LR_Haar.slurm, etc.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LR_ROOT = Path(__file__).resolve().parent

J2_HELPERS = (
    "build_B_mu_scale.py",
    "generate_cv_and_dv.py",
    "generate_displacements.py",
    "alpha_fiurasek.py",
    "hafnian_batched_statistics.py",
)

SPECS: list[dict] = [
    {
        "walrus_dir": "CLUSTER_J2_superposition",
        "lr_dir": "CLUSTER_J2_superposition_LR_Haar",
        "run_script": "run_parallel_J2_superposition_LR_Haar.py",
        "submit": "submit_J2_superposition_LR_Haar.slurm",
        "pbs_name": "pnbar_J2_superposition_LR_Haar",
        "log": "logs/j2_superposition_LR_Haar.log",
        "j": 2,
        "result_json": "Pnbar_superposition_J2_parallel_LR.json",
        "result_txt": "Pnbar_superposition_J2_parallel_LR.txt",
    },
    {
        "walrus_dir": "CLUSTER_J2_even_cat",
        "lr_dir": "CLUSTER_J2_even_cat_LR_Haar",
        "run_script": "run_parallel_J2_even_cat_LR_Haar.py",
        "submit": "submit_J2_even_cat_LR_Haar.slurm",
        "pbs_name": "pnbar_J2_even_cat_LR_Haar",
        "log": "logs/j2_even_cat_LR_Haar.log",
        "j": 2,
        "result_json": "Pnbar_J2_even_cat_LR.json",
        "result_txt": "Pnbar_J2_even_cat_LR.txt",
    },
    {
        "walrus_dir": "CLUSTER_J2_odd_cat",
        "lr_dir": "CLUSTER_J2_odd_cat_LR_Haar",
        "run_script": "run_parallel_J2_odd_cat_LR_Haar.py",
        "submit": "submit_J2_odd_cat_LR_Haar.slurm",
        "pbs_name": "pnbar_J2_odd_cat_LR_Haar",
        "log": "logs/j2_odd_cat_LR_Haar.log",
        "j": 2,
        "result_json": "Pnbar_J2_odd_cat_LR.json",
        "result_txt": "Pnbar_J2_odd_cat_LR.txt",
    },
    {
        "walrus_dir": "CLUSTER_J4_even_cat",
        "lr_dir": "CLUSTER_J4_even_cat_LR_Haar",
        "run_script": "run_parallel_J4_even_cat_LR_Haar.py",
        "submit": "submit_J4_even_cat_LR_Haar.slurm",
        "pbs_name": "pnbar_J4_even_cat_LR_Haar",
        "log": "logs/j4_even_cat_LR_Haar.log",
        "j": 4,
    },
    {
        "walrus_dir": "CLUSTER_J4_odd_cat",
        "lr_dir": "CLUSTER_J4_odd_cat_LR_Haar",
        "run_script": "run_parallel_J4_odd_cat_LR_Haar.py",
        "submit": "submit_J4_odd_cat_LR_Haar.slurm",
        "pbs_name": "pnbar_J4_odd_cat_LR_Haar",
        "log": "logs/j4_odd_cat_LR_Haar.log",
        "j": 4,
    },
]

LR_RUNNER_HEADER = '''\
"""Low-rank cluster runner — generated layout matches Walrus CLUSTER_* (engine LR only)."""
'''

SUBMIT_TEMPLATE = """#!/bin/bash
#PBS -N {pbs_name}
#PBS -j oe
#PBS -o {log}
#PBS -l walltime={walltime}
#PBS -l select=1:ncpus=8:mem={mem}

# Submit: mkdir -p logs output && qsub {submit}

set -euo pipefail
cd "${{PBS_O_WORKDIR:-$(pwd)}}"
mkdir -p logs output

NCPUS="${{PBS_NCPUS:-${{PBS_NP:-8}}}}"
PYTHON="${{PYTHON:-/storage/brno12-cerit/home/deeptisharma/.conda/envs/guassian/bin/python}}"

echo "HOST=$(hostname)  JOB=${{PBS_JOBID:-local}}  CPUS=${{NCPUS}}  ENGINE=low_rank"
"${{PYTHON}}" {run_script} --workers "${{NCPUS}}"
"""


def _patch_j2_config(text: str, result_json: str, result_txt: str) -> str:
    text = re.sub(
        r'^RESULT_JSON\s*=\s*".*?"',
        f'RESULT_JSON = "{result_json}"',
        text,
        count=1,
        flags=re.M,
    )
    text = re.sub(
        r'^RESULT_TXT\s*=\s*".*?"',
        f'RESULT_TXT = "{result_txt}"',
        text,
        count=1,
        flags=re.M,
    )
    return text


def _patch_j4_config(text: str, state: str) -> str:
    text = re.sub(
        r'^SUMMARY_JSON\s*=\s*".*?"',
        f'SUMMARY_JSON = "Pnbar_J4_{state}_LR_all_geometries.json"',
        text,
        count=1,
        flags=re.M,
    )
    text = re.sub(
        r'^SUMMARY_TXT\s*=\s*".*?"',
        f'SUMMARY_TXT = "Pnbar_J4_{state}_LR_all_geometries.txt"',
        text,
        count=1,
        flags=re.M,
    )
    text = text.replace(
        f"Pnbar_J4_{state}_geom",
        f"Pnbar_J4_{state}_LR_geom",
    )
    return text


def sync_spec(spec: dict) -> None:
    src = ROOT / spec["walrus_dir"]
    dst = LR_ROOT / spec["lr_dir"]
    if dst.exists():
        shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=shutil.ignore_patterns("output", "logs", "__pycache__", "*.pyc"))
    # Remove Walrus-only runners if present
    for p in dst.glob("run_parallel*.py"):
        if "_LR" not in p.name:
            p.unlink()
    for p in dst.glob("submit*.slurm"):
        if "_LR" not in p.name:
            p.unlink()

    cfg_path = dst / "config.py"
    cfg = cfg_path.read_text(encoding="utf-8")
    if spec["j"] == 2:
        cfg = _patch_j2_config(cfg, spec["result_json"], spec["result_txt"])
    else:
        state = "even_cat" if "even" in spec["lr_dir"] else "odd_cat"
        cfg = _patch_j4_config(cfg, state)
    cfg_path.write_text(cfg, encoding="utf-8")

    j2_body = (LR_ROOT / "_lr_run_parallel_j2.py").read_text()
    j4_body = (LR_ROOT / "_lr_run_parallel_j4.py").read_text()
    if spec["j"] == 2:
        (dst / spec["run_script"]).write_text(j2_body)
    else:
        (dst / spec["run_script"]).write_text(j4_body)

    walltime = "01:00:00" if spec["j"] == 2 else "04:00:00"
    mem = "8GB" if spec["j"] == 2 else "16GB"
    submit = SUBMIT_TEMPLATE.format(
        pbs_name=spec["pbs_name"],
        log=spec["log"],
        walltime=walltime,
        mem=mem,
        submit=spec["submit"],
        run_script=spec["run_script"],
    )
    (dst / spec["submit"]).write_text(submit)
    (dst / "output").mkdir(exist_ok=True)
    print(f"synced {spec['lr_dir']}")


def remove_legacy() -> None:
    for name in ("J2_even_cat", "J2_odd_cat", "J2_superposition", "_template"):
        p = LR_ROOT / name
        if p.exists():
            shutil.rmtree(p)


def main() -> None:
    remove_legacy()
    lr_src = ROOT / "probability_nbar_J2" / "lowrank" / "fast_loop_hafnian.py"
    lr_dst = LR_ROOT / "lowrank" / "fast_loop_hafnian.py"
    if lr_src.exists():
        shutil.copy2(lr_src, lr_dst)
    for spec in SPECS:
        sync_spec(spec)
    print("Done. codes_LR_Haar mirrors Walrus names with _LR_Haar suffix.")


if __name__ == "__main__":
    main()
