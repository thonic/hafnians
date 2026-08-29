# Low-rank cluster jobs — Haar interferometer (`codes_LR_Haar`)

> **Haar interferometer project:** Self-contained low-rank `P(n̄)` cluster jobs with a
> Haar-random passive unitary on signal modes (`INTERFEROMETER = "haar"` in each job's
> `config.py`). Set `INTERFEROMETER = "identity"` to skip mode mixing. The Fiurášek
> preparation, geometry construction, low-rank kernel, and hafnian evaluation are unchanged
> from the parent low-rank pipeline.

**Self-contained for Cerit:** rsync **only this folder** (e.g. `deepti/codes_LR_Haar/`).  
Low-rank math lives in **`codes_LR_Haar/lowrank/`** (not `probability_nbar_J2`).

Walrus reference jobs remain in repo-root `CLUSTER_J2_*` / `CLUSTER_J4_*`.  
Verified reference outputs from the sibling tree stay in **`codes_LR/`** only.

## Layout

```text
codes_LR_Haar/
  lowrank/                          # cluster kernel (Numba rank 1–4)
    fast_loop_hafnian.py
    lr_production_kernel.py
    cluster_reps.py
  ensemble_common.py
  helping_code_doc_Haar_parallel_architecture.md
  _lr_run_parallel_j2.py              # template → copied into each J=2 folder
  _lr_run_parallel_j4.py
  sync_from_walrus_cluster.py
  sync_lowrank_kernel.py              # refresh fast_loop from laptop P(n̄) tree
  validate_j2_all.py
  validate_j4_geometry_A.py
  # J=2 (6 jobs)
  CLUSTER_J2_superposition_LR_Haar/
  CLUSTER_J2_even_cat_LR_Haar/
  CLUSTER_J2_odd_cat_LR_Haar/
  CLUSTER_J2_Kerr_squeezed_Haar/
  CLUSTER_J2_Kerr_squeezed_kappa_1_Haar/
  # J=4 (6 jobs)
  CLUSTER_J4_superposition_LR_Haar/
  CLUSTER_J4_even_cat_LR_Haar/
  CLUSTER_J4_odd_cat_LR_Haar/
  CLUSTER_J4_Kerr_squeezed_Haar/
  CLUSTER_J4_Kerr_squeezed_kappa_1_Haar/
  Standard_Gaussian/
    CLUSTER_J2_Gaussian_squeezed_Haar/
    CLUSTER_J4_Gaussian_squeezed_Haar/
```

Each `CLUSTER_*_Haar/` folder is a self-contained cluster job:

- `run_parallel_J2_even_cat_LR_Haar.py`, `submit_J2_even_cat_LR_Haar.slurm`
- `config.py`: `INTERFEROMETER`, `HAAR_BASE_SEED`, `N_ENSEMBLE`, state parameters
- `output/Pnbar_*.json`

Runners add the parent tree to `sys.path` (`HERE.parent` from a state folder).

## Cerit rsync (one tree)

```bash
rsync -avz codes_LR_Haar/  user@zenith:~/deepti/codes_LR_Haar/
```

Then:

```bash
cd ~/deepti/codes_LR_Haar/CLUSTER_J2_superposition_LR_Haar
mkdir -p logs output && qsub submit_J2_superposition_LR_Haar.slurm
```

The Kerr-squeezed jobs use `κ=0.1`, `sinh(r)=1`
(`r=arcsinh(1)≈0.8814`), `INPUT_NMAX=5`, and output cutoff 6. Submit them
from their cluster folders in the same way:

```bash
cd ~/deepti/codes_LR_Haar/CLUSTER_J2_Kerr_squeezed_Haar
qsub submit_J2_Kerr_squeezed_LR_Haar.slurm
```

## After changing rank-4 / Takagi on laptop

```bash
python codes_LR_Haar/sync_lowrank_kernel.py
python codes_LR_Haar/sync_from_walrus_cluster.py   # if cluster runners changed
```

## Local validation

```bash
python codes_LR_Haar/validate_j2_all.py
python codes_LR_Haar/validate_j4_geometry_A.py
```

## Maintenance

- Edit parallel logic: `_lr_run_parallel_j2.py` / `_j4.py` → `sync_from_walrus_cluster.py`
- Edit Walrus physics in repo `CLUSTER_*` → sync Haar folders from Walrus templates
- Do **not** edit `codes_LR/` when working on Haar — keep the sibling reference tree frozen


## Haar ensemble (cluster)

All ensemble knobs live in **each job's `config.py` only**:

```python
HAAR_BASE_SEED = 20250810
N_ENSEMBLE = 1000         # production ensemble size
ENSEMBLE_BATCH_SIZE = 10
ZERO_SAVE_TOL = 1e-10
```

Submit via the wrapper (array range is derived from `config.py` — do not edit `#PBS -J` by hand):

```bash
cd CLUSTER_J2_even_cat_LR_Haar
mkdir -p logs output
./submit_ensemble.sh
# or locally (runs all N_ENSEMBLE from config):
python run_ensemble.py --workers 8
```

Outputs (folder name embeds `N_ENSEMBLE` and `HAAR_BASE_SEED`):

```text
output/ensemble/R10_base20250810/seed_000000.json    # pilot (if retained)
output/ensemble/R1000_base20250810/seed_000042.json  # production
```

Each JSON stores **sparse** `P(n̄)` (`ZERO_SAVE_TOL`) plus metadata (`seed`, `interferometer`, …).
Seed rule: `seed = HAAR_BASE_SEED + realization_index`.
Correlations / ensemble statistics are **not** computed on the cluster.

Pilot → production: change **only** `N_ENSEMBLE` in each `config.py`, then `./submit_ensemble.sh` again.
Production is currently set to `N_ENSEMBLE = 1000`.
