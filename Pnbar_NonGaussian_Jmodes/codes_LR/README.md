# Low-rank cluster jobs (`LR_Cluster_Codes`)

**Self-contained for Cerit:** rsync **only this folder** (`deepti/LR_Cluster_Codes/`).  
Low-rank math lives in **`LR_Cluster_Codes/lowrank/`** (not `probability_nbar_J2`).

Walrus reference jobs remain in repo-root `CLUSTER_J2_*` / `CLUSTER_J4_*`.

## Layout

```text
LR_Cluster_Codes/
  lowrank/                          # cluster kernel (Numba rank 1–4)
    fast_loop_hafnian.py
    lr_production_kernel.py
    cluster_reps.py
  _lr_run_parallel_j2.py              # template → copied into each J=2 folder
  _lr_run_parallel_j4.py
  sync_from_walrus_cluster.py
  sync_lowrank_kernel.py              # refresh fast_loop from laptop P(n̄) tree
  validate_j2_all.py
  validate_j4_geometry_A.py
  CLUSTER_J2_superposition_LR/
  CLUSTER_J2_even_cat_LR/
  CLUSTER_J2_odd_cat_LR/
  CLUSTER_J2_Kerr_squeezed/
  CLUSTER_J4_even_cat_LR/
  CLUSTER_J4_odd_cat_LR/
  CLUSTER_J4_Kerr_squeezed/
```

Each `CLUSTER_*_LR/` folder matches Walrus naming:

- `run_parallel_J2_even_cat_LR.py`, `submit_J2_even_cat_LR.slurm`
- `output/Pnbar_J2_even_cat_LR.json`

Runners add **`LR_Cluster_Codes`** to `sys.path` (`HERE.parent` from a state folder).

## Cerit rsync (one tree)

```bash
rsync -avz LR_Cluster_Codes/  user@zenith:~/deepti/LR_Cluster_Codes/
```

Then:

```bash
cd ~/deepti/LR_Cluster_Codes/CLUSTER_J2_superposition_LR
mkdir -p logs output && qsub submit_J2_superposition_LR.slurm
```

The Kerr-squeezed jobs use `κ=0.1`, `sinh(r)=1`
(`r=arcsinh(1)≈0.8814`), `INPUT_NMAX=5`, and output cutoff 6. Submit them
from their cluster folders in the same way:

```bash
cd ~/deepti/LR_Cluster_Codes/CLUSTER_J2_Kerr_squeezed
qsub submit_J2_Kerr_squeezed_LR.slurm
```

## After changing rank-4 / Takagi on laptop

```bash
python LR_Cluster_Codes/sync_lowrank_kernel.py
python LR_Cluster_Codes/sync_from_walrus_cluster.py   # if cluster runners changed
```

## Local validation

```bash
python LR_Cluster_Codes/validate_j2_all.py
python LR_Cluster_Codes/validate_j4_geometry_A.py
```

## Maintenance

- Edit parallel logic: `_lr_run_parallel_j2.py` / `_j4.py` → `sync_from_walrus_cluster.py`
- Edit Walrus physics in repo `CLUSTER_*` → sync LR folders from Walrus templates
