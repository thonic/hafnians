# Haar production code — structure & parallel architecture

**Purpose:** Reference for thesis/paper methods, cluster workflow, and “where does parallelism happen?”  
**Tree:** `codes_LR_Haar/` (do **not** edit frozen sibling `codes_LR/` for Haar work)  
**Last aligned with:** ensemble Option B (PBS array + ProcessPoolExecutor, no nested pools)

---

## How to use this document

1. **Before a cluster run** — read [Cluster workflow](#cluster-workflow-what-you-run) and [Ensemble parameters](#ensemble-parameters-configpy).
2. **Before production (1000 seeds)** — read [Pilot vs production](#pilot-vs-production) and [PBS array math](#pbs-array-level).
3. **For a methods section** — copy/adapt [Methods summary (paper-ready)](#methods-summary-paper-ready).
4. **When debugging slowness** — read [Efficiency notes](#efficiency-notes) and the ASCII timelines below.
5. **Ask Cursor** — “Explain X using `helping_code_doc_Haar_parallel_architecture.md`” to get answers grounded in this file.

---

## Repository layout (high level)

```text
codes_LR_Haar/
├── helping_code_doc_Haar_parallel_architecture.md   ← YOU ARE HERE
├── README.md
├── ensemble_common.py              # shared ensemble I/O, seeds, batch ranges
├── lowrank/                        # Numba low-rank hafnian kernel (cluster copy)
├── submit_ensemble_template.sh     # template copied into each CLUSTER_* folder
├── run_ensemble_j2_template.py
├── run_ensemble_j4_template.py
├── _lr_run_parallel_j2.py          # edit parallel logic here → sync to folders
├── _lr_run_parallel_j4.py
│
├── CLUSTER_J2_*_Haar/              # 5 non-Gaussian J=2 jobs
├── CLUSTER_J4_*_Haar/              # 5 non-Gaussian J=4 jobs
└── Standard_Gaussian/
    ├── CLUSTER_J2_Gaussian_squeezed_Haar/
    └── CLUSTER_J4_Gaussian_squeezed_Haar/
```

Each **`CLUSTER_*_Haar/`** folder is a **self-contained PBS job**. Typical contents:

| File | Role |
|------|------|
| `config.py` | Physics + **ensemble knobs** (`N_ENSEMBLE`, batch size, seeds) |
| `build_B_mu_scale.py` | Haar unitary, B matrix, μ, scale for one realization |
| `run_parallel_J*.py` | **Pattern-level** parallelism (`ProcessPoolExecutor`) |
| `run_ensemble.py` | **Batch driver**: loops seeds, calls `run_parallel_*.run()` |
| `submit_ensemble.sh` | **YOU RUN THIS** on cluster — reads `config.py`, calls `qsub` |
| `submit_ensemble.slurm` | PBS script — **do not run directly**; invoked by wrapper |
| `submit_J*_*.slurm` | Legacy single-shot job (one Haar draw, not full ensemble) |
| `output/ensemble/R{N}_base{seed}/` | Sparse JSON per realization |

**12 cluster folders** total (10 under root + 2 under `Standard_Gaussian/`):
J=2 and J=4 for Gaussian, even/odd cat, Kerr κ=0.1, Kerr κ=1, and superposition.

---

## Cluster workflow: what you run

### Correct command (ensemble production)

```bash
cd ~/deepti/codes_LR_Haar/CLUSTER_J2_even_cat_LR_Haar   # example folder
mkdir -p logs output
./submit_ensemble.sh
```

### Do **not** run this for ensemble submission

```bash
qsub submit_ensemble.slurm    # WRONG as primary entry — skips array logic
```

**Why:** `submit_ensemble.sh` reads `N_ENSEMBLE` and `ENSEMBLE_BATCH_SIZE` from `config.py`, computes how many PBS array tasks are needed, and calls `qsub` correctly:

- **1 batch** (e.g. pilot 10/10) → `qsub submit_ensemble.slurm` (no `-J`; OpenPBS rejects `-J 0-0`)
- **N batches** (e.g. 1000/10) → `qsub -J "0-99" submit_ensemble.slurm`

The `.slurm` file is the **PBS job body**; the `.sh` wrapper is the **submission interface**.

### Deploy to Cerit

```bash
rsync -avz codes_LR_Haar/  user@zenith:~/deepti/codes_LR_Haar/
```

### Local test (laptop, no PBS)

```bash
cd CLUSTER_J2_even_cat_LR_Haar
python run_ensemble.py --workers 8
# runs ALL N_ENSEMBLE seeds sequentially on your machine
```

### Single realization / legacy single-shot

```bash
qsub submit_J2_even_cat_LR_Haar.slurm
# OR
python run_parallel_J2_even_cat_LR_Haar.py --workers 8
```

---

## What PBS does here

PBS (Portable Batch System) is the cluster job scheduler on Cerit.

```text
You                          PBS scheduler                 Compute node
───                          ───────────────               ────────────
./submit_ensemble.sh  ──►    qsub registers job(s)  ──►  bash submit_ensemble.slurm
                             assigns node + CPUs            │
                                                             ▼
                                                      python run_ensemble.py
                                                             │
                                                             ▼
                                                      run_parallel_*.py (pool)
```

**Key PBS variables** (set automatically on the node):

| Variable | Meaning |
|----------|---------|
| `PBS_ARRAY_INDEX` | Which batch (0, 1, …, 99 for 100 batches) |
| `PBS_NCPUS` / `PBS_NP` | CPUs allocated (`select=1:ncpus=8` → 8) |
| `PBS_O_WORKDIR` | Job directory |
| `PBS_JOBID` | Unique job id |

**In `submit_ensemble.slurm`:**

```bash
#PBS -l select=1:ncpus=8:mem=8GB    # 1 node, 8 cores, 8 GB RAM
#PBS -l walltime=01:00:00           # max runtime per array task

python run_ensemble.py --batch-index "${PBS_ARRAY_INDEX}" --workers "${NCPUS}"
```

So: **PBS parallelizes batches across nodes**; **Python parallelizes patterns within one seed** on the allocated cores.

---

## Ensemble parameters (`config.py`)

Single source of truth in **each** cluster folder:

```python
HAAR_BASE_SEED = 20250810
N_ENSEMBLE = 1000             # production ensemble size
ENSEMBLE_BATCH_SIZE = 10       # realizations per PBS array task
ZERO_SAVE_TOL = 1e-10          # sparse storage threshold
INTERFEROMETER = "haar"
```

**Seed rule:** `seed = HAAR_BASE_SEED + realization_index`  
**Output dir:** `output/ensemble/R{N}_base{HAAR_BASE_SEED}/seed_XXXXXX.json`  
**J=4:** also `seed_XXXXXX_geomA.json`, `_geomB`, `_geomC` (geometries run sequentially per seed).

---

## Two-level parallel architecture (overview)

```text
LEVEL 1 — PBS job array (cluster / across nodes)
═══════════════════════════════════════════════════

  Batch 0          Batch 1          Batch 2     ...     Batch 99
 (seeds 0–9)    (seeds 10–19)   (seeds 20–29)         (seeds 990–999)
     │               │               │                      │
     ▼               ▼               ▼                      ▼
  Node A?         Node B?         Node C?                Node …?
  (can run simultaneously if scheduler allows)


LEVEL 2 — ProcessPoolExecutor (inside one PBS task / one node)
══════════════════════════════════════════════════════════════

  One realization (one Haar seed):
  ┌─ build B, μ, Takagi G  (serial, main process) ─┐
  └──────────────────────────────────────────────────┘
                          │
                          ▼
              ProcessPoolExecutor (8 workers)
              ┌─────────────────────────────┐
              │ Pattern worker 1  n̄=(0,0,…) │
              │ Pattern worker 2  n̄=(0,0,1) │
              │ Pattern worker 3  n̄=(0,0,2) │
              │ …                             │
              │ Pattern worker 8  n̄=…        │
              └─────────────────────────────┘
                          │
                          ▼
              normalize + sparse JSON write (serial)
```

**Explicitly NOT used:** nested multiprocessing (parallel seeds × parallel patterns in one process tree).

---

## Realization timeline inside one PBS job

### Pilot: `N_ENSEMBLE = 10`, `ENSEMBLE_BATCH_SIZE = 10`

→ **1 PBS job**, no array. All 10 seeds **sequentially** in that job:

```text
PBS job (batch_index = 0, 8 cores allocated)
│
├── Seed 0  (realization_index 0)
│       │
│       │   8 cores
│       │  ┌───────────────┐
│       │  │Pattern worker1│
│       │  │Pattern worker2│
│       │  │Pattern worker3│
│       │  │ ...           │
│       │  │Pattern worker8│
│       │  └───────────────┘
│       ▼
│   seed_000000.json
│
├── Seed 1  (index 1)     ← same 8-core pool, new ProcessPoolExecutor
│       │
│       │   8 cores
│       ▼
│   seed_000001.json
│
├── Seed 2
│       │
│       │   8 cores
│       ▼
│   ...
│
└── Seed 9
        ▼
    seed_000009.json
```

### Production: `N_ENSEMBLE = 1000`, `ENSEMBLE_BATCH_SIZE = 10`

→ **100 PBS array tasks** (`qsub -J "0-99"`). Each task runs **10 seeds sequentially**:

```text
PBS_ARRAY_INDEX = 0     PBS_ARRAY_INDEX = 1        ...    PBS_ARRAY_INDEX = 99
seeds 0–9               seeds 10–19                       seeds 990–999
(sequential)            (sequential)                      (sequential)
     │                       │                                  │
     ▼                       ▼                                  ▼
  may run on              may run on                         may run on
  node A                  node B                             node Z
  at the same time ───────────────────────────────────────────────►
                    (if queue has free nodes)
```

| Quantity | Value |
|----------|-------|
| PBS jobs submitted | **100** |
| Realizations per PBS job | **10** |
| Total realizations | 100 × 10 = **1000** |
| Cores used per realization | **8** (pattern pool) |
| Cores used across seeds inside one job | **still 8** (seeds are serial) |

---

## 1. Pattern-level parallelism (detail)

### Where `ProcessPoolExecutor` is created

File: `run_parallel_J*.py` in each cluster folder (e.g. `run_parallel_J2_Kerr_squeezed_LR_Haar.py`).

```python
with ProcessPoolExecutor(
    max_workers=n_workers,
    initializer=_init_worker,
    initargs=(G, mu, scale, MK_dict),
) as pool:
    for nbar, p in pool.map(_amplitude_prob, patterns, chunksize=chunksize):
        ...
```

### What runs in parallel

Function **`_amplitude_prob(nbar)`** — one photon-number pattern:

1. `reps_for_nbar` → repetition vector  
2. `loop_hafnian_lr_from_G(G, mu, reps)` → amplitude  
3. return `|amplitude|²`

### What is distributed

All patterns: `itertools.product(range(cutoff), repeat=J)`

| J | CUTOFF | Patterns per geometry |
|---|--------|-------------------------|
| 2 | 6 | 36 |
| 4 | 6 | 1296 |

### One-time serial work (per seed)

Before the pool (main process only):

- `build_B_mu_scale(..., haar_seed=seed)` — **Haar random unitary** for this realization  
- `takagi_factor_b_mat(b_mat)` — Takagi factor **G**

Workers receive **G, μ, scale, MK_dict** via `_init_worker` — Takagi is **not** recomputed per pattern.

After the pool (serial): normalization, `sparse_probability_list`, atomic JSON write.

### Are all 8 CPUs used?

- `submit_ensemble.slurm` passes `--workers "${NCPUS}"` (8 from `#PBS -l select=1:ncpus=8`).
- `resolve_workers()` in `run_parallel_*.py` honors that.
- **Utilization:** best for J=4 (1296 patterns); J=2 (36 patterns) still uses the pool but with coarser task granularity.

---

## 2. Realization-level execution

**Driver:** `run_ensemble.py`

```python
for idx in indices:
    seed = HAAR_BASE_SEED + idx
    runner.run(workers=..., haar_seed=seed, out_json=...)
```

- **One** `ProcessPoolExecutor` per seed (created inside `runner.run()`, destroyed when done).
- Realizations in a batch: **strictly sequential**.
- Valid existing JSON → **SKIP** (restart-safe).

**J=4 extra loop** (geometries A, B, C sequential per seed):

```python
for idx in indices:
    for g in geom_list:
        runner.run_geometry(g, ...)
```

---

## 3. PBS array level (formulas)

```text
N_BATCHES = ceil(N_ENSEMBLE / ENSEMBLE_BATCH_SIZE)
LAST_INDEX = N_BATCHES - 1
```

Batch `k` (`PBS_ARRAY_INDEX = k`) runs realization indices:

```text
start = k * ENSEMBLE_BATCH_SIZE
end   = min(start + ENSEMBLE_BATCH_SIZE, N_ENSEMBLE)
indices = range(start, end)
```

**Example (1000 / 10):**

| PBS_ARRAY_INDEX | Realization indices | Seeds |
|-----------------|---------------------|-------|
| 0 | 0–9 | base+0 … base+9 |
| 1 | 10–19 | … |
| … | … | … |
| 99 | 990–999 | … |

**Can multiple array jobs run at once?** Yes — each array element is an independent PBS job; the scheduler may place them on different nodes concurrently (subject to queue limits).

---

## 4. Overall architecture — confirmed

| Level | Mechanism | Parallelizes |
|-------|-----------|--------------|
| **1** | PBS job array (`./submit_ensemble.sh`) | Batches of Haar realizations across cluster nodes |
| **2** | `ProcessPoolExecutor` in `run_parallel_*.py` | Photon-number patterns within one realization |

**Not parallelized (by design):**

- Realizations within one batch  
- J=4 geometries A/B/C within one seed  
- Takagi / B,μ build within one seed  

From `run_ensemble.py` docstring: *“No nested multiprocessing.”*

---

## 5. Efficiency notes

**Why this design is good for large ensembles**

- Realizations are embarrassingly parallel → PBS array scales to many nodes.  
- Avoids nested pools → no 8×8 oversubscription on 8 cores.  
- One JSON per seed → cheap restart (skip valid files).  
- Predictable memory: one realization + one pool at a time per PBS job.

**Why nested multiprocessing would be worse**

- Oversubscription (many more processes than cores).  
- Higher RAM (each worker holds G, μ, …; multiple realizations × multiple pools).  
- Heavier fork/pickle overhead.  
- Harder to size jobs for the scheduler.

**Practical tuning:** if one batch is too slow for `walltime`, lower `ENSEMBLE_BATCH_SIZE` (more PBS tasks, fewer seeds per task). Outer parallelism increases; inner pattern parallelism unchanged.

---

## Call chain (one ensemble batch)

```text
./submit_ensemble.sh
    └── qsub [-J 0-(N-1)] submit_ensemble.slurm
            └── python run_ensemble.py --batch-index K --workers 8
                    └── for idx in batch_range:
                            └── run_parallel_*.run(haar_seed=base+idx)
                                    ├── build_B_mu_scale (Haar U)
                                    ├── takagi_factor_b_mat
                                    └── ProcessPoolExecutor.map(_amplitude_prob, patterns)
                                            └── loop_hafnian_lr_from_G (lowrank/)
                                    └── write seed_XXXXXX.json (sparse)
```

---

## Output files

```text
output/ensemble/R10_base20250810/
├── manifest.json           # run metadata (overwritten each batch start)
├── seed_000000.json
├── seed_000001.json
└── ...
```

Each `seed_*.json` contains:

- `seed`, `haar_base_seed`, `interferometer`, `state`, `J`, `cutoff`  
- `probabilities`: sparse `[[nbar, P], ...]` with `P >= ZERO_SAVE_TOL`  
- `sum_P`, timings, `zero_save_tol`  

**Correlations / ensemble statistics:** computed **offline** on laptop — not on cluster.

---

## Pilot vs production

| Setting | Pilot | Production |
|---------|-------|------------|
| `N_ENSEMBLE` | 10 | **1000** |
| `ENSEMBLE_BATCH_SIZE` | 10 | 10 (unchanged) |
| PBS mode | single job | array `0-99` |
| Change in code | — | **`config.py` only** |
| Submit | `./submit_ensemble.sh` | same |

---

## Methods summary (paper-ready)

> Each Haar ensemble realization corresponds to an independent draw of a random passive unitary on the signal modes, identified by a deterministic seed `HAAR_BASE_SEED + realization_index`. For a given realization, the low-rank Takagi factorization of the interferometer-transformed covariance is computed once; photon-number amplitudes P(n̄) are then evaluated at all cutoff patterns in parallel using a Python `ProcessPoolExecutor` over the allocated CPU cores. Probabilities are normalized over the full pattern set and stored sparsely (entries with P(n̄) ≥ 10⁻¹⁰ only). For large ensembles, realizations are grouped into batches of size `ENSEMBLE_BATCH_SIZE`; batches are submitted as a PBS job array so that independent batch jobs can run concurrently on different compute nodes. Realizations within a batch are executed sequentially, avoiding nested process pools. This yields two levels of parallelism—cluster-level batching over Haar realizations and node-level parallelization over photon-number patterns—without oversubscribing CPU resources.

---

## Quick FAQ

**Q: Which file has the Haar random matrix?**  
A: `build_B_mu_scale.py` — uses `scipy.stats.unitary_group` with `haar_seed`.

**Q: Where is parallelism?**  
A: Level 1 → PBS array in `submit_ensemble.sh` / `.slurm`. Level 2 → `ProcessPoolExecutor` in `run_parallel_*.py`.

**Q: What do I run on the cluster?**  
A: `./submit_ensemble.sh` from the desired `CLUSTER_*_Haar/` folder.

**Q: Why is my pilot (10 seeds) slow?**  
A: One PBS job runs all 10 seeds **one after another**; J=4 is ~100× heavier per seed than J=2. Production speed comes from **100 parallel PBS tasks**, not one giant job.

**Q: Can I resume after failure?**  
A: Yes. Re-run `./submit_ensemble.sh`; completed valid JSON files are skipped.

---

## File index (parallelism-related)

| File | Parallel role |
|------|----------------|
| `submit_ensemble.sh` | Computes batch count; submits PBS (array or single) |
| `submit_ensemble.slurm` | PBS resource request; launches `run_ensemble.py` |
| `run_ensemble.py` | Serial loop over seeds in one batch |
| `run_parallel_*.py` | **ProcessPoolExecutor** over patterns |
| `ensemble_common.py` | `batch_realization_range`, paths, sparse I/O |
| `config.py` | `N_ENSEMBLE`, `ENSEMBLE_BATCH_SIZE`, `HAAR_BASE_SEED` |
| `build_B_mu_scale.py` | Haar unitary + B, μ (serial, per seed) |
| `lowrank/lr_production_kernel.py` | Takagi + loop hafnian (called from workers) |
