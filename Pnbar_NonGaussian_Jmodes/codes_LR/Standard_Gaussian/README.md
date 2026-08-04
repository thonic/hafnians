# Standard Gaussian squeezed (The Walrus low-rank hafnian)

Self-contained cluster jobs for **standard Gaussian squeezed vacuum**
\(|S(r)\rangle\) with \(r=\operatorname{asinh}(1)\), using the same Fiurášek
pipeline, cutoffs, interferometer, and geometries as the even/odd-cat LR
jobs.

## Hafnian engine

- Uses the **ordinary low-rank hafnian** of \(A = GG^{T}\) (Walrus
  `low_rank_hafnian` semantics / Björklund App. C).
- Does **not** use the loop hafnian (`loop_hafnian_lr_from_G`).
- Takagi factor \(G\) of the Fiurášek \(B\) matrix is expanded by herald
  repetition vectors; displacements \(\mu\) are not used.
- Production calls the Numba ordinary-hafnian kernel (loop hafnian with
  \(\mu\equiv 0\)), which matches `thewalrus.hafnian(G@G.T)`. Stock
  `thewalrus.low_rank_hafnian` is SymPy-based and too slow for full
  pattern sweeps; v0.21 also has a `factorial2(-1)` bug.

## Layout

```text
Standard_Gaussian/
  walrus_low_rank_hafnian.py
  CLUSTER_J2_Gaussian_squeezed/
  CLUSTER_J4_Gaussian_squeezed/
```

## Submit

```bash
cd CLUSTER_J2_Gaussian_squeezed && mkdir -p logs output && qsub submit_J2_Gaussian_squeezed_LR.slurm
cd CLUSTER_J4_Gaussian_squeezed && mkdir -p logs output && qsub submit_J4_Gaussian_squeezed_LR.slurm
```
