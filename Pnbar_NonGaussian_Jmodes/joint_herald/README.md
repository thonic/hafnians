# joint_herald — herald numerator + P_herald

Self-contained folder for **joint heralded** calculations:

    ⟨n̂ | herald⟩ numerator, P_herald, ratio numerator/P_herald

**Not** P(n̄). For pure P(n̄) see `../probability_nbar_J2/` or `../CLUSTER_INPUT/INPUT_THEWALRUS/`.

## Run

```bash
cd joint_herald
python run_joint_heraldonly.py
python verify_network.py
python verify_secondary_state.py
```

Expected (superposition, K=1): P_herald ≈ 4×10⁻⁸, ratio ≈ 0.5.
