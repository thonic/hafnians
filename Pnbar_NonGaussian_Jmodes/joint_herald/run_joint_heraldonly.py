"""
Joint numerator + P_herald only.

Pure P(n̄)  →  ../probability_nbar_J2/
"""

import numpy as np

from code_joint_herald_fiurasek import compute_joint, generate_sigma_D
from state_definitions import CP_SUPERPOSITION
# from state_definitions import CP_SECONDARY

# =============================================================================
# CHANGE HERE
# =============================================================================

cp = CP_SUPERPOSITION
# cp = CP_SECONDARY                         # use K = 1 only

K = 1
apply_hadamard = True

# =============================================================================
# RUN
# =============================================================================

cp_norm = cp / np.linalg.norm(cp)
st = generate_sigma_D(cp=cp_norm)

joint = compute_joint(
    j=K,
    sigma_single=st.sigma,
    d_single=st.d,
    apply_hadamard=apply_hadamard and K > 1,
    include_p_herald=(K == 1),
)

print("=" * 70)
print("RUN SUMMARY")
print("=" * 70)
print(f"K (copies)        = {K}")
print(f"Σ network size    = ({4 * K}, {4 * K})")
print(f"cp (normalised)   = {np.round(st.cp, 6)}")
print(f"apply_hadamard    = {apply_hadamard and K > 1}")
print()
print(joint)
print("=" * 70)
