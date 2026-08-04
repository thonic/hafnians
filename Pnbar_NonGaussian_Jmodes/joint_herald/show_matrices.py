"""
Copy of hafnians/show_matrices(sigma,D).py — prints Σ and D for the Fiurásek state.

Run from this folder:
    python show_matrices.py
"""

import math
import numpy as np
from generate_cv_and_dv import generate_cv_and_dv

# --- Common setup code needed for both versions ---
def setup_matrices(n):
    nm = n - 1
    maxn = 2 * n
    am = np.zeros((maxn + 2, maxn + 2))
    ad = np.zeros((maxn + 2, maxn + 2))
    for j in range(maxn + 1):
        am[j][j + 1] = np.sqrt(j + 1)
        ad[j + 1][j] = np.sqrt(j + 1)
    sq = np.arcsinh(1)
    A = np.cosh(sq) * am + ad
    Am = {}
    for j in range(n + 1):
        Aj = np.linalg.matrix_power(A, j)
        Am[j] = Aj[:n, :n]
    cvac = np.zeros(n)
    cvac[0] = 1
    return Am, cvac, nm

# --- VERSION 1: YOUR ORIGINAL CODE from generate_displacements.py ---
def generate_alpha_original(cp, t=0.99999999):
    n = len(cp)
    Am, cvac, nm = setup_matrices(n)
    hvec = np.zeros(n, dtype="complex_")
    h = hvec
    cp2 = np.copy(cp)
    for j in range(n):
        h[j] = cp2[n - j - 1] / np.sqrt(math.factorial(nm - j))
        cp2 = cp2 - h[j] * Am[n - j - 1] @ cvac
    beta = np.roots(h)
    m = np.zeros((n - 1, n - 1))
    for j in range(n - 1, 0, -1):
        for k in range(n - 1, j - 1, -1):
            m[j - 1][k - 1] = np.power(t, (n - 1) - k)
    m_inv = np.linalg.inv(m)
    alpha = np.zeros(n, dtype="complex_")
    if n > 1:
        alpha[1:] = m_inv @ beta
    sq = np.arcsinh(1)
    s1, s2 = 0, 0
    for j in range(1, n):
        s1 = s1 + alpha[j] * np.power(t, n - j)
        s2 = s2 + np.conj(alpha[j]) * np.power(t, j - n)
    s1 = s1 * np.cosh(sq)
    sdiff = (s2 - s1) / np.cosh(sq)
    x = np.real(sdiff) / (np.power(t, n) - np.power(t, -n) / np.cosh(sq))
    y = np.imag(sdiff) / (np.power(t, n) + np.power(t, -n) / np.cosh(sq))
    alpha[0] = x + 1j * y
    return alpha

# --- VERSION 2: CORRECTED CODE based on the Fiurasek 2005 paper ---
def generate_alpha_corrected(cp, t=0.99999999):
    n = len(cp)
    Am, cvac, nm = setup_matrices(n)
    h = np.zeros(n, dtype="complex_")
    psi_j = np.copy(cp)
    for j in range(n):
        h[j] = psi_j[n - j - 1] / np.sqrt(math.factorial(nm - j))
        psi_j = psi_j - h[j] * (Am[n - j - 1] @ cvac)
    beta = np.roots(h)
    m = np.zeros((n - 1, n - 1))
    for j in range(n - 1, 0, -1):
        for k in range(n - 1, j - 1, -1):
            m[j - 1][k - 1] = np.power(t, (n - 1) - k)
    m_inv = np.linalg.inv(m)
    alpha = np.zeros(n, dtype="complex_")
    if n > 1:
        alpha[1:] = m_inv @ beta
    sq = np.arcsinh(1)
    s1, s2 = 0, 0
    for j in range(1, n):
        s1 = s1 + alpha[j] * np.power(t, n - j)
        s2 = s2 + np.conj(alpha[j]) * np.power(t, j - n)
    s1 = s1 * np.cosh(sq)
    sdiff = (s2 - s1) / np.cosh(sq)
    x = np.real(sdiff) / (np.power(t, n) - np.power(t, -n) / np.cosh(sq))
    y = np.imag(sdiff) / (np.power(t, n) + np.power(t, -n) / np.cosh(sq))
    alpha[0] = x + 1j * y
    return alpha

# --- Main Test Execution ---
# 1. Define the target state |0> + |1>
target_cp = np.array([1, 1], dtype=np.complex128)
target_cp_normalized = target_cp / np.sqrt(np.sum(np.abs(target_cp)**2))

print("=======================================================")
print("           ALPHA VECTOR COMPARISON TEST")
print("=======================================================")
print("Target State: (|0> + |1>) / sqrt(2)")
print("-------------------------------------------------------")

# 2. Run YOUR ORIGINAL version
alpha_original = generate_alpha_original(target_cp_normalized)
print("\nALPHA from YOUR ORIGINAL code:")
print(np.round(alpha_original, 8))
print("(This is the source of the incorrect D vector you are seeing)")

# 3. Run the CORRECTED version
alpha_corrected = generate_alpha_corrected(target_cp_normalized)
print("\nALPHA from the CORRECTED code (based on paper):")
print(np.round(alpha_corrected, 8))
print("(This produces the correct D vector)")
print("-------------------------------------------------------")

# 4. Now, calculate and print the final matrices using the CORRECTED alpha
print("\n=======================================================")
print("   FINAL MATRICES from CORRECTED ALPHA")
print("=======================================================")

cv_final, dv_final = generate_cv_and_dv(alpha_corrected, K=1, M=2, N=2, single_mode=True)

print("\nCORRECT Covariance Matrix (SIGMA):")
print(np.round(cv_final, 8))

print("\nCORRECT Displacement Vector (D):")
print(np.round(dv_final, 8))
print("-------------------------------------------------------")