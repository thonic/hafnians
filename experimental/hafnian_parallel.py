from joblib import Parallel, delayed
import numpy as np
from thewalrus.decompositions import takagi
from sympy import symbols, expand, Poly
from scipy.special import factorial2
import time

def compute_term_value(key, coeff, G):
    """Compute the value of a term given its exponents, coefficient, and matrix G."""
    # Adjust exponents by subtracting 1
    adjusted_exponents = [exp - 1 for exp in key]

    # Compute double factorial for adjusted exponents
    double_factorials = [factorial2(exp) if exp >= 0 else 1 for exp in adjusted_exponents]

    # Multiply the double factorials and coefficient
    term_value = coeff
    for df in double_factorials:
        term_value *= df

    return term_value

def generate_hafnian(G):
    n, r = G.shape  # Get the dimensions of G

    # If n is odd, return 0 as specified
    if n % 2 != 0:
        return 0

    # Define symbolic variables x0, x1, ..., xr-1
    x = symbols(f"x0:{r}")
    poly = 1.0

    # Construct the multivariate polynomial
    for i in range(n):
        term = sum(G[i, j] * x[j] for j in range(r))
        poly = expand(poly * term)

    # Convert polynomial to dictionary
    p = Poly(poly, x)  # Polynomial with respect to all x variables
    terms = p.as_dict()

    # Filter terms where sum of exponents is equal to n and each exponent is even
    valid_terms = {k: v for k, v in terms.items() if sum(k) == n and all(exp % 2 == 0 for exp in k)}

    # Parallel computation of term values
    results = Parallel(n_jobs=-1)(
        delayed(compute_term_value)(key, coeff, G) for key, coeff in valid_terms.items()
    )

    # Sum the results
    return sum(results)

# Parameters
n = 12  # Reduced matrix size for debugging
r = 6  # Reduced rank for debugging

def generate_random_symmetric_matrix(n, r):
    """Generate a symmetric matrix of size n x n with rank r."""
    A = np.random.randn(n, n)
    A = (A + A.T) / 2  # Make symmetric
    U, S, Vt = np.linalg.svd(A)
    S[r:] = 0  # Set singular values beyond rank r to zero
    return U @ np.diag(S) @ Vt  # Reconstruct A

# Generate A and G matrices
start_time = time.time()
while True:
    A = generate_random_symmetric_matrix(n, r)
    if np.linalg.matrix_rank(A) != r:
        continue
    singular_values, U = takagi(A, svd_order=True, rtol=1e-16)
    G = U[:, :r] @ np.diag(np.sqrt(singular_values[:r]))
    G = np.real_if_close(G, tol=1e-10)
    if np.allclose(A, G @ G.T, atol=1e-8):
        break

# Manual Hafnian Calculation
start_manual = time.time()
hafnian_manual = generate_hafnian(G)
end_manual = time.time()

print(f"Hafnian__manual: {hafnian_manual}")
print(f"Manual computation completed in {end_manual - start_manual:.2f} seconds.")

# Validate using The Walrus
from thewalrus import hafnian
start_walrus = time.time()
hafnian_walrus = hafnian(A)
end_walrus = time.time()

print(f"Hafnian__walrus: {hafnian_walrus}")
print(f"The Walrus computation completed in {end_walrus - start_walrus:.2f} seconds.")
