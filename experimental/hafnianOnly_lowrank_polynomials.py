import numpy as np
from thewalrus.decompositions import takagi
from thewalrus import hafnian
import time
import numpy as np


def factorial2(m):
    """Double factorial of a nonnegative integer m: m * (m-2) * (m-4) * ..."""
    if m <= 1:
        return 1
    return m * factorial2(m - 2)


def multiply_polynomials(polyA, polyB):
    """
    Multiply two multivariate polynomials represented as dictionaries:
      polyA[ (e0,e1,...,er-1) ] = coefficientA
      polyB[ (f0,f1,...,fr-1) ] = coefficientB
    Returns new_poly, also a dict with exponent tuples -> numeric coefficient.
    """
    new_poly = {}
    for expA, coeffA in polyA.items():
        for expB, coeffB in polyB.items():
            # Add exponents componentwise
            exp_new = tuple(a + b for a, b in zip(expA, expB))
            coeff_new = coeffA * coeffB
            # Accumulate into new_poly
            new_poly[exp_new] = new_poly.get(exp_new, 0.0) + coeff_new
    return new_poly


def generate_hafnian(G):
    """
    Parameters:
      G: shape (n, r)
         We form the polynomial P(x_0,...,x_{r-1}) = ∏_{i=1..n}(Σ_{j=1..r} G[i,j]*x_j).
         Then we look for terms whose total exponent == n and each exponent is even,
         do a double-factorial adjustment, and sum up the result.

    Returns:
      Numeric hafnian-like value as a float.
    """
    n, r = G.shape

    # If n is odd, the result is 0 (as in your original code)
    if n % 2 != 0:
        return 0.0

    # Start with the constant polynomial "1", i.e. exponents (0,0,...,0) -> 1.0
    poly = {(0,) * r: 1.0}

    # Multiply in each row's linear term
    for i in range(n):
        # Build a dict for the polynomial row_poly = sum_j G[i,j] * x_j
        row_poly = {}
        for j in range(r):
            coeff = G[i, j]
            if abs(coeff) > 1e-15:  # skip near-zero entries
                exp_list = [0] * r
                exp_list[j] = 1
                row_poly[tuple(exp_list)] = coeff

        # Multiply our current polynomial by this row's polynomial
        poly = multiply_polynomials(poly, row_poly)

    # Filter terms where sum(exponents) == n and all exponents are even
    result = 0.0
    for exps, coeff in poly.items():
        if sum(exps) == n and all(e % 2 == 0 for e in exps):
            # Double-factorial part: each exponent e => factorial2(e-1) if e>=1
            tmp = coeff
            for e in exps:
                if e >= 1:
                    tmp *= factorial2(e - 1)
            result += tmp

    return result


def generate_random_symmetric_matrix(n, r):
    """Generate a symmetric matrix of size n x n with rank r."""
    A = np.random.randn(n, n)
    A = (A + A.T) / 2  # Make symmetric
    U, S, Vt = np.linalg.svd(A)
    S[r:] = 0  # Set singular values beyond rank r to zero
    return U @ np.diag(S) @ Vt  # Reconstruct A


# Parameters
n = 12  # Matrix size
r = 6  # Desired rank

start_time = time.time()

while True:
    # Step 1: Generate a symmetric matrix A with rank r
    A = generate_random_symmetric_matrix(n, r)

    # Step 2: Verify the rank of A
    if np.linalg.matrix_rank(A) != r:
        continue  # Retry if rank does not match

    # Step 3: Perform Takagi decomposition
    singular_values, U = takagi(A, svd_order=True)

    # Step 4: Construct G using U and square roots of singular values
    G = U[:, :r] @ np.diag(np.sqrt(singular_values[:r]))

    # Step 5: Ensure G is real (handle small imaginary parts due to precision)
    G = np.real_if_close(G, tol=1e-10)

    # Step 6: Verify A ≈ G @ G.T
    if np.allclose(A, G @ G.T, atol=1e-8):
        np.savez("hafnianOnly_2.npz", A=A, G=G)  # Save matrices to a file
        break  # Exit loop if condition is satisfied

datatest = np.load("hafnianOnly_2.npz")
A = datatest["A"]
G = datatest["G"]
# print(f"G matrix: {G}")

# Generate polynomial and calculate the result
my_haf = generate_hafnian(G)
print("Hafnian__code:", my_haf)

# Output timing
print(f"Process completed in {time.time() - start_time:.2f} seconds.")

# Calculate loop Hafnian using thewalrus
start_time2 = time.time()
# thewalrus_haf = loop_hafnian(A = A_reconstructed, D=Mu, reps=None, glynn=True)
hafnian_thewalrus = hafnian(A)
print(f"Hafnian__thewalrus: {hafnian_thewalrus}")
# Adjusted to match size of n
end_time2 = time.time()
time_taken2 = end_time2 - start_time2
print(f" Time_LHaf={time_taken2} seconds")
