import numpy as np
import numpy as np
from thewalrus.decompositions import takagi
import time
##### only thewalrus wont work, you will need to w
from sympy import symbols, expand, factorial2, Poly
from sympy import re


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
        #print(f"poly: {poly}")  # Output the constant term symbolically

    # Convert polynomial to dictionary
    p = Poly(poly, x)  # Polynomial with respect to all x variables
    #print(f"polynomial: {p}")  # Output the constant term symbolically
    terms = p.as_dict()
    #print(f"terms: {terms}")  # Output the constant term symbolically

    # Filter terms where sum of exponents is equal to n and each exponent is even
    valid_terms = {k: v for k, v in terms.items() if sum(k) == n and all(exp % 2 == 0 for exp in k)}

    # Initialize the result
    result = 0

    # Perform double factorial computation for each valid term
    for key, coeff in valid_terms.items():
        #print(f"key: {key}, coeff: {coeff}")  # Output the term (exponent tuple and coefficient)

        # Adjust exponents by subtracting 1
        adjusted_exponents = [exp - 1 for exp in key]
        #print(f"adjusted_exponents: {adjusted_exponents}")  # Output the adjusted exponents

        # Compute double factorial for adjusted exponents
        double_factorials = [factorial2(exp) if exp >= 0 else 1 for exp in adjusted_exponents]
        #print(f"double_factorials: {double_factorials}")  # Output the double factorials

        # Multiply the double factorials and coefficient
        term_value = coeff
        for df in double_factorials:
            term_value *= df

        # Add to the result
        result += term_value

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
    singular_values, U = takagi(A, svd_order=True, rtol=1e-16)

    # Step 4: Construct G using U and square roots of singular values
    G = U[:, :r] @ np.diag(np.sqrt(singular_values[:r]))

    # Step 5: Ensure G is real (handle small imaginary parts due to precision)
    G = np.real_if_close(G, tol=1e-10)

    # Step 6: Verify A ≈ G @ G.T
    if np.allclose(A, G @ G.T, atol=1e-8):
        np.savez("hafnianOnly_2.npz", A=A, G=G)  # Save matrices to a file
        break  # Exit loop if condition is satisfied

datatest = np.load("hafnianOnly_2.npz")
A = datatest['A']
G = datatest['G']
#print(f"G matrix: {G}")

# Generate polynomial and calculate the result
final_result = generate_hafnian(G)
print("Hafnian__code:", final_result)

# Output timing
print(f"Process completed in {time.time() - start_time:.2f} seconds.")

from thewalrus import hafnian
# Calculate loop Hafnian using thewalrus
start_time2 = time.time()
# thewalrus_haf = loop_hafnian(A = A_reconstructed, D=Mu, reps=None, glynn=True)
hafnian_thewalrus = hafnian(A)
print(f"Hafnian__thewalrus: {hafnian_thewalrus}")
# Adjusted to match size of n
end_time2 = time.time()
time_taken2 = end_time2 - start_time2
print(f" Time_LHaf={time_taken2} seconds")