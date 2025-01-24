import numpy as np
from thewalrus.decompositions import takagi
from thewalrus import hafnian
import time
from functools import lru_cache


def generate_hafnian_no_sympy_bitmask(G):
    """
    Compute the hafnian of A = G G^T using a bitmask-based DP recursion.
    This typically outperforms explicit enumeration.
    """
    n, r = G.shape
    if n % 2 != 0:
        return 0.0

    # Precompute A = G G^T
    A = G @ G.T

    @lru_cache(None)
    def haf_bitmask(mask):
        # If mask == 0, it means no indices remain, so the matching is complete
        if mask == 0:
            return 1.0

        # Extract the lowest set bit (lowest index in the subset)
        i = (mask & -mask).bit_length() - 1

        res = 0.0
        # Remove i from the subset
        mask_without_i = mask ^ (1 << i)

        # Try pairing i with each j in mask_without_i
        sub = mask_without_i
        while sub != 0:
            # Extract the lowest set bit from sub
            j = (sub & -sub).bit_length() - 1
            # Remove j from sub
            sub = sub ^ (1 << j)

            # Add A[i,j] times the recursion with i and j removed
            res += A[i, j] * haf_bitmask(mask_without_i ^ (1 << j))

        return res

    # Initially, all n indices (bits) are set: mask = (1 << n) - 1
    full_mask = (1 << n) - 1
    return haf_bitmask(full_mask)


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
my_haf = generate_hafnian_no_sympy_bitmask(G)
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
