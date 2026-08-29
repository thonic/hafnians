# Copyright 2019 Xanadu Quantum Technologies Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
Algorithms for hafnians of low-rank matrices.
"""

from functools import lru_cache
from itertools import product
import numpy as np
import time
from sympy import symbols, expand
from scipy.special import factorial2
from thewalrus.decompositions import takagi


@lru_cache(maxsize=1000000)
def partitions(r, n):
    r"""Returns a list of lists with the r-partitions of the integer :math:`n`, i.e. the r-tuples of non-negative
    integers such that their sum is precisely :math:`n`.

    Note that there are :math:`n + r - 1 \choose r-1` such partitions.

    Args:
        r (int): number of partitions.
        n (int): integer to be partitioned.

    Returns:
        list: r-partitions of n.
    """
    if r == 1:
        return [[n]]

    new_combos = []
    for first_val in range(n + 1):
        rest = partitions(r - 1, n - first_val)
        new = [p[0] + p[1] for p in product([[first_val]], rest)]
        new_combos += new
    return new_combos


def low_rank_hafnian(G):
    r"""Returns the hafnian of the low rank matrix :math:`\bm{A} = \bm{G} \bm{G}^T` where :math:`\bm{G}` is rectangular of size
    :math:`n \times r`  with :math:`r \leq n`.

    Note that the rank of :math:`\bm{A}` is precisely :math:`r`.

    The hafnian is calculated using the algorithm described in Appendix C of
    *A faster hafnian formula for complex matrices and its benchmarking on a supercomputer*,
    :cite:`bjorklund2018faster`.

    Args:
        G (array): factorization of the low rank matrix A = G @ G.T.

    Returns:
        (complex): hafnian of A.
    """
    n, r = G.shape
    if n % 2 != 0:
        return 0
    if r == 1:
        return factorial2(n - 1) * np.prod(G)

    # Create symbolic variables
    x = symbols("x0:" + str(r))

    # Start with the polynomial
    poly = 1
    for k in range(n):
        term = 0
        for j in range(r):
            term += G[k, j] * x[j]
        poly = expand(poly * term)
    # print(f"poly: {poly}")  # Output the constant term symbolically

    # Generate the r-partitions
    comb = partitions(r, n // 2)
    haf_val = 0.0

    # Loop over all partitions and compute the hafnian
    for c in comb:
        monomial = 1
        facts = 1
        # Loop over each element of the partition `c`
        for i, pi in enumerate(c):
            # Construct the monomial for the partition `c`
            monomial *= x[i] ** (2 * pi)
            # print(f"monomial so far: {monomial}")

            # Correct factorial calculation, now calculating factorials for non-zero pi values
            if pi > 0:
                facts *= factorial2(2 * pi - 1)
            else:
                facts *= 1  # No factorial is needed for pi == 0

        # After constructing the monomial, extract the coefficient
        coeff = poly.coeff(monomial)
        # print(f"Extracted coefficient for {monomial}: {coeff}")

        # Add the weighted coefficient to the hafnian value
        haf_val += complex(coeff * facts)

    return haf_val


# Example usage:
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

start_time1 = time.time()

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
        np.savez("hafnian_lowrank_thewalrus.npz", A=A, G=G)  # Save matrices to a file
        break  # Exit loop if condition is satisfied

datatest = np.load("hafnian_lowrank_thewalrus.npz")
A = datatest["A"]
G = datatest["G"]

# Generate polynomial and calculate the result
final_result = low_rank_hafnian(G)
print("Hafnian__code:", final_result)
end_time1 = time.time()
time_taken1 = end_time1 - start_time1
print(f" Time_LHaf={time_taken1} seconds")

# Compare with thewalrus
from thewalrus import hafnian

start_time2 = time.time()
hafnian_thewalrus = hafnian(A)
print(f"Hafnian__thewalrus: {hafnian_thewalrus}")
end_time2 = time.time()
time_taken2 = end_time2 - start_time2
print(f" Time_LHaf={time_taken2} seconds")
