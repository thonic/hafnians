import numpy as np
from collections import defaultdict
from thewalrus.decompositions import takagi
import time
import cProfile, pstats, io
from pstats import SortKey

# Double factorial computation

# Initialize the profiler
"""pr = cProfile.Profile()
pr.enable()
import numpy as np
import time
from collections import defaultdict"""


# Function to generate symmetric matrix with exact rank r
def generate_random_symmetric_matrix(n, r):
    """Generate a symmetric matrix of size n x n with rank r."""
    A = np.random.randn(n, n)
    A = (A + A.T) / 2  # Make symmetric
    U, S, Vt = np.linalg.svd(A)
    S[r:] = 0  # Set singular values beyond rank r to zero
    return U @ np.diag(S) @ Vt  # Reconstruct A

# Parameters
n = 20  # Matrix size
r = 10  # Desired rank
print(f"n, r: {n, r}")
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

#print(f"A matrix: {A}")
#print(f"G matrix: {G}")

# Double factorial function optimized
def double_factorial(n):
    return 1 if n <= 0 else np.prod(range(n, 0, -2))
# Vectorized function to compute the weighted sum
def compute_weighted_sum(even_terms):
    return sum(
        coeff * np.prod([double_factorial(p - 1) for p in power]) for power, coeff in even_terms.items()
    )
#Step 2: Collect terms
terms = {(): 1}  # Initialize with neutral term
for row in G:
    new_terms = defaultdict(int)
    for power, coeff in terms.items():
        for i, g in enumerate(row):
            # Ensure power is extended to length r
            new_power = tuple((power[j] if j < len(power) else 0) + (1 if j == i else 0) for j in range(r))
            new_terms[new_power] += coeff * g
    terms = new_terms

# Step 3: Filter even powers
even_terms = {k: v for k, v in terms.items() if all(p % 2 == 0 for p in k)}

# Step 4: Compute the weighted sum
weighted_sum = compute_weighted_sum(even_terms)


# Final result
print("\nStep 4: Final Computation:")
print(f"The weighted sum based on filtered terms is: {weighted_sum}")

# Measure time taken
end_time = time.time()
timetaken = end_time - start_time
print(f"Time taken: {timetaken:.4f} seconds")



"""################## End profiling
pr.disable()
s = io.StringIO()
sortby = SortKey.CUMULATIVE #When Sorted by CUMULATIVE Time (SortKey.CUMULATIVE):
#The first line will show the function that has spent the most cumulative time (i.e., the total time spent in that function and all functions it called).

ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()     # printing all lines
ps.print_stats(10)    # You can print only the top N=10 lines of the result.
# This helps in focusing on the most time-consuming functions.
print(s.getvalue())"""

#########################################################################

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

