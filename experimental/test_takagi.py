import numpy as np
from thewalrus.decompositions import takagi
import time
from sympy import symbols, expand, factorial2, Poly

def generate_polynomial_and_calculate(G, mu):

    # If n is odd, return 0 as specified
    if n % 2 != 0:
        return 0

    # Define symbolic variables x0, x1, ..., xr-1
    x = symbols(f"x0:{r}")
    poly = 1

    # Construct the multivariate polynomial
    for i in range(n):
        term = sum(G[i, j] * x[j] for j in range(r)) + mu[i]
        poly = expand(poly * term)

    # Convert polynomial to dictionary
    p = Poly(poly, x)  # Polynomial with respect to all x variables
    #print(f"polynomial: {p}")  # Output the constant term symbolically
    terms = p.as_dict()

    # Filter terms with even powers only
    even_terms = {k: v for k, v in terms.items() if all(exp % 2 == 0 for exp in k)}

    # Initialize the result
    result = 0

    # Perform double factorial computation for each term
    for key, coeff in even_terms.items():
        # Adjust exponents by subtracting 1
        adjusted_exponents = [exp - 1 for exp in key]

        # Compute double factorial for adjusted exponents
        double_factorials = [factorial2(exp) if exp >= 0 else 1 for exp in adjusted_exponents]

        # Multiply the double factorials and coefficient
        term_value = coeff
        for df in double_factorials:
            term_value *= df

        # Add to the result
        result += term_value

    return result


# Example usage
if __name__ == "__main__":
    # Define the matrix G and vector mu
    # Step 1: Generate a random symmetric matrix A

    """n = 4  # Size of the matrix
    r = 4
    # A = np.random.rand(n, n)  # Generate a random matrix
    A = np.array([[1, 2, 3, 4],
                  [2, 5, 6, 7],
                  [3, 6, 8, 9],
                  [4, 7, 9, 10]])
    A = (A + A.T) / 2  # Make it symmetric"""


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
            np.savez("test_takagi_2.npz", A=A, G=G)  # Save matrices to a file
            break  # Exit loop if condition is satisfied

    # Output timing
    print(f"Process completed in {time.time() - start_time:.2f} seconds.")

    mu = np.ones(n)  # Adjusted to match size of n

    # Generate polynomial and calculate the result
    final_result = generate_polynomial_and_calculate(G, mu)
    print("Low rank loop Hafnian:", final_result)

    end_time = time.time()  # Start timing

    timetaken = end_time - start_time
    print("Time taken :", timetaken)



    from thewalrus import loop_hafnian

    # Calculate loop Hafnian using thewalrus
    start_time2 = time.time()
    # thewalrus_haf = loop_hafnian(A = A_reconstructed, D=Mu, reps=None, glynn=True)
    loophafnian = loop_hafnian(A=A, D=mu, reps=None, glynn=True)
    print(f"Loop Hafnian (thewalrus): {loophafnian}")
    # Adjusted to match size of n
    end_time2 = time.time()
    time_taken2 = end_time2 - start_time2
    print(f" Time_LHaf={time_taken2} seconds")

 # Calculate loop Hafnian using thewalrus
    from thewalrus import hafnian
    # thewalrus_haf = loop_hafnian(A = A_reconstructed, D=Mu, reps=None, glynn=True)
    hafn = hafnian(A)
    print(f"Hafnian (thewalrus): {hafn}")
    # Adjusted to match size of n

"""
# Example usage:
A = np.array([[1, 2, 3, 4],
                  [2, 5, 6, 7],
                  [3, 6, 8, 9],
                  [4, 7, 9, 10]])
r = 4  # Desired rank
G = takagi(A, r)
print("Low-rank matrix G from Takagi decomposition:")
print(G)
"""