import sympy as sp
from typing import Optional, Dict
from sympy.matrices.common import NonInvertibleMatrixError
from thewalrus import hafnian as hafnian_thewalrus
import numpy as np

# --------------------------
# Define sigma and M outside
# --------------------------
r = sp.symbols('r', real=True)
I2 = sp.eye(2)

C = sp.cosh(2*r) * I2
S = sp.sinh(2*r) * I2

# σ = [[C, S], [S, C]]
sigma = sp.Matrix([[C, S],
                   [S, C]])

# M = diag(I2, -I2)
M = sp.Matrix([[I2,            sp.zeros(2)],
               [sp.zeros(2),  -I2        ]])

# ---------------------------------------
# Function: symbolic hafnian calculator
# ---------------------------------------
def hafnian_symbolic(A):
    """
    Compute hafnian of a symbolic matrix A using recursion.
    
    For a matrix A (n×n where n is even), hafnian is:
    haf(A) = sum over all perfect matchings of product of matched elements
    
    Args:
        A: sympy Matrix (n×n, symmetric)
    
    Returns:
        Symbolic expression for hafnian
    """
    n = A.shape[0]
    
    # Base case: empty matrix
    if n == 0:
        return sp.Integer(1)
    
    # Hafnian is zero for odd-sized matrices
    if n % 2 != 0:
        return sp.Integer(0)
    
    # Base case: 2×2 matrix
    if n == 2:
        return A[0, 1]
    
    # Recursive case: pair first element with each other element
    result = sp.Integer(0)
    
    for j in range(1, n):
        # Pair element 0 with element j
        # Create submatrix by removing rows/columns 0 and j
        indices = [i for i in range(1, n) if i != j]
        
        if len(indices) == 0:
            submatrix_haf = sp.Integer(1)
        else:
            submatrix = A[indices, indices]
            submatrix_haf = hafnian_symbolic(submatrix)
        
        result += A[0, j] * submatrix_haf
    
    return sp.simplify(result)

# ---------------------------------------
# Function: derive blocks from σ and M
# ---------------------------------------
def derive_M_blocks_from_inputs(sigma_in: sp.Matrix,
                                M_in: sp.Matrix,
                                block_dim: int = 2):
    """
    Given:
      - sigma_in: 4x4 block matrix [[C,S],[S,C]]
      - M_in:     4x4 correction matrix diag(I2, -I2)
    Returns:
      M_alpha, M_beta, M_alpha_beta (each 2x2), derived from:
         sigma_eff = sigma_in - M_in
         -1/2 v† sigma_eff v = -α† M_alpha α - β† M_beta β
                               + α† M_alpha_beta β + β† M_alpha_beta α
      so:
         M_alpha      = + (1/2) * sigma_eff(aa)
         M_beta       = + (1/2) * sigma_eff(bb)
         M_alpha_beta = - (1/2) * sigma_eff(ab)

      Additionally returns sigma_in and M_in for convenience.
    """
    # Effective matrix
    sigma_eff = sp.simplify(sigma_in - M_in)

    # Slice 2x2 blocks (aa, ab, bb). Support two representations:
    #  1) 4x4 explicit matrices (use slicing)
    #  2) 2x2 matrix of 2x2 blocks (use [i,j] block access)
    if sigma_eff.shape == (2, 2) and isinstance(sigma_eff[0, 0], sp.MatrixBase):
        aa = sigma_eff[0, 0]
        ab = sigma_eff[0, 1]
        bb = sigma_eff[1, 1]
    else:
        aa = sigma_eff[:block_dim, :block_dim]
        ab = sigma_eff[:block_dim, block_dim:2*block_dim]
        bb = sigma_eff[block_dim:2*block_dim, block_dim:2*block_dim]

    # Identify the coefficient blocks (choose + sign for cross term)
    M_alpha      = sp.simplify( aa / 2 )
    M_beta       = sp.simplify( bb / 2 )
    M_alpha_beta = sp.simplify( ab / 2 )

    return M_alpha, M_beta, M_alpha_beta, sigma_eff, sigma_in, M_in

# ---------------------------------------
# Function: compute B, C, and gamma
# ---------------------------------------
def compute_B_C_gamma(M_alpha_beta: sp.Matrix,
                      M_beta: sp.Matrix,
                      alpha: sp.Matrix,
                      zeta: sp.Matrix):
    """
    Compute B, C_row, and exponent for displacement terms.
    
    Returns: B, C_row, exponent
    """
    half = sp.Rational(1, 2)
    i = sp.I
    B = M_alpha_beta * alpha - i * half * zeta
    C_row = (alpha.H * M_alpha_beta.H) - i * half * zeta.H
    M_beta_inv = M_beta.inv()
    term1 = (C_row * M_beta_inv * B)[0]
    diff = C_row - B.H
    term2 = sp.Rational(1, 4) * (diff * M_beta_inv * diff.H)[0]
    exponent = sp.simplify(term1 + term2)
    return B, C_row, exponent

# ---------------------------------------
# Function: build Mat for Z=[zeta; zeta*] and its hafnian
# ---------------------------------------
def compute_mat_and_hafnian(M_alpha: sp.Matrix,
                            M_beta: sp.Matrix,
                            M_alpha_beta: sp.Matrix,
                            alpha: sp.Matrix,
                            zeta: sp.Matrix):
    """
    Build Mat from exp(ζ†ζ) · exp(-α†M_alpha α) · exp(exponent),
    compute hafnian(Mat), and normalize by sqrt(det(M_beta)).
    
    Expands into: 1/2 Z^T Mat Z where Z = [alpha, alpha*, zeta, zeta*]
    
    Returns: Mat, haf_value, det_M_beta, result
    """
    # Get exponent from compute_B_C_gamma
    B, C_row, exponent_gamma = compute_B_C_gamma(M_alpha_beta, M_beta, alpha, zeta)
    
    # Build total exponent = ζ†ζ - α†M_alpha α + exponent_gamma
    zeta_dag_zeta = sp.expand((zeta.H * zeta)[0,0])
    alpha_dag_M_alpha_alpha = sp.expand(-(alpha.H * M_alpha * alpha)[0,0])
    total_exponent = sp.expand(zeta_dag_zeta + alpha_dag_M_alpha_alpha + exponent_gamma)
    
    # Build Z = [alpha; zeta] - concatenate the mode vectors
    Z = sp.Matrix.vstack(alpha, zeta)
    n = len(Z)  # Derive dimension from Z
    
    # Extract Mat from total_exponent: total_exponent = 1/2 Z^T Mat Z
    # Mat[i,j] = 2 * coefficient of Z[i] * conjugate(Z[j]) in total_exponent
    Mat = sp.zeros(n, n)
    
    for i in range(n):
        for j in range(n):
            coeff = total_exponent.coeff(Z[i] * sp.conjugate(Z[j]))
            Mat[i, j] = 2 * coeff if coeff is not None else 0
    
    Mat = sp.simplify(Mat)
    
    # Calculate hafnian(Mat) using general symbolic hafnian
    haf_value = hafnian_symbolic(Mat)
    
    det_M_beta = sp.simplify(M_beta.det())
    
    # Final result = hafnian(Mat) / sqrt(det(M_beta))
    result = sp.simplify(haf_value / sp.sqrt(det_M_beta))
    
    return Mat, haf_value, det_M_beta, result


if __name__ == "__main__":
    M_alpha, M_beta, M_alpha_beta, *_ = derive_M_blocks_from_inputs(sigma, M)

    # Symbolic vectors (alpha, zeta) 
    # For 2-mode system: a1, a2 are independent; z1, z2 are independent
    # For clarity: a2 could be conjugate(a1), z2 could be conjugate(z1) for 1-mode
    # Current: 2-mode system with independent variables
    a1, a2, z1, z2 = sp.symbols('a1 a2 z1 z2', complex=True)
    alpha_vec = sp.Matrix([a1, a2])
    zeta_vec = sp.Matrix([z1, z2])
    
    # Alternative for 1-mode with explicit conjugates (uncomment to use):
    # a1, z1 = sp.symbols('a1 z1', complex=True)
    # alpha_vec = sp.Matrix([a1, sp.conjugate(a1)])
    # zeta_vec = sp.Matrix([z1, sp.conjugate(z1)])

    # ==============================
    # SYMBOLIC RESULTS
    # ==============================
    print("="*70)
    print("SYMBOLIC CALCULATION")
    print("="*70)
    
    Mat, haf_value, det_M_beta, result = compute_mat_and_hafnian(
        M_alpha, M_beta, M_alpha_beta, alpha_vec, zeta_vec
    )

    print("\nMat =")
    sp.pprint(Mat)

    print("\nHafnian(Mat) =")
    sp.pprint(haf_value)
    
    print("\ndet_M_beta =")
    sp.pprint(det_M_beta)
    
    print("\nResult: Hafnian(Mat) / sqrt(det(M_beta)) =")
    sp.pprint(result)

    # ==============================
    # NUMERIC RESULTS (r = 0.5)
    # ==============================
    print("\n" + "="*70)
    print("NUMERIC CALCULATION (r = 0.5)")
    print("="*70)
    
    # Substitution value (only r is needed, Mat doesn't depend on a1, a2, z1, z2)
    r_val = 0.5
    subs_dict = {r: r_val}
    
    # Evaluate Mat numerically
    Mat_numeric = Mat.subs(subs_dict).evalf()
    Mat_np = np.array(Mat_numeric.tolist(), dtype=complex)
    
    # Hafnian using thewalrus
    haf_numeric = hafnian_thewalrus(Mat_np)
    
    # det_M_beta (need sqrt for normalization)
    det_M_beta_numeric = complex(det_M_beta.subs(subs_dict).evalf())
    norm_numeric = complex(sp.sqrt(det_M_beta.subs(subs_dict)).evalf())
    
    # Final result
    result_numeric = haf_numeric / norm_numeric
    
    print(f"\nMat (numeric) =")
    print(Mat_np)
    
    print(f"\nHafnian(Mat) [thewalrus] = {haf_numeric}")
    print(f"det_M_beta = {det_M_beta_numeric}")
    print(f"sqrt(det_M_beta) = {norm_numeric}")
    print(f"\nResult: Hafnian/sqrt(det(M_beta)) = {result_numeric}")
    
    # Verify symbolic matches numeric
    result_symbolic_eval = complex(result.subs(subs_dict).evalf())
    print(f"\nVerification: Symbolic result evaluated = {result_symbolic_eval}")
    print(f"Match: {np.isclose(result_numeric, result_symbolic_eval)}")

# This module is self-contained and exposes the following functions:
# - hafnian_symbolic(A) - General symbolic hafnian calculator
# - derive_M_blocks_from_inputs(sigma, M) - Extract M_alpha, M_beta, M_alpha_beta
# - compute_B_C_gamma(M_alpha_beta, M_beta, alpha, zeta) - Compute displacement terms
# - compute_mat_and_hafnian(M_alpha, M_beta, M_alpha_beta, alpha, zeta) - Main calculation

