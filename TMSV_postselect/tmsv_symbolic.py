import sympy as sp
from typing import Optional, Dict
from sympy.matrices.common import NonInvertibleMatrixError
from thewalrus import hafnian as hafnian

# --------------------------
# Define sigma and M outside
# --------------------------
r = sp.symbols('r', real=True)
I2 = sp.eye(2)

sigma_x = sp.Matrix([[0, 1],
                     [1, 0]])  # Pauli X
# σ = [[C, S], [S, C]]
C = sp.cosh(2 * r) * I2
S = sp.sinh(2 * r) * sigma_x  # <-- replaced I2 with Pauli X

sigma = sp.Matrix([[C, S],
                   [S, C]])

# M = diag(I2, -I2)
M = sp.Matrix([[I2, sp.zeros(2)],
               [sp.zeros(2), -I2]])


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
        ab = sigma_eff[:block_dim, block_dim:2 * block_dim]
        bb = sigma_eff[block_dim:2 * block_dim, block_dim:2 * block_dim]

    # Identify the coefficient blocks (choose + sign for cross term)
    M_alpha = sp.simplify(aa / 2)
    M_beta = sp.simplify(bb / 2)
    M_alpha_beta = sp.simplify(ab / 2)

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

    Returns: Mat, haf_value, normalization, result
    """
    # Get exponent from compute_B_C_gamma
    B, C_row, exponent_gamma = compute_B_C_gamma(M_alpha_beta, M_beta, alpha, zeta)

    # Build total exponent = ζ†ζ - α†M_alpha α + exponent_gamma
    zeta_dag_zeta = sp.expand((zeta.H * zeta)[0, 0])
    alpha_dag_M_alpha_alpha = sp.expand(-(alpha.H * M_alpha * alpha)[0, 0])
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

    # Calculate hafnian(Mat) for [[A, K], [K^T, F]] structure
    # hafnian = sum over all perfect matchings
    # For 2-mode: hafnian = K[0,0]*K[1,1] + K[0,1]*K[1,0]
    m = len(alpha)  # Number of alpha modes
    K = Mat[:m, m:]  # Extract cross-coupling block
    haf_value = sp.simplify(K[0, 0] * K[1, 1] + K[0, 1] * K[1, 0])

    # Normalization = sqrt(det(M_beta))
    normalization = sp.sqrt(sp.simplify(M_beta.det()))

    # Final result = hafnian(Mat) / sqrt(det(M_beta))
    result = sp.simplify(haf_value / normalization)

    return Mat, haf_value, normalization, result


if __name__ == "__main__":
    M_alpha, M_beta, M_alpha_beta, *_ = derive_M_blocks_from_inputs(sigma, M)

    # Symbolic 2-vectors (alpha, zeta)
    a1, a2, z1, z2 = sp.symbols('a1 a2 z1 z2', complex=True)
    alpha_vec = sp.Matrix([a1, a2])
    zeta_vec = sp.Matrix([z1, z2])

    # Compute Mat and hafnian
    Mat, haf_value, normalization, result = compute_mat_and_hafnian(
        M_alpha, M_beta, M_alpha_beta, alpha_vec, zeta_vec
    )

    print("Mat =")
    sp.pprint(Mat)

    print("\nHafnian(Mat) =")
    sp.pprint(haf_value)

    print("\nNormalization sqrt(det(M_beta)) =")
    sp.pprint(normalization)

    print("\nResult: Hafnian(Mat) / sqrt(det(M_beta)) =")
    sp.pprint(result)

    # --- numeric sanity check ---
    vals = {
        r: 0.3,
        a1: 0.2+0.1j, a2: -0.15+0.05j,
        z1: 0.1-0.2j, z2: -0.05+0.1j,
    }
    print("\nNumeric result:")
    print(sp.N(result.subs(vals)))

# No top-level prints or side effects. This module exposes the following:
# - derive_M_blocks_from_inputs(sigma, M)
# - compute_B_C_gamma(M_alpha_beta, M_beta, alpha, zeta)
# - compute_mat_and_hafnian(M_alpha, M_beta, M_alpha_beta, alpha, zeta)  [symbolic Mat]
# - mat_and_hafnian_from_blocks(M_alpha, M_beta, M_alpha_beta, subs=None)  [Mat + thewalrus hafnian]

