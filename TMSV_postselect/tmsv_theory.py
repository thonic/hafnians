import sympy as sp
from typing import Optional, Dict
from sympy.matrices.common import NonInvertibleMatrixError
from thewalrus import hafnian as hafnian

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
    half = sp.Rational(1, 2)
    i = sp.I
    B = M_alpha_beta * alpha - i * half * zeta
    C_row = (alpha.H * M_alpha_beta.H) - i * half * zeta.H
    M_beta_inv = M_beta.inv()
    term1 = (C_row * M_beta_inv * B)[0]
    diff = C_row - B.H
    term2 = sp.Rational(1, 4) * (diff * M_beta_inv * diff.H)[0]
    exponent = sp.simplify(term1 + term2)
    return B, C_row, sp.exp(exponent)

# ---------------------------------------
# Function: build Mat for Z=[zeta; zeta*] and its hafnian
# ---------------------------------------
def compute_mat_and_hafnian(M_alpha: sp.Matrix,
                            M_beta: sp.Matrix,
                            M_alpha_beta: sp.Matrix,
                            alpha: sp.Matrix,
                            zeta: sp.Matrix):
    """
    Build the quadratic exponent matrix for Z = [alpha; zeta]:
        1/2 Z^† Mat_herm Z

    with blocks
        A = 2( M_alpha_beta^† M_beta^{-1} M_alpha_beta - M_alpha ),
        F = 2 I2,
        G = -i M_alpha_beta^† M_beta^{-1}.

    We also return the complex-symmetric Mat_sym = [[0,K],[K^T,0]] with
        K = M_beta^{-1} M_alpha_beta
    which is suitable for hafnian(Mat_sym).
    """
    I2 = sp.eye(2)
    D = M_beta.inv()

    A = sp.simplify(2 * (M_alpha_beta.H * D * M_alpha_beta - M_alpha))
    F = 2 * I2
    G = sp.simplify(sp.I * M_alpha_beta.H * D)

    Mat_herm = sp.BlockMatrix([[A, G], [G.H, F]]).as_explicit()

    # For 1/2 Z^T Mat Z form, off-diagonal should be i * K with no minus sign
    K = sp.simplify(sp.I * D * M_alpha_beta)
    Mat_sym = sp.BlockMatrix([[sp.zeros(2), K], [K.T, sp.zeros(2)]]).as_explicit()

    haf = sp.Symbol("hafnian(Mat_sym)")
    return Mat_herm, Mat_sym, haf, K


def mat_and_hafnian_from_blocks(M_alpha: sp.Matrix,
                                M_beta: sp.Matrix,
                                M_alpha_beta: sp.Matrix,
                                subs: Optional[Dict] = None):
    """Build Mat for Z=[alpha; zeta] and hafnian input without eliminating alpha.

    Returns:
      Mat_herm = [[A, G], [G^†, F]] with
        A = 2( M_ab^† M_beta^{-1} M_ab - M_alpha ),
        F = 2 I2,
        G = -i M_ab^† M_beta^{-1},
      Mat_sym = [[0, K], [K^T, 0]] with K = M_beta^{-1} M_ab,
      haf_value via The Walrus if `subs` yields numeric matrix, and K.
    """
    def _ensure_2x2(mat):
        return mat[0, 0] if isinstance(mat[0, 0], sp.MatrixBase) else mat

    I2_loc = sp.eye(2)
    Ma = _ensure_2x2(M_alpha)
    Mb = _ensure_2x2(M_beta)
    Mab = _ensure_2x2(M_alpha_beta)

    D = Mb.inv()
    A = sp.simplify(2 * (Mab.H * D * Mab - Ma))
    F = 2 * I2_loc
    G = sp.simplify(-sp.I * Mab.H * D)
    Mat_herm = sp.BlockMatrix([[A, G], [G.H, F]]).as_explicit()

    K = sp.simplify(D * Mab)
    Mat_sym = sp.BlockMatrix([[sp.zeros(2), K], [K.T, sp.zeros(2)]]).as_explicit()

    haf_value = None
    try:
        haf_value = hafnian_with_walrus(Mat_sym, subs=subs)
    except Exception:
        pass

    return Mat_herm, Mat_sym, haf_value, K


def hafnian_with_walrus(Mat_sym: sp.Matrix, subs: Optional[Dict] = None):
    """Compute the hafnian using The Walrus for a complex-symmetric matrix.

    Args:
        Mat_sym: 2n x 2n complex-symmetric sympy Matrix (e.g., [[0,K],[K^T,0]]).
        subs:    Optional dict of substitutions to evaluate symbols numerically
                 (e.g., {r: 0.3}). Required if Mat_sym contains symbols.

    Returns:
        Complex numeric hafnian as returned by The Walrus.
    """
    M_eval = Mat_sym
    if subs:
        M_eval = M_eval.subs(subs)
    M_eval = sp.N(M_eval)

    if len(M_eval.free_symbols) != 0:
        raise ValueError("Mat_sym contains symbols; provide numeric substitutions via `subs`.")

    import numpy as np
    M_np = np.array(M_eval.tolist(), dtype=complex)
    return hafnian(M_np)


def _hafnian_4x4_from_K_symbolic(K: sp.Matrix) -> sp.Expr:
    """Symbolic hafnian for Mat_sym = [[0,K],[K^T,0]] (4x4 case).

    haf(Mat_sym) = K00*K11 + K01*K10
    """
    return sp.simplify(K[0, 0] * K[1, 1] + K[0, 1] * K[1, 0])


def hafnian_symbolic_in_CS() -> sp.Expr:
    """Symbolic hafnian for Mat built from generic scalars C and S.

    Uses the closed-form K for the scalar-block model with no simplification
    by the identity C^2 - S^2 = 1.
    """
    C, S = sp.symbols('C S', complex=True)
    k = 1 - (S**2) / (2 * (C + 1) * (S**2 - (C**2 - 1)))
    return (k**2).expand()


def simple_mat_and_hafnian_from_blocks(M_alpha: sp.Matrix,
                                       M_beta: sp.Matrix,
                                       M_alpha_beta: sp.Matrix):
    """Simplest Mat construction that yields nontrivial r-dependence.

    Define K = M_beta^{-1} M_alpha_beta. Then set Mat_sym = [[0,K],[K^T,0]].
    Returns (Mat_sym, haf_expr, K) with haf_expr symbolic via 4x4 closed form.
    """
    def _ensure_2x2(mat):
        return mat[0, 0] if isinstance(mat[0, 0], sp.MatrixBase) else mat

    Ma = _ensure_2x2(M_alpha)
    Mb = _ensure_2x2(M_beta)
    Mab = _ensure_2x2(M_alpha_beta)

    K = sp.simplify(sp.I * Mb.inv() * Mab)
    Mat_sym = sp.BlockMatrix([[sp.zeros(2), K],
                              [K.T,        sp.zeros(2)]]).as_explicit()
    haf_expr = _hafnian_4x4_from_K_symbolic(K)
    return Mat_sym, haf_expr, K


def det_M_beta(M_beta: sp.Matrix) -> sp.Expr:
    """Determinant of the 2x2 block M_beta (supports block or explicit)."""
    Mb = M_beta[0, 0] if isinstance(M_beta[0, 0], sp.MatrixBase) else M_beta
    return sp.simplify(Mb.det())


if __name__ == "__main__":
    # Build blocks, construct Mat_herm for 1/2 Z^† Mat Z (Z=[alpha; zeta]),
    # and print hafnian(Mat_sym)/det(M_beta) symbolically.
    M_alpha, M_beta, M_alpha_beta, *_ = derive_M_blocks_from_inputs(sigma, M)

    # symbolic 2-vectors (alpha, zeta)
    a1, a2, z1, z2 = sp.symbols('a1 a2 z1 z2', complex=True)
    alpha_vec = sp.Matrix([a1, a2])
    zeta_vec = sp.Matrix([z1, z2])

    Mat_herm, Mat_sym, _, K = compute_mat_and_hafnian(
        M_alpha, M_beta, M_alpha_beta, alpha_vec, zeta_vec
    )

    # For 1/2 Z^T Mat Z, use complex-symmetric Mat (Mat_sym)
    print("Mat =")
    sp.pprint(Mat_sym)

    det_beta = det_M_beta(M_beta)
    haf_expr = _hafnian_4x4_from_K_symbolic(K)
    print("\nHafnian(Mat)/det(M_beta) =")
    sp.pprint(sp.simplify(haf_expr / sp.sqrt(det_beta)))

# No top-level prints or side effects. This module exposes the following:
# - derive_M_blocks_from_inputs(sigma, M)
# - compute_B_C_gamma(M_alpha_beta, M_beta, alpha, zeta)
# - compute_mat_and_hafnian(M_alpha, M_beta, M_alpha_beta, alpha, zeta)  [symbolic Mat]
# - mat_and_hafnian_from_blocks(M_alpha, M_beta, M_alpha_beta, subs=None)  [Mat + thewalrus hafnian]
