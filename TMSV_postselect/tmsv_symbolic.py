import sympy as sp
from typing import Optional, Dict
from sympy.matrices.common import NonInvertibleMatrixError
from thewalrus import hafnian as hafnian_thewalrus


r = sp.symbols('r', real=True)
I2 = sp.eye(2)

C = sp.cosh(2*r) * I2
S = sp.sinh(2*r) * I2

# ---- Define multiple 2x2 unitaries ----
# Identity unitary (no transformation)
u_I = sp.eye(2)

# Hadamard unitary
u_h = (1/sp.sqrt(2)) * sp.Matrix([[1, 1],
                                  [1, -1]])

# 50/50 Beam Splitter with phase (Hadamard + π/2 phase shift)
u_bs = (1/sp.sqrt(2)) * sp.Matrix([[1, sp.I],
                                   [sp.I, 1]])

# Asymmetric beam splitter (30 degree angle = π/6) - will give different result!
theta_val = sp.pi/6  # 30 degrees
u_asym = sp.Matrix([[sp.cos(theta_val), sp.sin(theta_val)],
                    [-sp.sin(theta_val), sp.cos(theta_val)]])



# ---- Select which unitary to use ----
use_unitary = True  # Set to False to skip transformation (original behavior)
u = u_h  # Options: u_I, u_h, u_bs (same as identity!), u_phase (same as identity!), u_swap, u_asym (complicated)

# ---- σ in v_p ordering (alpha, alpha*, beta, beta*): [[C, S], [S, C]] as 4x4 ----
sigma_p_block = sp.BlockMatrix([[C, S],
                                [S, C]])
sigma_p = sigma_p_block.as_explicit()

if use_unitary:
    # build big U in v-ordering (alpha, beta, alpha*, beta*) ----
    U_v = sp.BlockMatrix([[u,              sp.ZeroMatrix(2,2)],
                          [sp.ZeroMatrix(2,2),    u.conjugate()]])
    
    # permutation P: v_p = P * v  (maps (alpha, beta, alpha*, beta*) -> (alpha, alpha*, beta, beta*))
    P = sp.Matrix([[1, 0, 0, 0],
                   [0, 0, 1, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 1]])
    
    # permute the unitary into v_p ordering ----
    U_p = P * U_v.as_explicit() * P.T
    
    # active state transformation in v_p ordering: σ'_p = U_p^† σ_p U_p ----
    sigma = sp.simplify(U_p.H * sigma_p * U_p)
else:
    # Original behavior: no transformation
    sigma = sigma_p

# M = diag(I2, -I2)
M_block = sp.BlockMatrix([[I2,                sp.ZeroMatrix(2,2)],
                          [sp.ZeroMatrix(2,2), -I2              ]])
M = M_block.as_explicit()

# ---------------------------------------
# Function: symbolic hafnian calculator
# ---------------------------------------
def hafnian_symbolic(A):
    """
        Symbolic expression for hafnian
    """
    n = A.shape[0]
    
    # empty matrix
    if n == 0:
        return sp.Integer(1)
    
    # Hafnian is zero for odd-sized matrices
    if n % 2 != 0:
        return sp.Integer(0)
    
    # 2×2 matrix
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


# M blocks from σ and M

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

    
    if sigma_eff.shape == (2, 2) and isinstance(sigma_eff[0, 0], sp.MatrixBase):
        aa = sigma_eff[0, 0]
        ab = sigma_eff[0, 1]
        bb = sigma_eff[1, 1]
    else:
        aa = sigma_eff[:block_dim, :block_dim]
        ab = sigma_eff[:block_dim, block_dim:2*block_dim]
        bb = sigma_eff[block_dim:2*block_dim, block_dim:2*block_dim]

    
    M_alpha      = sp.simplify( aa / 2 )
    M_beta       = sp.simplify( bb / 2 )
    M_alpha_beta = sp.simplify( ab / 2 )

    return M_alpha, M_beta, M_alpha_beta, sigma_eff, sigma_in, M_in


# compute B, C, and gamma

def compute_B_C_gamma(M_alpha_beta: sp.Matrix,
                      M_beta: sp.Matrix,
                      alpha: sp.Matrix,
                      zeta: sp.Matrix):
   
    half = sp.Rational(1, 2)
    i = sp.I
    
    B = M_alpha_beta * alpha - i * half * zeta
    C_row = (alpha.H * M_alpha_beta.H) - i * half * zeta.H
    
    M_beta_inv = M_beta.inv()
    
    # gamma = exp[C M_β⁻¹ B + ¼ (C – B†) M_β⁻¹ (C – B†)†], exponent of this gamma we are computing
    term1 = (C_row * M_beta_inv * B)[0, 0]  # C M_β⁻¹ B
    diff = C_row - B.H  # (C – B†)
    term2 = sp.Rational(1, 4) * (diff * M_beta_inv * diff.H)[0, 0]  # ¼ (C – B†) M_β⁻¹ (C – B†)†
    
    exponent = sp.simplify(term1 + term2)
    return B, C_row, exponent


#  build Mat for Z=[zeta; zeta*] and compute its hafnian, where exp(exponent) is gamma

def compute_mat_and_hafnian(M_alpha: sp.Matrix,
                            M_beta: sp.Matrix,
                            M_alpha_beta: sp.Matrix,
                            alpha: sp.Matrix,
                            zeta: sp.Matrix):
    """
    Build Mat from exp(ζ†ζ) · exp(-alpha† M_alpha alpha) · exp(exponent),
    compute hafnian(Mat), and divide by sqrt(det(M_beta)) according to theory.
    
    Expands into: 1/2 Z^T Mat Z where Z = [alpha, alpha*, zeta, zeta*]
    
    Returns: Mat, haf_value, det_M_beta, result
    """
    
    B, C_row, exponent_gamma = compute_B_C_gamma(M_alpha_beta, M_beta, alpha, zeta)
    
    
    zeta_dag_zeta = sp.expand((zeta.H * zeta)[0,0])
    alpha_dag_M_alpha_alpha = sp.expand(-(alpha.H * M_alpha * alpha)[0,0])
    total_exponent = sp.expand(zeta_dag_zeta + alpha_dag_M_alpha_alpha + exponent_gamma)
    
    Z = sp.Matrix.vstack(alpha, zeta)
    n = len(Z)  # Derive dimension from Z
    
    Mat = sp.zeros(n, n)
    expanded_for_extraction = sp.expand(2 * total_exponent)
    
    for i in range(n):
        for j in range(n):
            coeff = expanded_for_extraction.coeff(Z[i] * sp.conjugate(Z[j]))
            Mat[i, j] = coeff if coeff is not None else 0
    
    Mat = sp.simplify(Mat)
    ## Calculate hafnian(Mat) for [[0, K], [K^T, F]] structure
    ## hafnian = sum over all perfect matchings
    ## For 2-mode: hafnian = K[0,0]*K[1,1] + K[0,1]*K[1,0]
    #m = len(alpha)  # Number of alpha modes
    #K = Mat[:m, m:]  # Extract cross-coupling block
    #haf_value = sp.simplify(K[0,0] * K[1,1] + K[0,1] * K[1,0])
    #det_M_beta = sp.simplify(M_beta.det())
    #result = sp.simplify(haf_value / sp.sqrt(det_M_beta))
    
    # Calculate hafnian(Mat) using general symbolic hafnian
    haf_value = hafnian_symbolic(Mat)
    det_M_beta = sp.simplify(M_beta.det()) 
    
    # Final result = hafnian(Mat) / sqrt(det(M_beta))
    result = sp.simplify(haf_value / sp.sqrt(det_M_beta))
    
    return Mat, haf_value, det_M_beta, result


if __name__ == "__main__":
    M_alpha, M_beta, M_alpha_beta, *_ = derive_M_blocks_from_inputs(sigma, M)

    # Use independent variables for 2-mode system (working version)
    a1, a2, z1, z2 = sp.symbols('a1 a2 z1 z2', complex=True)
    alpha_vec = sp.Matrix([a1, a2])
    zeta_vec = sp.Matrix([z1, z2])

    # Determine which configuration is being used
    if not use_unitary:
        config_name = "No transformation (original)"
    else:
        if u.equals(u_h):
            unitary_name = "Hadamard (u_h)"
        elif u.equals(u_I):
            unitary_name = "Identity (u_I)"
        elif u.equals(u_bs):
            unitary_name = "50/50 Beam Splitter with phase (u_bs) - SAME AS IDENTITY"
        elif u.equals(u_asym):
            unitary_name = "Asymmetric beam splitter (u_asym, θ=π/6=30°) - COMPLICATED"
        else:
            unitary_name = "Custom unitary"
        config_name = f"Unitary applied: {unitary_name}"
    
    print("="*70)
    print(f"2-mode system: alpha = [a1, a2], zeta = [z1, z2]")
    print(config_name)
    print("="*70)
    
    # Print transformed sigma for diagnostic
    print("\nTransformed σ =")
    sp.pprint(sigma)

    # Compute Mat and hafnian
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

