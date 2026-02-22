import sympy as sp
#import time


# ============================================================
# basic vector formation (alpha, beta structure will follow def complex_vector(name, n): to form desired vector ordering)
# ============================================================

def complex_vector(name, n):
    vec = sp.Matrix(sp.symbols(f'{name}1:{n+1}'))
    vec_star = sp.Matrix(sp.symbols(f'{name}1:{n+1}_star'))
    return vec, vec_star


def phase_space_vector(n, m):
    alpha, alpha_star = complex_vector("alpha", n)
    beta, beta_star = complex_vector("beta", m)
    v = sp.Matrix.vstack(alpha, beta, alpha_star, beta_star)
    return v, alpha, beta, alpha_star, beta_star


def dagger_v(alpha, beta, alpha_star, beta_star):
    return sp.Matrix(
        list(alpha_star) + list(beta_star) + list(alpha) + list(beta)
    ).T


# ============================================================
# J matrix
# ============================================================

def symplectic_J(n, m):
    I_n = sp.eye(n)
    I_m = sp.eye(m)

    Z_nn = sp.zeros(n)
    Z_mm = sp.zeros(m)
    Z_nm = sp.zeros(n, m)
    Z_mn = sp.zeros(m, n)

    return sp.BlockMatrix([
        [Z_nn, Z_nm,  I_n,  Z_nm],
        [Z_mn, Z_mm,  Z_mn, I_m ],
        [-I_n, Z_nm,  Z_nn, Z_nm],
        [Z_mn, -I_m,  Z_mn, Z_mm]
    ]).as_explicit()

# ============================================================
# two_mode_squeezed_covariance matrix
# ============================================================
def two_mode_squeezed_covariance(n, m, r_val):

    c2 = sp.cosh(2*r_val) / 2
    s2 = sp.sinh(2*r_val) / 2

    k = min(n, m)

    I_n = sp.eye(n)
    I_m = sp.eye(m)

    Z_nn = sp.zeros(n)
    Z_mm = sp.zeros(m)
    Z_nm = sp.zeros(n, m)
    Z_mn = sp.zeros(m, n)

    C = sp.zeros(n, m)
    for i in range(k):
        C[i, i] = 1

    return sp.BlockMatrix([
        [c2*I_n,   Z_nm,     Z_nn,     s2*C     ],
        [Z_mn,     c2*I_m,   s2*C.T,   Z_mm     ],
        [Z_nn,     s2*C,     c2*I_n,   Z_nm     ],
        [s2*C.T,   Z_mm,     Z_mn,     c2*I_m   ]
    ]).as_explicit()


# ============================================================
# displacement vector (numeric substitution early)
# ============================================================

def displacement_subs(n, m,
                      d_alpha_vals=0,
                      d_beta_vals=0):

    d_alpha_syms = sp.symbols(f'd_alpha1:{n+1}')
    d_beta_syms  = sp.symbols(f'd_beta1:{m+1}')

    subs = {}

    # -------- alpha --------
    if isinstance(d_alpha_vals, (int, float, complex, sp.Number)):
        for s in d_alpha_syms:
            subs[s] = d_alpha_vals
    else:
        if len(d_alpha_vals) != n:
            raise ValueError("Length of d_alpha_vals must be n")
        for s, val in zip(d_alpha_syms, d_alpha_vals):
            subs[s] = val

    # -------- beta --------
    if isinstance(d_beta_vals, (int, float, complex, sp.Number)):
        for s in d_beta_syms:
            subs[s] = d_beta_vals
    else:
        if len(d_beta_vals) != m:
            raise ValueError("Length of d_beta_vals must be m")
        for s, val in zip(d_beta_syms, d_beta_vals):
            subs[s] = val

    return subs



# ============================================================
# Extract M, B, C matrices and vectors to get * gamma *
# ============================================================

def extract_matrix_MBC(u_total, alpha, alpha_star):

    n = alpha.rows

    M = sp.Matrix(n, n, lambda i, j:
                  -sp.diff(u_total, alpha[j], alpha_star[i]))

    u1 = u_total + (alpha_star.T * M * alpha)[0]

    subs0 = {alpha[k]: 0 for k in range(n)}
    subs0.update({alpha_star[k]: 0 for k in range(n)})

    B = sp.Matrix([
        sp.diff(u1, alpha[j]).subs(subs0)
        for j in range(n)
    ])

    C = sp.Matrix([
        sp.diff(u1, alpha_star[i]).subs(subs0)
        for i in range(n)
    ])

    u_rest = u1 - (B.T * alpha)[0] - (alpha_star.T * C)[0]

    return M, B, C, u_rest


# ============================================================
# gamma
# ============================================================

def gamma_from_MBC(M, B, C):
    detM = M.det()
    MinvC = M.LUsolve(C)
    return sp.exp((B.T * MinvC)[0]) / detM


# ============================================================
# FULL NUMERIC PIPELINE
# ============================================================

def evaluate_numeric(n, m, j,
                     d_alpha_vals=0,
                     d_beta_vals=0):

    import time
    t0 = time.perf_counter()

    # --------------------------------------------------------
    # 1) Fix squeezing parameter numerically
    # --------------------------------------------------------
    r_val = sp.asinh(1)

    # --------------------------------------------------------
    # 2) Build phase-space variables
    # --------------------------------------------------------
    v, alpha, beta, alpha_star, beta_star = phase_space_vector(n, m)
    v_dag = dagger_v(alpha, beta, alpha_star, beta_star)

    # --------------------------------------------------------
    # 3) Build covariance and J
    # --------------------------------------------------------
    Sigma = two_mode_squeezed_covariance(n, m, r_val)
    J = symplectic_J(n, m)

    # --------------------------------------------------------
    # 4) Build displacement symbols properly
    #    Star components are true conjugates
    # --------------------------------------------------------
    d_alpha_syms = sp.symbols(f'd_alpha1:{n+1}')
    d_beta_syms  = sp.symbols(f'd_beta1:{m+1}')

    D = sp.Matrix(
        list(d_alpha_syms)
        + list(d_beta_syms)
        + [sp.conjugate(x) for x in d_alpha_syms]
        + [sp.conjugate(x) for x in d_beta_syms]
    )

    # --------------------------------------------------------
    # 5) Build ζ variables
    # --------------------------------------------------------
    zeta, zeta_star = complex_vector("zeta", n)

    # --------------------------------------------------------
    # 6) Construct full exponent EXACTLY as symbolic version
    # --------------------------------------------------------
    quad = -sp.Rational(1, 2) * (v_dag * Sigma * v)[0]

    ordering = (
        -sp.Rational(1, 2) * (alpha_star.T * alpha)[0]
        + sp.Rational(1, 2) * (beta_star.T * beta)[0]
    )

    zeta_linear = (alpha_star.T * zeta)[0] - (zeta_star.T * alpha)[0]
    zeta_sq = (zeta_star.T * zeta)[0]

    linear_disp = (v.T * (J * D))[0]

    u_total = quad + ordering + zeta_linear + zeta_sq + linear_disp

    # --------------------------------------------------------
    # 7) Substitute displacement NUMERICALLY
    # --------------------------------------------------------
    disp_subs = {}

    # alpha displacement
    if isinstance(d_alpha_vals, (int, float, complex, sp.Number)):
        for s in d_alpha_syms:
            disp_subs[s] = d_alpha_vals
    else:
        if len(d_alpha_vals) != n:
            raise ValueError("Length of d_alpha_vals must equal n")
        for s, val in zip(d_alpha_syms, d_alpha_vals):
            disp_subs[s] = val

    # beta displacement
    if isinstance(d_beta_vals, (int, float, complex, sp.Number)):
        for s in d_beta_syms:
            disp_subs[s] = d_beta_vals
    else:
        if len(d_beta_vals) != m:
            raise ValueError("Length of d_beta_vals must equal m")
        for s, val in zip(d_beta_syms, d_beta_vals):
            disp_subs[s] = val

    u_total = u_total.subs(disp_subs)

    # --------------------------------------------------------
    # 8) Eliminate alpha via derivative definition
    # --------------------------------------------------------
    M, B, C, u_rest = extract_matrix_MBC(
        u_total, alpha, alpha_star
    )

    gamma = gamma_from_MBC(M, B, C)

    effective_expr = sp.exp(u_rest) * gamma

    # --------------------------------------------------------
    # 9) Apply derivatives:
    #    ∂_beta ∂_beta*
    #    ∂^j_zeta ∂^j_zeta*
    # --------------------------------------------------------
    vars_to_diff = list(beta) + list(beta_star)

    for i in range(n):
        vars_to_diff += [zeta[i]] * j
        vars_to_diff += [zeta_star[i]] * j

    result = sp.diff(effective_expr, *vars_to_diff)

    # Apply (-1)^m due to d/d(-beta*)
    result *= (-1) ** len(beta_star)

    # --------------------------------------------------------
    # 10) Set β, ζ = 0
    # --------------------------------------------------------
    subs_zero = {
        x: 0 for x in
        list(beta) + list(beta_star)
        + list(zeta) + list(zeta_star)
    }

    result = result.subs(subs_zero)

    t1 = time.perf_counter()

    result_numeric = sp.N(result)

    print("Result =", result_numeric)
    print("Runtime =", t1 - t0)

    return result_numeric



# ============================================================
# substitute no of alpha,beta modes as (n, m) and no of photons in projection operator as j..
# with displacement vector component values
# ============================================================

if __name__ == "__main__":

    evaluate_numeric(
        n=1,
        m=1,
        j=1,
        d_alpha_vals=1.0,
        d_beta_vals=[1.0]
    )


# ============================================================
#
#         *  Algorithm  *
#
# ============================================================
#

# evaluate_numeric(n, m, j, d_alpha_vals, d_beta_vals)
#
# Computes  <j | ρ | j>  using the exact generating-function
# method derived symbolically.
#
# Exponent constructed as:
#
#   u_total =
#     -1/2 v† Σ v
#     + vᵀ J D
#     -1/2 α*α
#     +1/2 β*β
#     + α*ζ - ζ*α
#     + ζ*ζ
#
# Steps:
#   1) Substitute numeric displacement.
#   2) Extract α-quadratic form:
#          -α* M α + Bᵀα + α*C
#      via derivative definitions.
#   3) Perform Gaussian integration:
#          γ = exp(Bᵀ M⁻¹ C) / det(M)
#   4) Differentiate remaining expression w.r.t.
#          β, (-β*), ζ^j, ζ*^j
#   5) Set β=ζ=0 → numeric result.
#
# Algorithm matches the symbolic derivation exactly.
# ============================================================
