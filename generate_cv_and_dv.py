from scipy.linalg import block_diag

from generate_displacements import *


def delete_vec(m, v):
    x = m
    for j in range(len(v) - 1, -1, -1):
        if v[j] == 1:
            x = np.delete(x, j)
    return x


def delete_cv(m, v):
    x = m
    for j in range(len(v) - 1, -1, -1):
        if v[j] == 1:
            x = np.delete(x, j, 0)
            x = np.delete(x, j, 1)
    return x


def reshape_mat(m, MK_dict):
    new_vec = np.zeros(len(MK_dict), dtype=int)
    c1, N = 0, 0
    for key in MK_dict.keys():
        state = MK_dict[key]
        M_state = state
        new_vec[key] = c1
        c1 += M_state
        N += M_state
    
    new_mat = np.identity(N, dtype="complex_")
    for i in range(len(MK_dict)):
        c1 = new_vec[i]
        for j in range(len(MK_dict)):
            c2 = new_vec[j]
            new_mat[c1][c2] = m[i][j]

    return new_mat


def create_bs(U, n, m, N):
    if n == m:
        return False
    y = np.identity(N, dtype="complex_")
    y[n][n] = U[0][0]
    y[m][m] = U[1][1]
    y[n][m] = U[0][1]
    y[m][n] = U[1][0]
    return y


def generate_cv_and_dv(alpha, K, M, N, single_mode=False):
    t = 0.99999999
    sq = np.arcsinh(1)

    if single_mode:
        # Covariance Matrix for single Mode
        cm = np.identity(M)
        cm[0][0] = np.cosh(sq)
        sm = np.zeros((M, M))
        sm[0][0] = np.sinh(sq)

        # Squeezing Matrix and Inverse Squeezing Matrix
        s = np.block([[cm, sm], [sm, cm]])
        s2 = np.block([[cm, -sm], [-sm, cm]])

        # Initial Squeezing Transform
        cv = 0.5 * s @ np.conjugate(np.transpose(s))
        phi = np.arccos(t)
        totbs = np.identity(M)
        U = [[np.cos(phi), np.sin(phi)], [np.sin(phi), -np.cos(phi)]]

        cv1 = cv
        for j in range(1, M):
            bs = create_bs(U, 0, j, M)
            bt = block_diag(bs, np.conjugate(bs))
            totbs = bs @ totbs
            cv = bt @ cv @ np.conjugate(np.transpose(bt))

        # final squeezing transform
        cv1 = s2 @ cv @ np.conjugate(np.transpose(s2))

        # displacement vector
        dv = np.zeros(M)
        dv[0] = dv[0] + alpha[0]
        for j in range(1, len(alpha)):
            bs = create_bs(U, 0, j, M)
            dv = bs @ dv
            dv[0] = dv[0] + alpha[j]
        dv1 = dv
        # big dv vector
        bdv1 = np.block([dv1, np.conjugate(dv1)])

        return cv1, bdv1

    # Full CV matrix for N=KxM modes
    cm = np.identity(N, dtype="complex_")
    sm = np.zeros((N, N), dtype="complex_")
    c = 0
    for j in range(K):
        cm[c][c] = np.cosh(sq)
        sm[c][c] = np.sinh(sq)
        c = c + M

    s = np.block([[cm, sm], [sm, cm]])
    s2 = np.block([[cm, -sm], [-sm, cm]])
    cv = 0.5 * s @ np.conjugate(np.transpose(s))
    phi = np.arccos(t)
    U = np.array([[np.cos(phi), np.sin(phi)], [np.sin(phi), -np.cos(phi)]], dtype="complex_")

    # displacement vector
    dv = np.zeros(M, dtype="complex_")
    dv[0] = dv[0] + alpha[0]
    bs1 = np.identity(M, dtype="complex_")
    for j in range(1, M):
        bs = create_bs(U, 0, j, M)
        bs1 = bs @ bs1
        dv = bs @ dv
        dv[0] = dv[0] + alpha[j]

    totbs = bs1
    bigdv = dv
    for j in range(K - 1):
        bigdv = np.block([bigdv, dv])
    totdv = np.block([bigdv, np.conjugate(bigdv)])

    for j in range(K - 1):
        totbs = block_diag(totbs, bs1)

    totbs = block_diag(totbs, np.conjugate(totbs))
    cv = totbs @ cv @ np.conjugate(np.transpose(totbs))
    cv = s2 @ cv @ np.conjugate(np.transpose(s2))

    return cv, totdv


def generate_u_cv_and_dv_udag(cv, dv, MK_dict):
    K = len(MK_dict)
    uint = np.fft.fft(np.eye(K))/np.sqrt(K)
    bigu = reshape_mat(uint, MK_dict)
    bigt = block_diag(bigu, np.conjugate(bigu))
    cv = bigt @ cv @ np.conjugate(np.transpose(bigt))
    dv = bigt @ dv

    return cv, dv


def reduced_cv_and_dv(cv, dv):
    # reduced CV (delete row & colummn) and DV (delete row) for herald modes
    rcv = delete_cv(cv, [1])
    rdv = delete_vec(dv, [1])
    rbv = np.block([rdv, np.conjugate((rdv))])

    return rcv, rdv
