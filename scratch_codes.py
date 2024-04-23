import numpy as np

def rearrange_cv_and_dv(cvs, dvs, count, size):
    covariance_matrix = np.zeros((size * 2, size * 2), dtype="complex_")
    displacement_vector = np.zeros(size * 2, dtype="complex_")
    counter, i, j, m, ocount = 0, 0, 0, 0, 0

    while i <= 2 * size - 1:
        xinc = yinc = ssize = int(len(cvs[counter % count]) / 2)
        # dim = size / ssize
        dim = count
        ocount = ocount % dim
        icount, n, j = 0, 0, 0
        displacement_vector[i : i + xinc] = dvs[counter % count][m : m + xinc]
        while j <= 2 * size - 1:
            icount = icount % dim
            if ocount == icount:
                covariance_matrix[i : i + xinc, j : j + yinc] = cvs[counter % count][m : m + xinc, n : n + yinc]
                n += yinc
            j += yinc
            icount += 1
        ocount += 1
        i += xinc
        counter += 1
        if counter % count == 0:
            m += xinc
    return covariance_matrix, displacement_vector


def find_probabilities(states, cutoff, is_gaussian, same_type_nongaussian, surface_map=True, NG_with_vaccum = True):
    if not is_gaussian: 
        #Driver code for different type of input states
        beta, kappa, M = 1.5, 0.5, 0
        # states = (np.array([1, 1, 1]), np.array([1, 1, 1]), np.array([1, 1, 1]), np.array([1, 1, 1]))
        K = len(states)
        cvs, dvs, count, size = {}, {}, 0, 0
        for cp in states:
            size += len(cp)
            cp = cp / np.sqrt((np.conj(np.transpose(cp)) @ cp))
            M = len(cp)
            k = 1
            alpha = generate_alpha(cp)
            cv, dv = generate_cv_and_dv(alpha, k, M, M * k, single_mode=True)
            cvs[count] = cv
            dvs[count] = dv
            count += 1
        if count > 2:
            surface_map = False
        covariance_matrix, displacement_vector = new_rearrange_cv_and_dv(cvs, dvs, count, size)
        N = M * K
        covariance_matrix, displacement_vector = generate_u_cv_and_dv_udag(covariance_matrix, displacement_vector, K, M, N)
        deletion_array = np.zeros(2 * N)
        for i in range(2 * N):
            if i % M == 0:
                deletion_array[i] = 1
        print(deletion_array)
        reduced_cv = delete_cv(covariance_matrix, deletion_array)
        reduced_dv = delete_vec(displacement_vector, deletion_array)
        prob_herald_hafnian = probability(reduced_cv, reduced_dv, cutoff=2)
        n = (1,) * (K * (M - 1))
        prob_herald = slice_probabilities(prob_herald_hafnian, n)
        prob_hafnian_nbar = non_gaussian_probability(covariance_matrix, displacement_vector, cutoff)
    else:
        # Driver code for gaussian states
        K, M = 2, 1
        count = K
        a1, a2, r = 3, -3, 1.5
        a1c = np.conj(a1)
        a2c = np.conj(a2)
        coshr = np.cosh(r)
        sinhr = np.sinh(r)
        tanhr = np.tanh(r)
        sinh2r = np.sinh(2 * r)
        cosh2r = np.cosh(2 * r)
        covariance_matrix = 0.5 * np.array(
            [[cosh2r, 0, 0, -sinh2r], [0, cosh2r, -sinh2r, 0], [0, -sinh2r, cosh2r, 0], [-sinh2r, 0, 0, cosh2r]]
        )
        displacement_vector = np.array([a1, a2, a1c, a2c])
        prob_herald = 1
        prob_hafnian_nbar = probability(covariance_matrix, displacement_vector, cutoff)


def reshape_mat(m, M, K, N):
    if M * K != N:
        return False
    # if len(m) != K:
    #     return False

    new_mat = np.identity(N, dtype="complex_")
    c1, j1 = 0, 0
    for j1 in range(len(m)):
        c2 = 0
        for j2 in range(len(m)):
            new_mat[c1][c2] = m[j1][j2]
            c2 = c2 + M
        c1 = c1 + M
    return new_mat


def generate_u_cv_and_dv_udag(cv, dv, K, M, N):
    # cv,dv = generate_cv_and_dv(alpha, K, M, N)
    # random interferometer, this can be chosen to be anything - 50/50 beamsplitter
    # uint = [[1, 1], [-1, 1]] / np.sqrt(2)
    uint = np.fft.fft(np.eye(K))/np.sqrt(K)

    # reshape changes the 2x2 matrix to a larger matrix, KMxKM, between the modes specified in the argument
    bigu = reshape_mat(uint, M, K, N)
    bigt = block_diag(bigu, np.conjugate(bigu))
    cv = bigt @ cv @ np.conjugate(np.transpose(bigt))
    dv = bigt @ dv

    return cv, dv


def post_select_on_herald_modes(hafnian, n, K, M):
    # n is a tuple (n1,n2,...ns) indicating a particular pattern of photons in the system modes
    idx = []
    #for i in range(2 * K):
    for i in range(K):
        idx += [n[i % K]] + [1] * (M - 1)
    
    return hafnian[tuple(idx)] * np.conj(hafnian[tuple(idx)])
