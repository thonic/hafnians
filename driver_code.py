import cmath
import itertools
import math
import pprint

import numpy as np
from scipy.linalg import block_diag

from generate_cv_and_dv import *
from generate_displacements import *
from hafnian_batched_statistics import *


def rearrange_cv_and_dv(cvs, dvs, count, size):
    covariance_matrix = np.zeros((size * 2, size * 2), dtype="complex_")
    displacement_vector = np.zeros(size * 2, dtype="complex_")
    counter, i, j, m, ocount = 0, 0, 0, 0, 0
    xinc = yinc = ssize = int(len(cvs[counter % count]) / 2)
    dim = size / ssize
    while i < 2 * size - 1:
        ocount = ocount % dim
        icount, n, j = 0, 0, 0
        displacement_vector[i : i + xinc] = dvs[counter % count][m : m + xinc]
        while j < 2 * size - 1:
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


# This code is designed to just take inputs for generating covariance matrices and displacement vectors and feed into the hafnian batched.
cutoff = 3
is_gaussian = False
surface_map = True
same_type_nongaussian = True


if not is_gaussian:
    beta, kappa, M = 1.5, 0.5, 0
    states = (np.array([1, 1, 1, 1]), np.array([1, 1, 1, 1]))
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
    covariance_matrix, displacement_vector = rearrange_cv_and_dv(cvs, dvs, count, size)
    N = M * K
    covariance_matrix, displacement_vector = generate_u_cv_and_dv_udag(covariance_matrix, displacement_vector, K, M, N)
    deletion_array = np.zeros(2 * N)
    for i in range(2 * N):
        if i % M == 0:
            deletion_array[i] = 1
    reduced_cv = delete_cv(covariance_matrix, deletion_array)
    reduced_dv = delete_vec(displacement_vector, deletion_array)
    prob_herald_hafnian = probability(reduced_cv, reduced_dv, cutoff=cutoff)
    n = (1,) * (K * (M - 1))
    prob_herald = slice_probabilities(prob_herald_hafnian, n)
    prob_hafnian_nbar = non_gaussian_probability(covariance_matrix, displacement_vector, cutoff, prob_herald)


elif same_type_nongaussian:
    # Driver code for multiple but same non gaussian states
    beta = 1.5
    kappa = 0.5
    # cp = compute_coherent_coefficients(beta, kappa, n=4)
    cp = [1, 1]
    cp = cp / np.sqrt((np.conj(np.transpose(cp)) @ cp))
    M = len(cp)
    count = K = 5
    N = M * K
    if count > 2:
        surface_map = False
    alpha = generate_alpha(cp)
    cv, dv = generate_cv_and_dv(alpha, K, M, N)
    covariance_matrix, displacement_vector = generate_u_cv_and_dv_udag(cv, dv, K, M, N)
    deletion_array = np.zeros(2 * N)
    for i in range(2 * N):
        if i % M == 0:
            deletion_array[i] = 1
    reduced_cv = delete_cv(covariance_matrix, deletion_array)
    reduced_dv = delete_vec(displacement_vector, deletion_array)
    prob_herald_hafnian = probability(reduced_cv, reduced_dv, cutoff=cutoff)
    n = (1,) * (K * (M - 1))
    prob_herald = slice_probabilities(prob_herald_hafnian, n)
    prob_hafnian_nbar = non_gaussian_probability(covariance_matrix, displacement_vector, cutoff, prob_herald)


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
    prob_hafnian_nbar = probability(covariance_matrix, displacement_vector, cutoff)


if surface_map:
    # Feed the covariance matrix into the hafnian_batched_function and generate plots for 2 modes, for higher modes, see below!
    n1 = np.arange(0, cutoff, 1)
    n2 = np.arange(0, cutoff, 1)
    x, y = np.meshgrid(n1, n2)
    prob = np.array(
        [
            post_select_on_herald_modes(prob_hafnian_nbar, (n1, n2), K, M).real
            for n1, n2 in zip(np.ravel(x), np.ravel(y))
        ]
    )
    prob_to_display = {}
    ranges = []
    for i in range(count):
        ranges.append(range(0, cutoff))
    for xs in itertools.product(*ranges):
        index = tuple(xs)
        prob_to_display[index] = post_select_on_herald_modes(prob_hafnian_nbar, index, K, M).real
    pprint.pprint(prob_to_display)

    # # to plot a surface map
    # z = prob.reshape(x.shape)
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    # ax.set_zlim(0, 1)
    # ax.zaxis.set_major_locator(LinearLocator(10))
    # ax.zaxis.set_major_formatter('{x:.02f}')
    # plt.show()

    # to plot a histogram
    x = x.flatten()
    y = y.flatten()
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.bar3d(x, y, np.zeros(len(prob)), 1, 1, prob, color="lightgreen", ec="black")
    ax.set_zlim(0, 1)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter("{x:.02f}")
    plt.xlabel("n1")
    plt.ylabel("n2")
    plt.show()


else:
    prob_to_display = {}
    ranges = []
    for i in range(count):
        ranges.append(range(0, cutoff))
    for xs in itertools.product(*ranges):
        index = tuple(xs)
        prob_to_display[index] = post_select_on_herald_modes(prob_hafnian_nbar, index, K, M).real
    pprint.pprint(prob_to_display)
