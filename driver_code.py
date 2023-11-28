import cmath
import math
import pprint

import numpy as np
from scipy.linalg import block_diag

from generate_cv_and_dv import *
from generate_displacements import *
from hafnian_batched_statistics import *

# This code is designed to just take inputs for generating covariance matrices and displacement vectors and feed into the hafnian batched.
is_gaussian = True
surface_map = False
cutoff = 3
test = True

if not is_gaussian:
    # Driver Code for non gaussian states generation of same type K * M -> generate_cv_And_dv.py and generate_displacements.py
    beta = 1.5
    kappa = 0.5    
    M = 0   
    states = (compute_coherent_coefficients(beta, kappa, n=2),compute_coherent_coefficients(-beta, kappa, n=2))
    K = len(states)
    covariance_matrix = np.zeros((8,8),dtype="complex_")
    displacement_vector = np.zeros(8,dtype="complex_")

    # cp = np.array([1,1,1], dtype = "complex_")
    for cp in states:
        cp = cp/np.sqrt((np.conj(np.transpose(cp))@cp))
        print(cp)
        M = len(cp)
        k = 1
        alpha = generate_alpha(cp)
        cv, dv = generate_cv_and_dv(alpha,k,M,M*k,single_mode=True)
        # if covariance_matrix is None and displacement_vector is None:
        covariance_matrix[0:2,0:2] = covariance_matrix[2:4,2:4] = cv[0:2,0:2]
        covariance_matrix[0:2,4:6] = covariance_matrix[2:4,6:8] = cv[0:2,2:4]
        covariance_matrix[4:6,0:2] = covariance_matrix[6:8,2:4] = cv[2:4,0:2]
        covariance_matrix[4:6,4:6] = covariance_matrix[6:8,6:8] = cv[2:4,2:4]

        displacement_vector[0:2] = dv[0:2]
        displacement_vector[2:4] = dv[0:2]
        displacement_vector[4:6] = dv[2:4]
        displacement_vector[6:8] = dv[2:4]
        
        # else:
        #     covariance_matrix = block_diag(covariance_matrix,cv)
        #     displacement_vector = np.block([displacement_vector,dv])

        print("\nCOVARIANCE MATRIX X= \n",covariance_matrix)
        print("\nDISPLACEMENT VECTOR X= \n",displacement_vector)
    
    # we will be calculating prob herald for getting one photon in all the herald modes by reducing the cv and dv.
    N = M*K
    covariance_matrix, displacement_vector = generate_u_cv_and_dv_udag(covariance_matrix,displacement_vector,K,M,N)
    # print("\nCOVARIANCE MATRIX = \n",covariance_matrix)
    # print("\nDISPLACEMENT VECTOR = \n",displacement_vector)
    deletion_array = np.zeros(2*N)
    for i in range(2*N):
        if i%M == 0:
            deletion_array[i] = 1

    reduced_cv = delete_cv(covariance_matrix, deletion_array)
    reduced_dv = delete_vec(displacement_vector, deletion_array)
    print("\nReduced COVARIANCE MATRIX = \n",reduced_cv)
    print("\nReduced DISPLACEMENT VECTOR = \n",reduced_dv)
    prob_herald_hafnian = probability(reduced_cv, reduced_dv, cutoff=cutoff)
    n = (1,) * (K*(M-1))
    prob_herald = slice_probabilities(prob_herald_hafnian, n)
    prob_hafnian_nbar = non_gaussian_probability(covariance_matrix, displacement_vector, cutoff, prob_herald)

elif test:
    beta = 1.5
    kappa = 0.5     
    # cp = compute_coherent_coefficients(beta, kappa, n=4)
    cp = [1,1]
    cp = cp/np.sqrt((np.conj(np.transpose(cp))@cp))    
    M = len(cp)
    K = 4
    N = M*K
    alpha = generate_alpha(cp)
    cv, dv = generate_cv_and_dv(alpha,K,M,N)
    print("shape =",cv.shape,dv.shape)
    covariance_matrix, displacement_vector = generate_u_cv_and_dv_udag(cv,dv,K,M,N)
    print("\nCOVARIANCE MATRIX = \n",covariance_matrix)
    print("\nDISPLACEMENT VECTOR = \n",displacement_vector)
    deletion_array = np.zeros(2*N)
    for i in range(2*N):
        if i%M == 0:
            deletion_array[i] = 1

    reduced_cv = delete_cv(covariance_matrix, deletion_array)
    reduced_dv = delete_vec(displacement_vector, deletion_array)
    prob_herald_hafnian = probability(reduced_cv, reduced_dv, cutoff=cutoff)
    n = (1,) * (K*(M-1))
    prob_herald = slice_probabilities(prob_herald_hafnian, n)
    prob_hafnian_nbar = non_gaussian_probability(covariance_matrix, displacement_vector, cutoff, prob_herald)



    
else:
    # Driver code for gaussian states, specifically two mode squeezed states from hafnian_batched_statistics.py
    K = 2
    M = 1
    a1 = 3
    a2 = -3
    r  = 1.5
    a1c = np.conj(a1)
    a2c = np.conj(a2)
    coshr  = np.cosh(r)
    sinhr  = np.sinh(r)
    tanhr  = np.tanh(r)
    sinh2r = np.sinh(2*r)
    cosh2r = np.cosh(2*r)
    covariance_matrix = 0.5*np.array([[cosh2r, 0, 0, -sinh2r], [0, cosh2r, -sinh2r, 0], [0, -sinh2r, cosh2r, 0], [-sinh2r, 0, 0, cosh2r]])
    displacement_vector = np.array([a1, a2, a1c, a2c])
    prob_hafnian_nbar = probability(covariance_matrix, displacement_vector, cutoff)

if surface_map:
    # Feed the covariance matrix into the hafnian_batched_function and generate plots for 2 modes, for higher modes, see below!
    n1 = np.arange(0, cutoff, 1)
    n2 = np.arange(0, cutoff, 1)
    x, y = np.meshgrid(n1, n2)
    prob = np.array([post_select_on_herald_modes(prob_hafnian_nbar,(n1,n2),K,M).real for n1,n2 in zip(np.ravel(x), np.ravel(y))])
    prob_to_display = {}
    for i in n1:
        for j in n2:
            prob_to_display[(i,j)] = post_select_on_herald_modes(prob_hafnian_nbar,(i,j),K,M).real
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
    ax.bar3d(x, y, np.zeros(len(prob)), 1, 1, prob, color = "lightgreen", ec="black")
    ax.set_zlim(0, 1)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter('{x:.02f}')
    plt.xlabel ("n1")
    plt.ylabel ("n2")
    plt.show()

else:
    # nbar represents [n1,n2,n3......,nm] no of photons in each output modes arranged in a tuple
    nbar = (0,0) 
    n1 = np.arange(0, cutoff, 1)
    n2 = np.arange(0, cutoff, 1)
    n3 = np.arange(0, cutoff, 1)
    n4 = np.arange(0, cutoff, 1)

    prob_to_display = {}
    for i in n1:
        for j in n2:
            for k in n3:
                for l in n4:
                    prob_to_display[(i,j,k,l)] = post_select_on_herald_modes(prob_hafnian_nbar,(i,j,k,l),K=4,M=2).real
    pprint.pprint(prob_to_display)

    # prob_nbar = post_select_on_herald_modes(prob_hafnian_nbar, nbar, K=4, M=2).real
    # print("probability of given nbar is", prob_nbar)