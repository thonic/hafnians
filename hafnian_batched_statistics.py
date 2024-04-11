import cmath
import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy.linalg import block_diag

from thewalrus import hafnian_batched
from thewalrus.quantum import Xmat


def slice_probabilities(hafnian, n):
    # n is a tuple (n1,n2,...ns) indicating a particular pattern of photons
    idx = []
    for i in range(2*len(n)):
        idx += [n[i % len(n)]]

    return hafnian[tuple(idx)] 

# logic to post select 1 photon from every herald modes
def post_select_on_herald_modes(hafnian, n, K, M):
    # n is a tuple (n1,n2,...ns) indicating a particular pattern of photons in the system modes
    idx = []
    #for i in range(2 * K):
    for i in range(K):
        idx += [n[i % K]] + [1] * (M - 1)
    
    return hafnian[tuple(idx)] * np.conj(hafnian[tuple(idx)])


def probability(covariance_matrix, displacement_vector, cutoff):
    N = len(covariance_matrix) // 2
    Q = covariance_matrix + np.identity(2 * N) / 2
    A = Xmat(N) @ (np.identity(2 * N) - np.linalg.inv(Q))
    replacement_vector = np.linalg.inv(Q) @ displacement_vector
    B = A[0:N,0:N]

    detq = np.linalg.det(Q)
    exp_factor = -0.5 * np.transpose(np.conj(displacement_vector)) @ np.linalg.inv(Q) @ displacement_vector
    exp_norm = cmath.exp(exp_factor)
    
    hafnians = hafnian_batched(A, cutoff, mu=replacement_vector, renorm=True)
    hafnians = hafnians / cmath.sqrt(detq)
    hafnians = hafnians * exp_norm

    return hafnians


def non_gaussian_probability(covariance_matrix, displacement_vector, cutoff):
    N = len(covariance_matrix) // 2
    Q = covariance_matrix + np.identity(2 * N) / 2
    A = Xmat(N) @ (np.identity(2 * N) - np.linalg.inv(Q))
    B = A[0:N,0:N]
    replacement_vector = np.linalg.inv(Q) @ displacement_vector
    
    detq = np.linalg.det(Q)
    exp_factor = -0.5 * np.transpose(np.conj(displacement_vector)) @ np.linalg.inv(Q) @ displacement_vector
    exp_norm = cmath.exp(exp_factor)
    
    hafnians = hafnian_batched(B, cutoff, mu=replacement_vector[0:N], renorm=True)
    hafnians = hafnians / cmath.sqrt(cmath.sqrt(detq))
    hafnians = hafnians * cmath.sqrt(exp_norm)
    
    return hafnians
