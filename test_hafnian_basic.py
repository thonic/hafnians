import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from scipy.linalg import block_diag
from thewalrus import hafnian_batched
from thewalrus.quantum import Amat, Qmat


r  = 1.5
a1 = 3
a2 = 3
a1c = np.conj(a1)
a2c = np.conj(a2)
coshr  = np.cosh(r)
sinhr  = np.sinh(r)
tanhr  = np.tanh(r)
sinh2r = np.sinh(2*r)
cosh2r = np.cosh(2*r)

cutoff = 2
covariance_matrix = 0.5*np.array([[cosh2r, 0, 0, -sinh2r], [0, cosh2r, -sinh2r, 0], [0, -sinh2r, cosh2r, 0], [-sinh2r, 0, 0, cosh2r]])
displacement_vector = np.array([a1, a2, a1c, a2c])
Q = Qmat(covariance_matrix, hbar=1)
detq = np.linalg.det(Q)
A = Amat(covariance_matrix, hbar=1)
Q_inv = np.linalg.inv(Q)
exponential_factor = -0.5 * np.conj(np.transpose(displacement_vector))@Q_inv@displacement_vector
exponential_factor = math.exp(exponential_factor)

hafnians = hafnian_batched(A, cutoff, mu=displacement_vector, renorm=True)
hafnians = hafnians/np.sqrt(detq)
hafnians = hafnians*exponential_factor

print(hafnians[0][0][0][0])
print(hafnians[0][1][0][1])
print(hafnians[1][0][1][0])
print(hafnians[1][1][1][1])