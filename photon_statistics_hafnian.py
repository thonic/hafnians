import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import block_diag
from numpy.linalg import norm
from thewalrus import hafnian_batched
from thewalrus.quantum import Amat, Qmat


def slice_probabilities(hafnians, n1, n2):
    return hafnians[n1][n2][n1][n2]

def probability(covariance_matrix, displacement_vector, cutoff):
    Q = Qmat(covariance_matrix, hbar=1)
    A = Amat(covariance_matrix, hbar=1)
    detq = np.linalg.det(Q)
    hafnians = hafnian_batched(A, cutoff, mu=displacement_vector, renorm=True)
    hafnians = hafnians/np.sqrt(detq)
    return hafnians


#  Reproducing results of photon statistics using the hafnian
r  = 1.5
a1 = 3.0
a2 = 3.0 
a1c = np.conj(a1)
a2c = np.conj(a2)
coshr  = np.cosh(r)
sinhr  = np.sinh(r)
tanhr  = np.tanh(r)
sinh2r = np.sinh(2*r)
cosh2r = np.cosh(2*r)

# Two mode squeezed states
cutoff = 50
covariance_matrix = 0.5*np.array([[cosh2r, 0, 0, -sinh2r], [0, cosh2r, -sinh2r, 0], [0, -sinh2r, cosh2r, 0], [-sinh2r, 0, 0, cosh2r]])
displacement_vector = np.array([coshr*a1-sinhr*a2c, coshr*a2-sinhr*a1c, coshr*a1c-sinhr*a2, coshr*a2c - sinhr*a1])
# displacement_norm = norm(displacement_vector, 2)
# displacement_vector = displacement_vector/displacement_norm
print(displacement_vector)
hafnians = probability(covariance_matrix, displacement_vector, cutoff)
n1 = np.arange(0, cutoff, 1)
n2 = np.arange(0, cutoff, 1)
x, y = np.meshgrid(n1, n2)
prob = np.array([slice_probabilities(hafnians, n1, n2) for n1,n2 in zip(np.ravel(x), np.ravel(y))])
z = prob.reshape(x.shape)

# Plotting a Surface Map for the Probabilities
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_zlim(0, 1)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter('{x:.02f}')
plt.show()
