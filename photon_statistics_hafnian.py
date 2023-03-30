import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator

from thewalrus import loop_hafnian
from thewalrus.quantum import Amat, Qmat, Xmat


def probability(covariance_matrix, displacement_vector, pattern):
    N = len(covariance_matrix)//2
    Q = covariance_matrix + np.identity(2*N)/2
    A = Xmat(N)@(np.identity(2*N)-np.linalg.inv(Q))
    replacement_vector = np.linalg.inv(Q)@displacement_vector
    
    detq = np.linalg.det(Q)
    exp_factor = -0.5*np.transpose(np.conj(displacement_vector))@np.linalg.inv(Q)@displacement_vector
    exp_norm = math.exp(exp_factor)
    factorial_norm = 1
    for i in pattern[:N]:
        factorial_norm = factorial_norm*math.factorial(i)
    
    hafnians = loop_hafnian(A, replacement_vector, pattern, True)
    hafnians = hafnians/(math.sqrt(detq)*factorial_norm)
    hafnians = hafnians*exp_norm

    return hafnians
    

#  Reproducing results of photon statistics using the hafnian
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

# Two mode squeezed states
cutoff = 15
covariance_matrix = 0.5*np.array([[cosh2r, 0, 0, -sinh2r], [0, cosh2r, -sinh2r, 0], [0, -sinh2r, cosh2r, 0], [-sinh2r, 0, 0, cosh2r]])
displacement_vector = np.array([a1, a2, a1c, a2c])

n1 = np.arange(0, cutoff, 1)
n2 = np.arange(0, cutoff, 1)
x, y = np.meshgrid(n1, n2)
prob = np.array([probability(covariance_matrix, displacement_vector, [n1,n2,n1,n2]).real for n1,n2 in zip(np.ravel(x), np.ravel(y))])
z = prob.reshape(x.shape)

# Plotting a Surface Map for the Probabilities
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_zlim(0, 0.1)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter('{x:.02f}')
plt.show()
