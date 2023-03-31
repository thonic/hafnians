# inputs: covariance, displacement, which photons to measure
import numpy as np
from thewalrus.thewalrus import hafnian_batched
from thewalrus.thewalrus.quantum import Amat, Qmat

# simple covariance matrix, two modes
covariance_matrix = np.identity(4)
mu = np.array([2] * 4)

Q = Qmat(covariance_matrix, hbar=1)
detq = np.linalg.det(Q)

cutoff = 2
A = Amat(covariance_matrix, hbar=1)

hafnians = hafnian_batched(A, cutoff, mu=mu, renorm=True)
hafnians = hafnians / np.sqrt(detq)

# print(hafnians)
print(hafnians[0][0][0][0])
print(hafnians[0][1][0][1])
print(hafnians[1][0][1][0])
print(hafnians[1][1][1][1])
