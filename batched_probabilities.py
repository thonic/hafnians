# inputs: covariance, displacement, which photons to measure
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import block_diag
from thewalrus.thewalrus import hafnian, hafnian_batched, loop_hafnian
from thewalrus.thewalrus.quantum import (Amat, Qmat, density_matrix,
                               density_matrix_element, reduced_gaussian)
from thewalrus.thewalrus.samples import generate_hafnian_sample, hafnian_sample_state

# logging.basicConfig(level=logging.DEBUG)
r = 0.5
sinh2r = np.sinh(2*r)
cosh2r = np.cosh(2*r)

# Interferometer matrix
fourier_matrix  = np.fft.fft(np.eye(2))/np.sqrt(2)
# fourier_matrix = np.array([[1, 1j],[1j, 1]])/np.sqrt(2)

fourier_conjugate = np.conjugate(fourier_matrix)
fourier_transpose = np.transpose(fourier_matrix)
fourier_dagger = np.transpose(fourier_conjugate)
interferometer = block_diag(fourier_matrix, fourier_conjugate)
interferometer_dagger = block_diag(fourier_dagger, fourier_transpose)

# simple covariance matrix, two modes
covariance_matrix = np.array([[cosh2r, 0, 0, -sinh2r], [0, cosh2r, -sinh2r, 0], [0, -sinh2r, cosh2r, 0], [-sinh2r, 0, 0, cosh2r]])
covariance_matrix = 0.5*interferometer @ covariance_matrix @ interferometer_dagger

Q = Qmat(covariance_matrix, hbar=1)
detq = np.linalg.det(Q)

cutoff = 3
A = Amat(covariance_matrix, hbar=1)

# see Eq. (6) in papers/boson_sampling.pdf
hafnians = hafnian_batched(A, cutoff, mu=None, renorm=True)
hafnians = hafnians/np.sqrt(detq)

# print(hafnians)
print(hafnians[0][0][0][0])
print(hafnians[0][1][0][1])
print(hafnians[1][0][1][0])
print(hafnians[1][1][1][1])
print(hafnians[0][2][0][2])
print(hafnians[2][0][2][0])
print(hafnians[1][2][1][2])
print(hafnians[2][1][2][1])
print(hafnians[2][2][2][2])

if cutoff == 2:
    # for 0 or 1 photon present in the output mode
    data2D = [[hafnians[0][0][0][0], hafnians[0][1][0][1]],[hafnians[1][0][1][0], hafnians[1][1][1][1]]]
else:
# for 0, 1 or 2 photons present in the output mode
    data2D = [[hafnians[0][0][0][0], hafnians[0][1][0][1], hafnians[0][2][0][2]],
              [hafnians[1][0][1][0], hafnians[1][1][1][1], hafnians[1][2][1][2]],
              [hafnians[2][0][2][0], hafnians[2][1][2][1], hafnians[2][2][2][2]]]

data_array = np.array(data2D)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x_data, y_data = np.meshgrid( np.arange(data_array.shape[1]), np.arange(data_array.shape[0]) )
x_data = x_data.flatten()
y_data = y_data.flatten()
z_data = data_array.flatten()
ax.bar3d(x_data, y_data, np.zeros(len(z_data)), 1, 1, z_data )
plt.show()