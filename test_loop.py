import argparse
from datetime import datetime
from thewalrus.reference import hafnian as haf_ref
from thewalrus import hafnian
import numpy as np
from thewalrus.quantum import Amat, Qmat, reduced_gaussian, density_matrix_element, density_matrix
from thewalrus.samples import hafnian_sample_state, generate_hafnian_sample
from thewalrus import hafnian, loop_hafnian
from time import time
import logging
from pathlib import Path

from blackbird_code.subtraction import state_generation


# covariance_matrix = np.array([[np.sqrt(4), 0], [0, np.sqrt(2)]])
#covariance_matrix = np.array([[0, -0.5, 0, 0], [-0.5, 0, 0, 0], [0, 0, 0, -0.5], [0, 0, -0.5, 0]])

covariance_matrix = np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])

D = np.array([3, 3, 3, 3])


#covariance_matrix = np.array([[0, 0, 0, 1,0, 0], [0, 0, 0, 1, 1, 0], [0, 0, 0, 1, 1, 1], [1, 1, 1, 0, 0, 0],[0, 1, 1, 0, 1, 0],[0, 0, 1, 0, 0, 0]])
#D = np.array([1, 0, 0, 0, 1, 0])


At = np.array([[1,0,0,1,0,0],
             [0,0,0,1,1,0],
             [0,0,0,1,1,1],
             [1,1,1,0,0,0],
             [0,1,1,0,1,0],
             [0,0,1,0,0,0]])

At2 = np.array([[1.5, -0.5, 0, 0], [-0.5, 1.5, 0, 0], [0, 0, 1.5, -0.5], [0, 0, -0.5, 1.5]])

#haf = loop_hafnian(covariance_matrix, D=D)

haf= haf_ref(At2, loop=True)

print("done")