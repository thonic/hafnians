import argparse
from datetime import datetime

import numpy as np
from thewalrus.quantum import Amat, Qmat, reduced_gaussian, density_matrix_element, density_matrix
from thewalrus.samples import hafnian_sample_state, generate_hafnian_sample
from thewalrus import hafnian, loop_hafnian
from time import time
import logging
from pathlib import Path

from blackbird_code.subtraction import state_generation


covariance_matrix = np.array([[np.sqrt(4), 0], [0, np.sqrt(2)]])
D = np.array([2, 2])
haf = loop_hafnian(covariance_matrix, D=D)

print("done")