# inputs: covariance, displacement, which photons to measure

import argparse
from datetime import datetime

import numpy as np
from thewalrus.quantum import Amat, Qmat, reduced_gaussian, density_matrix_element, density_matrix
from thewalrus.samples import hafnian_sample_state, generate_hafnian_sample
from thewalrus import hafnian, loop_hafnian, hafnian_batched
from time import time
import logging
from pathlib import Path

logging.basicConfig(level=logging.DEBUG)

# -----------------------------
# Configuration
# -----------------------------

# simple covariance matrix, two modes
covariance_matrix = np.array([[1.5000, -1.4142], [-1.4142, 1.5000]])

# displacement vector for testing
displacement_vector = np.array([1.0, 1.0])

# cutoff in number of photons
cutoff = np.array([1, 2])

# -----------------------------
# Calculation
# -----------------------------

A = Amat(covariance_matrix, hbar=1)

# see Eq. (6) in papers/boson_sampling.pdf
hafnians = hafnian_batched(A, cutoff, mu=displacement_vector, renorm=True)

print(hafnians)
