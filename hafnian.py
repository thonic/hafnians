# 20x20 haffnian (look for loop haffnian)
# inputs: covariance, displacement, which photons to measure

import argparse
import numpy as np
from thewalrus.samples import hafnian_sample_state
from time import time

# parser = argparse.ArgumentParser(description="Calculate samples from Hafnian of a Gaussian State")
#
# parser.add_argument(
#     "--covariance_file",
#     dest="covariance_file",
#     type=str,
#     help="Path to the CSV file with to input " "Covariance Matrix.",
#     required=True,
# )
# parser.add_argument(
#     "--mean_file", dest="mean_file", type=str, help="Path to the CSV file with the mean vector.", required=True
# )
# parser.add_argument(
#     "--n_samples", dest="n_samples", type=int, help="How many samples should be generated.", required=True
# )
# args = parser.parse_args()
#
#
# covariance_matrix = np.genfromtxt(args.covariance_file, delimiter=",")
# mean = np.genfromtxt(args.mean_file, delimiter=",")
# n_samples = args.n_samples

n_photons = 10
n_samples = 10
covariance_matrix = np.identity(2 * n_photons, dtype=np.int64)
np.savetxt("covariance_matrix.csv", covariance_matrix, delimiter=",")
mean = np.array([1] * 2 * n_photons, dtype=np.float64)
np.savetxt("mean.csv", mean, delimiter=",")

generated_samples = hafnian_sample_state(
    cov=covariance_matrix,
    samples=n_samples,
    mean=mean,
    hbar=1,
    cutoff=10,
    max_photons=100,
    approx=False,
    parallel=True,
)

np.savetxt("samples.csv", generated_samples, delimiter=",")
