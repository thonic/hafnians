# inputs: covariance, displacement, which photons to measure

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

logging.basicConfig(level=logging.DEBUG)

parser = argparse.ArgumentParser(description="Calculate samples from Hafnian of a Gaussian State")

parser.add_argument(
    "--covariance_file",
    dest="covariance_file",
    type=str,
    help="Path to the CSV file with the input Covariance Matrix.",
    default="covariance_matrix.csv",
)

parser.add_argument(
    "--mean_file", dest="mean_file", type=str, help="Path to the CSV file with the mean vector.", default="mean.csv"
)

parser.add_argument(
    "--post_select_file",
    dest="post_select_file",
    type=str,
    help="Path to the CSV file with the indices to consider",
    default="post_select.csv",
)
parser.add_argument(
    "--job_name",
    dest="job_name",
    type=str,
    help="Part of the filename with the results.",
    default="hafnian_calc",
)
args = parser.parse_args()


job_name = args.job_name
dateTimeObj = datetime.now()
output_fname = f"{dateTimeObj.strftime('%Y-%m-%d %H.%M.%S')} {job_name}.txt"
output_fname = output_fname.replace(" ", "_")
logging.info(f"Result will be saved to {output_fname}.")

output_path = Path("results")
output_path.mkdir(parents=True, exist_ok=True)
output_path = output_path / output_fname

covariance_matrix = np.genfromtxt(args.covariance_file, delimiter=",")
mean = np.genfromtxt(args.mean_file, delimiter=",")
post_select = np.atleast_1d(np.genfromtxt(args.post_select_file, delimiter=",", dtype=int))


# check the inputs are valid
assert (
    covariance_matrix.shape[0] == covariance_matrix.shape[1]
), f"Covariance matrix is not a square matrix (got {covariance_matrix.shape[0]}x{covariance_matrix.shape[1]})."
assert (
    covariance_matrix.shape[0] == mean.shape[0]
), f"Covariance matrix and the mean vector dimensions do not match (got {covariance_matrix.shape[0]}x{covariance_matrix.shape[1]} and {mean.shape[0]})."
assert (
    post_select.shape[0] == mean.shape[0] / 2
), f"Mean vector and post select dimensions do not match (got {mean.shape[0]} and {post_select.shape[0]})."


post_select = {k: v for k, v in enumerate(post_select) if v != -1}

# dm = density_matrix(mean, covariance_matrix, post_select, hbar=1, normalize=True)
A = Amat(covariance_matrix, hbar=1)
haf = loop_hafnian(A, D=mean, reps=None, glynn=True)

headers = ["hafnian"]
with open(output_path, 'a') as csvfile:
    for h, r in zip(headers, (haf, )):
        np.savetxt(csvfile, np.atleast_1d(r), "%.5e", delimiter=', ', newline='\n', header=h, footer='')
