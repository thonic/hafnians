import numpy as np
import pytest

from driver_code import rearrange_cv_and_dv


def test_rearrange_basic():
    # Basic test to check if the function rearranges a simple set of CVs and DVs correctly
    cvs = {
        0:np.array([[1, 1], [1, 1]]),    
        1:np.array([[-1, -1], [-1, -1]])
    }
    dvs = {
        0 : np.array([0.5, 0.5]),
        1 : np.array([-0.5, -0.5])
    }
    count = 2
    size = 2
    expected_covariance_matrix = np.array(
        [
            [1.0 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, -1.0 + 0.0j, 0.0 + 0.0j, -1.0 + 0.0j],
            [1.0 + 0.0j, 0.0 + 0.0j, 1.0 + 0.0j, 0.0 + 0.0j],
            [0.0 + 0.0j, -1.0 + 0.0j, 0.0 + 0.0j, -1.0 + 0.0j],
        ],
        dtype="complex_",
    )
    expected_displacement_vector = np.array([0.5 + 0.0j, -0.5 + 0.0j, 0.5 + 0.0j, -0.5 + 0.0j], dtype="complex_")
    covariance_matrix, displacement_vector = rearrange_cv_and_dv(cvs, dvs, count, size)
    np.testing.assert_array_equal(covariance_matrix, expected_covariance_matrix)
    np.testing.assert_array_equal(displacement_vector, expected_displacement_vector)
