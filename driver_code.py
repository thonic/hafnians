import cmath
import itertools
import math
import pprint

import numpy as np
from scipy.linalg import block_diag

from generate_cv_and_dv import *
from generate_displacements import *
from hafnian_batched_statistics import *


def diag_block_mat(tList):
    blkXsize = [ tbl.shape[1] for tbl in tList ]
    outBlocks = []
    for i, tbl in enumerate(tList):
        tBefore = np.zeros((tbl.shape[0], sum(blkXsize[:i])))
        tAfter = np.zeros((tbl.shape[0], sum(blkXsize[i+1:])))
        outBlocks.append(np.hstack([tBefore, tbl, tAfter]))
    return np.vstack(outBlocks)


def rearrange_cv_and_dv(cvs, dvs, count, size):
    cv1, cv2, cv3, cv4, dv1, dv2 = [], [], [], [], [], []
    for key in cvs.keys():
        cv = cvs[key]
        dv = dvs[key]
        dim = len(cv)
        n = int(dim/2)
        temp1 = cv[0:n, 0:n]
        temp2 = cv[0:n, n:dim]
        temp3 = cv[n:dim, 0:n]
        temp4 = cv[n:dim, n:dim]
        demp1 = dv[0:n]
        demp2 = dv[n:dim]
        cv1.append(temp1)
        cv2.append(temp2)
        cv3.append(temp3)
        cv4.append(temp4)
        dv1.append(demp1)
        dv2.append(demp2)
    block_cv1 = diag_block_mat(cv1)
    block_cv2 = diag_block_mat(cv2)
    block_cv3 = diag_block_mat(cv3)
    block_cv4 = diag_block_mat(cv4)
    covariance_matrix = np.block([[block_cv1,block_cv2],[block_cv3,block_cv4]])
    displacement_vector = np.concatenate(dv1+dv2)
    
    return covariance_matrix, displacement_vector

def find_probabilities(states, cutoff, uint, surface_map):
    MK_dict = {}
    N, cvs, dvs, count = 0, {}, {}, 0
    for cp in states:
        k = 1
        m = len(cp)
        MK_dict[count] = m
        N += m
        cp = cp / np.sqrt((np.conj(np.transpose(cp)) @ cp))
        alpha = generate_alpha(cp)
        cv, dv = generate_cv_and_dv(alpha, k, m, m * k, single_mode=True)
        cvs[count] = cv
        dvs[count] = dv
        count += 1
    
    covariance_matrix, displacement_vector = rearrange_cv_and_dv(cvs, dvs, count, size=N)
    covariance_matrix, displacement_vector = generate_u_cv_and_dv_udag(covariance_matrix, displacement_vector, MK_dict, uint)
        
    d, del_array, n = 0, np.zeros(N), []
    for key in MK_dict:
        m = MK_dict[key]
        del_array[d] = 1
        for i in range(m-1):
            n.append(1)
        d += m
    deletion_array = np.concatenate((del_array,del_array))
    reduced_cv = delete_cv(covariance_matrix, deletion_array)
    reduced_dv = delete_vec(displacement_vector, deletion_array)
    
    prob_herald_hafnian = probability(reduced_cv, reduced_dv, cutoff=2)
    prob_herald = slice_probabilities(prob_herald_hafnian, tuple(n))
    prob_hafnian_nbar = non_gaussian_probability(covariance_matrix, displacement_vector, cutoff)


    if surface_map:
        n1 = np.arange(0, cutoff, 1)
        n2 = np.arange(0, cutoff, 1)
        x, y = np.meshgrid(n1, n2)
        prob = np.array(
            [
                post_select_on_herald_modes(prob_hafnian_nbar, (n1, n2), MK_dict).real/prob_herald.real
                for n1, n2 in zip(np.ravel(x), np.ravel(y))
            ]
        )  
        x = x.flatten()
        y = y.flatten()
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        ax.bar3d(x, y, np.zeros(len(prob)), 1, 1, prob, color="lightgreen", ec="black")
        ax.set_zlim(0, 1)
        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter("{x:.02f}")
        plt.xlabel("n1")
        plt.ylabel("n2")
        plt.show()


    prob_to_display = {}
    total_probability, non_zero_probabilities = 0, 0
    ranges = []
    for i in range(count):
        ranges.append(range(0, cutoff))
    for xs in itertools.product(*ranges):
        index = tuple(xs)
        instance_probability = post_select_on_herald_modes(prob_hafnian_nbar, index, MK_dict).real/prob_herald.real
        if(instance_probability >= 0.0001):
            non_zero_probabilities += 1
            total_probability += instance_probability
            prob_to_display[index] = instance_probability
    # pprint.pprint(prob_to_display)
    # print("Number of non-zero probability states = ", non_zero_probabilities)
    # print("Total Probability = ", total_probability)
    
    return prob_to_display



# # This code is designed to just take inputs for generating covariance matrices and displacement vectors and feed into the hafnian batched.
# cutoff = 3
# surface_map = True
# states = (np.array([1,1]),np.array([1,1]),np.array([1]),np.array([1]),np.array([1]),np.array([1]),np.array([1]),np.array([1]),np.array([1]),np.array([1]),np.array([1]))
# K = len(states)
# if (K!= 2):
#     surface_map = False 
# uint = np.fft.fft(np.eye(K))/np.sqrt(K)
# prob_to_display = find_probabilities(states, cutoff, uint, surface_map)