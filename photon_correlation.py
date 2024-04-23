import math
import pprint

import numpy as np

from driver_code import find_probabilities

cutoff = 4
surface_map = True
states = (np.array([1,1]), np.array([1, 1]), np.array([1, 1]))
if(len(states)!=2):
    surface_map = False

prob_to_display = find_probabilities(states, cutoff, surface_map)
# pprint.pprint(prob_to_display)


# p(ni,nj) from p(n1,n2,n3,....nk)
def marginal_probability2(prob, i, n_i, j, n_j):
    probability = 0
    for key in prob.keys():
        if(key[i] == n_i and key[j] == n_j):
            probability += prob[key]

    return probability


# p(ni) from p(n1,n2,n3,....nk)
def marginal_probability1(prob, i, n_i):
    probability = 0
    for key in prob.keys():
        if(key[i] == n_i):
            probability += prob[key]

    return probability


# This function finds the correlation between mode i and j  ===  <n_i n_j>
def correlation_function(prob, i, j):
    correlation = {}
    total_correlation = 0
    for n_i in range(cutoff):
        for n_j in range(cutoff):
            correlation_value = 0
            for key in prob.keys():
                if(key[i] == n_i and key[j] == n_j and (n_i,n_j) not in correlation.keys()):
                    marginal_probability = marginal_probability2(prob, i, n_i, j, n_j)
                    correlation_value += n_i*n_j*marginal_probability
                    correlation[(n_i,n_j)] = correlation_value
    for key in correlation.keys():
        total_correlation += correlation[key]

    return total_correlation


# This function finds the mean photons in mode i  ===  <n_i>
def mean_function(prob, i):
    mean = {}
    total_mean = 0
    for n_i in range(cutoff):
        mean_value = 0
        for key in prob.keys():
            if(key[i] == n_i and n_i not in mean.keys()):
                marginal_probability = marginal_probability1(prob, i, n_i)
                mean_value += n_i*marginal_probability
                mean[n_i] = mean_value
    for key in mean.keys():
        total_mean += mean[key]
    
    return total_mean


# We need the value of <ni * nj> - <ni><nj> 
def ij_correlation_value(prob, i, j):
    return correlation_function(prob,i,j) - mean_function(prob,i)*mean_function(prob,j)


# i and j are mode numbers and they start from 0 -> len(states)-1
i = 1
j = 2
print("Correlation value between mode number {} and {} is {}".format(i,j,ij_correlation_value(prob_to_display,i,j)))