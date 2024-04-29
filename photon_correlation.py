import math
import pprint

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import unitary_group

from driver_code import find_probabilities


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


# we will compute the correlation matrix in this segment
def calculate_correlation_matrix(states, prob_to_display):
    N = len(states)
    correlation_matrix = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            correlation_matrix[i][j] = ij_correlation_value(prob_to_display,i,j)

    return correlation_matrix


# cutoff = 5
# surface_map = True
# states = (np.array([1,1]), np.array([1, 1]), np.array([1, 1]), np.array([1, 1]), np.array([1, 1]), np.array([1, 1]))
# K = len(states)
# if(K!=2):
#     surface_map = False
# uint = np.fft.fft(np.eye(K))/np.sqrt(K)
# prob_to_display = find_probabilities(states, cutoff, uint, surface_map)


# Generating 100 random unitaries and plotting a histogram wrt to the C1,2 values.
cutoff = 3
states = (np.array([1,1]),np.array([1,1]))
K = len(states)
correlation_values = []
for _ in range(1000):
    uint = unitary_group.rvs(K)
    probability = find_probabilities(states, cutoff, uint, surface_map = False)
    correlation_01 = ij_correlation_value(probability, 0, 1)
    correlation_values.append(correlation_01)

plt.hist(correlation_values, density=True, bins=150, edgecolor='black', linewidth=1)
plt.title('M = 0')
plt.xlabel('$C_{12}$')
plt.ylabel('number of incidences')
plt.show()
