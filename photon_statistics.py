import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from mpl_toolkits.mplot3d import Axes3D
from scipy.linalg import block_diag
from scipy.special import genlaguerre


#Implementation of the paper photon statistics for two mode squeezed states
def probability(n1, n2, a1, a2, r):
    p = min(n1,n2)
    q = max(n1,n2)
    tanhr = np.tanh(r)
    coshr = np.cosh(r)
    u1 = a1 + np.conj(a2)*tanhr
    u2 = a2 + np.conj(a1)*tanhr
    u = (u1*u2)/tanhr
    L = genlaguerre(p,q-p)(u)

    # Amplitude C(n1,n2)
    t = -1*tanhr
    pf = math.factorial(p)
    qf = math.factorial(q)
    earg = -1*(np.conj(a1)*u1 + np.conj(a2)*u2)/2
    exp = math.exp(earg)
    c = (np.power(t,p)/coshr) * np.power((pf/qf),0.5) * np.power(u1,n1-p) * np.power(u2,n2-p)* L * exp

    #Probability Distribution P(n,n)
    if n1 == n2:
        n = n1
        absL = np.power(abs(L),2)
        earg2 = -1*(np.power(abs(a1),2) + np.power(abs(a2),2) + tanhr*(a1*a2 + np.conj(a1)*np.conj(a2)))
        prob = (np.power(tanhr,2*n)/(coshr*coshr)) * absL * math.exp(earg2)
    #Probability Distribution P(n1,n2)
    else:
        prob = np.power(abs(c),2)
    return prob

#two mode squeezed state sample configuration.
diagonal = False
a1 = 3
a2 = 3
r = 1.5

if diagonal:
    # Two Dimensional n1 = n2 plot
    x = [i for i in range(50)]
    prob = []
    for n in range(50):
        n1 = n
        n2 = n
        prob.append(probability(n1, n2, a1, a2, r))
    fig = plt.figure()
    ax = plt.axes()
    ax.plot(x, prob)
    plt.show()

else:
    # Three Dimensional plot for different values of n1, n2
    n1 = np.arange(0, 2, 1)
    n2 = np.arange(0, 2, 1)
    x, y = np.meshgrid(n1, n2)
    prob = np.array([probability(n1,n2,a1,a2,r) for n1,n2 in zip(np.ravel(x), np.ravel(y))])
    z = prob.reshape(x.shape)
    print(z)
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.set_zlim(0, 1)
    ax.zaxis.set_major_locator(LinearLocator(10))
    # A StrMethodFormatter is used automatically
    ax.zaxis.set_major_formatter('{x:.02f}')
    plt.show()



    






