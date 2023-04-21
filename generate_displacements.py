# Implementation of CERF paper 2005 to generate displacement parameters alpha 

import numpy as np
import math

# Input coefficients of the desired state in fock basis
cp = np.array([1,1,1,1])
cp = cp/np.sqrt((np.conj(np.transpose(cp))@cp))

# Beam Splitter Parameters
t = 0.99999999
phi = np.arccos(t)
sq = np.arcsinh(1)
n = len(cp)

#  Max photon number
nm=n-1
N=2*n
maxn=N
am=np.zeros((maxn+2,maxn+2))
ad=np.zeros((maxn+2,maxn+2))
cvac=np.zeros(n)
cvac[0] = 1

# Filling a and a+ operators
for j in range(maxn+1):
    am[j][j+1] = np.sqrt(j+1)
    ad[j+1][j] = np.sqrt(j+1)

A = np.cosh(sq)*am + ad
Am = {}
for j in range(n+1):
    Aj = np.linalg.matrix_power(A,j)
    Am[j] = Aj[:n,:n]

hvec = np.zeros(n)
h = hvec
cp2 = cp
for j in range(n):
    h[j] = cp2[n-j-1]/np.sqrt(math.factorial(nm-j))
    cp2 = cp2 - h[j]*Am[n-j-1]@cvac

beta = np.roots(h)
m = np.zeros((n-1,n-1))
for j in range(n-1,0,-1):
    for k in range(n-1,j-1,-1):
        m[j-1][k-1] = np.power(t,(n-1)-k)

m_inv = np.linalg.inv(m)
alpha = np.zeros(n)
alpha[1:] = m_inv@beta

s1=0
s2=0
for j in range(1,n):
    s1=s1+alpha[j]*np.power(t,n-j)
    s2=s2+np.conj(alpha[j])*np.power(t,j-n)

s1=s1*np.cosh(sq)
sdiff = (s2-s1)/np.cosh(sq)
x = np.real(sdiff)/(np.power(t,n)-np.power(t,-n)/np.cosh(sq))
y = np.imag(sdiff)/(np.power(t,n)+np.power(t,-n)/np.cosh(sq) )
alpha[0]=x+1j*y

print(alpha)




