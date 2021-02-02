#!/usr/bin/env python3

from hydrogen  import *
from qmc_stats import *

norm_gauss = 1./(2.*np.pi)**(1.5)
def gaussian(r):
    return norm_gauss * np.exp(-np.dot(r,r)*0.5)

def MonteCarlo(a,nmax):
    E = 0.
    N = 0.
    for istep in range(nmax):
        r = np.random.normal(loc=0., scale=1.0, size=(3))
        w = psi(a,r)
        w = w*w / gaussian(r)
        N += w
        E += w * e_loc(a,r)
    return E/N

a = 1.2
nmax = 100000
X = [MonteCarlo(a,nmax) for i in range(30)]
E, deltaE = ave_error(X)
print(f"E = {E} +/- {deltaE}")
