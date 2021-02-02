#!/usr/bin/env python3

from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a, nmax):
     energy = 0.
     normalization = 0.

     for istep in range(nmax):
          r = np.random.uniform(-5., 5., (3))

          w = psi(a,r)
          w = w*w

          energy        += w * e_loc(a,r)
          normalization += w

     return energy / normalization

a    = 1.2
nmax = 100000

X = [MonteCarlo(a,nmax) for i in range(30)]
E, deltaE = ave_error(X)

print(f"E = {E} +/- {deltaE}")
