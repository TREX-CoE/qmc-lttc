from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a, nmax):
     E = 0.
     N = 0.
     for istep in range(nmax):
          r = np.random.uniform(-5., 5., (3))
          w = psi(a,r)
          w = w*w
          N += w
          E += w * e_loc(a,r)
   return E/N

a = 0.9
nmax = 100000
X = [MonteCarlo(a,nmax) for i in range(30)]
E, deltaE = ave_error(X)
print(f"E = {E} +/- {deltaE}")
