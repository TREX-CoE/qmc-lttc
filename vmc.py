from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a,tau,nmax):
    sq_tau = np.sqrt(tau)

    # Initialization
    E = 0.
    N = 0.
    r_old = np.random.normal(loc=0., scale=1.0, size=(3))

    for istep in range(nmax):
        d_old = drift(a,r_old)
        chi = np.random.normal(loc=0., scale=1.0, size=(3))
        r_new = r_old + tau * d_old + chi*sq_tau
        N += 1.
        E += e_loc(a,r_new)
        r_old = r_new
    return E/N


a = 0.9
nmax = 100000
tau = 0.2
X = [MonteCarlo(a,tau,nmax) for i in range(30)]
E, deltaE = ave_error(X)
print(f"E = {E} +/- {deltaE}")
