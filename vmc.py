from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a,dt,nmax):
    sq_dt = np.sqrt(dt)

    # Initialization
    E = 0.
    N = 0.
    r_old = np.random.normal(loc=0., scale=1.0, size=(3))

    for istep in range(nmax):
        d_old = drift(a,r_old)
        chi = np.random.normal(loc=0., scale=1.0, size=(3))
        r_new = r_old + dt * d_old + chi*sq_dt
        N += 1.
        E += e_loc(a,r_new)
        r_old = r_new
    return E/N


a = 0.9
nmax = 100000
dt = 0.2
X = [MonteCarlo(a,dt,nmax) for i in range(30)]
E, deltaE = ave_error(X)
print(f"E = {E} +/- {deltaE}")
