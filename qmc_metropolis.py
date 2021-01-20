from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a,nmax,tau):
    E = 0.
    N = 0.
    N_accep = 0.
    r_old = np.random.uniform(-tau, tau, (3))
    psi_old = psi(a,r_old)
    for istep in range(nmax):
        r_new = r_old + np.random.uniform(-tau,tau,(3))
        psi_new = psi(a,r_new)
        ratio = (psi_new / psi_old)**2
        v = np.random.uniform(0,1,(1))
        if v < ratio:
            N_accep += 1.
            r_old = r_new
            psi_old = psi_new
        N += 1.
        E += e_loc(a,r_old)
    return E/N, N_accep/N

a = 0.9
nmax = 100000
tau = 1.3
X0 = [ MonteCarlo(a,nmax,tau) for i in range(30)]
X = [ x for x, _ in X0 ]
A = [ x for _, x in X0 ]
E, deltaE = ave_error(X)
A, deltaA = ave_error(A)
print(f"E = {E} +/- {deltaE}")
print(f"A = {A} +/- {deltaA}")
