from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a,nmax,tau):
    energy = 0.
    N_accep = 0
    r_old = np.random.uniform(-tau, tau, (3))
    psi_old = psi(a,r_old)
    for istep in range(nmax):
        r_new = r_old + np.random.uniform(-tau,tau,(3))
        psi_new = psi(a,r_new)
        ratio = (psi_new / psi_old)**2
        v = np.random.uniform()
        if v <= ratio:
            N_accep += 1
            r_old = r_new
            psi_old = psi_new
        energy += e_loc(a,r_old)
    return energy/nmax, N_accep/nmax

# Run simulation
a = 0.9
nmax = 100000
tau = 1.3
X0 = [ MonteCarlo(a,nmax,tau) for i in range(30)]

# Energy
X = [ x for (x, _) in X0 ]
E, deltaE = ave_error(X)
print(f"E = {E} +/- {deltaE}")

# Acceptance rate
X = [ x for (_, x) in X0 ]
A, deltaA = ave_error(X)
print(f"A = {A} +/- {deltaA}")
