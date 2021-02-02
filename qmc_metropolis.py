#!/usr/bin/env python3

from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a,nmax,dt):
    energy  = 0.
    N_accep = 0

    r_old = np.random.uniform(-dt, dt, (3))
    psi_old = psi(a,r_old)

    for istep in range(nmax):
        energy += e_loc(a,r_old)

        r_new = r_old + np.random.uniform(-dt,dt,(3))
        psi_new = psi(a,r_new)

        ratio = (psi_new / psi_old)**2

        if np.random.uniform() <= ratio:
            N_accep += 1

            r_old   = r_new
            psi_old = psi_new

    return energy/nmax, N_accep/nmax

# Run simulation
a    = 1.2
nmax = 100000
dt   = 1.0

X0 = [ MonteCarlo(a,nmax,dt) for i in range(30)]

# Energy
X = [ x for (x, _) in X0 ]
E, deltaE = ave_error(X)
print(f"E = {E} +/- {deltaE}")

# Acceptance rate
X = [ x for (_, x) in X0 ]
A, deltaA = ave_error(X)
print(f"A = {A} +/- {deltaA}")
