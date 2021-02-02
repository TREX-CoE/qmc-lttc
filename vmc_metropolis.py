#!/usr/bin/env python3

from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a,nmax,dt):
    sq_dt = np.sqrt(dt)

    energy  = 0.
    N_accep = 0

    r_old   = np.random.normal(loc=0., scale=1.0, size=(3))
    d_old   = drift(a,r_old)
    d2_old  = np.dot(d_old,d_old)
    psi_old = psi(a,r_old)

    for istep in range(nmax):
        chi = np.random.normal(loc=0., scale=1.0, size=(3))

        energy += e_loc(a,r_old)

        r_new   = r_old + dt * d_old + sq_dt * chi
        d_new   = drift(a,r_new)
        d2_new  = np.dot(d_new,d_new)
        psi_new = psi(a,r_new)

        # Metropolis
        prod    = np.dot((d_new + d_old), (r_new - r_old))
        argexpo = 0.5 * (d2_new - d2_old)*dt + prod

        q = psi_new / psi_old
        q = np.exp(-argexpo) * q*q

        if np.random.uniform() <= q:
            N_accep += 1

            r_old   = r_new
            d_old   = d_new
            d2_old  = d2_new
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
