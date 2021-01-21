from hydrogen  import *
from qmc_stats import *

def MonteCarlo(a,tau,nmax):
    E = 0.
    N = 0.
    accep_rate = 0.
    sq_tau = np.sqrt(tau)
    r_old = np.random.normal(loc=0., scale=1.0, size=(3))
    d_old = drift(a,r_old)
    d2_old = np.dot(d_old,d_old)
    psi_old = psi(a,r_old)
    for istep in range(nmax):
        chi = np.random.normal(loc=0., scale=1.0, size=(3))
        r_new = r_old + tau * d_old + sq_tau * chi
        d_new = drift(a,r_new)
        d2_new = np.dot(d_new,d_new)
        psi_new = psi(a,r_new)
        # Metropolis
        prod = np.dot((d_new + d_old), (r_new - r_old))
        argexpo = 0.5 * (d2_new - d2_old)*tau + prod
        q = psi_new / psi_old
        q = np.exp(-argexpo) * q*q
        if np.random.uniform() < q:
            accep_rate += 1.
            r_old = r_new
            d_old = d_new
            d2_old = d2_new
            psi_old = psi_new
        N += 1.
        E += e_loc(a,r_old)
    return E/N, accep_rate/N


a = 0.9
nmax = 100000
tau = 1.0
X = [MonteCarlo(a,tau,nmax) for i in range(30)]
E, deltaE = ave_error([x[0] for x in X])
A, deltaA = ave_error([x[1] for x in X])
print(f"E = {E} +/- {deltaE}\nA = {A} +/- {deltaA}")
