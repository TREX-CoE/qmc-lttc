import numpy as np

def potential(r):
    r2 = np.dot(r,r)
    assert (r2 > 0.)
    distance = np.sqrt(r2)
    
    return -1. / distance

def psi(a, r):
    return np.exp(-a*np.sqrt(np.dot(r,r)))

def kinetic(a,r):
    r2 = np.dot(r,r)
    assert (r2 > 0.)
    distance = np.sqrt(r2)

    return a * (1./distance - 0.5 * a)

def e_loc(a,r):
    return kinetic(a,r) + potential(r)

def drift(a,r):
   ar_inv = -a/np.sqrt(np.dot(r,r))
   return r * ar_inv
