import numpy as np

def potential(r):
    return -1. / np.sqrt(np.dot(r,r))

def psi(a, r):
    return np.exp(-a*np.sqrt(np.dot(r,r)))

def kinetic(a,r):
    return -0.5 * (a**2 - (2.*a)/np.sqrt(np.dot(r,r)))

def e_loc(a,r):
    return kinetic(a,r) + potential(r)
