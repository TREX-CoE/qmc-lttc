#!/usr/bin/env python3

import numpy as np
from hydrogen import e_loc, psi

interval = np.linspace(-5,5,num=50)
delta = (interval[1]-interval[0])**3

r = np.array([0.,0.,0.])

for a in [0.1, 0.2, 0.5, 0.9, 1., 1.5, 2.]:
    E    = 0.
    norm = 0.

    for x in interval:
        r[0] = x
        for y in interval:
            r[1] = y
            for z in interval:
                r[2] = z

                w = psi(a,r)
                w = w * w * delta

                E    += w * e_loc(a,r)
                norm += w 

    E = E / norm
    print(f"a = {a} \t E = {E}")
