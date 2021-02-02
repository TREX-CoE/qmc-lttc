#!/usr/bin/env python3

from math import sqrt
def ave_error(arr):
    M = len(arr)
    assert(M>0)

    if M == 1:
        average = arr[0]
        error   = 0.

    else:
        average = sum(arr)/M
        variance = 1./(M-1) * sum( [ (x - average)**2 for x in arr ] )
        error = sqrt(variance/M)

    return (average, error)
