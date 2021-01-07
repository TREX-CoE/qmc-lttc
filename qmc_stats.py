from math import sqrt
def ave_error(arr):
    M = len(arr)
    assert (M>1)
    average = sum(arr)/M
    variance = 1./(M-1) * sum( [ (x - average)**2 for x in arr ] )
    return (average, sqrt(variance/M))
