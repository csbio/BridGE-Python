from scipy.stats import hypergeom
import numpy as np
import matplotlib.pyplot as plt

# HYGETEST computes Hypergeometric cumulative distribution for
# the given inputs.
#
# INPUTS:
# n - is the population size
# d - is the number of draws
# k - is the number of successes
# m - is the number of success states in the population
#
# OUTPUTS:
# logpv - negative log10(p-value)
# pv - p-value

def hygetest(n,d,k,m):
    pv = 1 - hypergeom.cdf(k-1,n,m,d)
    #logpv = -np.log10(pv)
    return pv
