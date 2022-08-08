import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform

def bpmsim(BPMind1x, BPMind2x, BPMind1y, BPMind2y):

    t1 = len(BPMind1x)
    t2 = len(BPMind1y)

    test = BPMind1x.equals(BPMind1y)==True and BPMind2x.equals(BPMind2y)==True

    BPMind1x = np.transpose(np.tile(BPMind1x,(t2,1)))
    BPMind2x = np.transpose(np.tile(BPMind2x,(t2,1)))

    BPMind1y = np.tile(BPMind1y,(t1,1))
    BPMind2y = np.tile(BPMind2y,(t1,1))

    if test:
        BPMind1x = cellsquareform(BPMind1x)
        BPMind2x = cellsquareform(BPMind2x)
        BPMind1y = cellsquareform(BPMind1y)
        BPMind2y = cellsquareform(BPMind2y)

    result_shape = BPMind1x.shape
    BPMind1x = BPMind1x.flatten()
    BPMind2x = BPMind2x.flatten()
    BPMind1y = BPMind1y.flatten()
    BPMind2y = BPMind2y.flatten()

    mm = list(map(lambda x,y,z,r : min(len(x)*len(y),len(z)*len(r)),BPMind1x,BPMind2x,BPMind1y,BPMind2y))

    ind1 = list(map(lambda x,y : list(set(x) & set(y)), BPMind1x,BPMind1y))
    ind2 = list(map(lambda x,y : list(set(x) & set(y)), BPMind2x,BPMind2y))

    simtmp1 = list(map(lambda x,y : len(x)*len(y), ind1,ind2))
    simtmp1 = np.divide(simtmp1, mm)

    ind1 = list(map(lambda x,y : list(set(x) & set(y)), BPMind1x,BPMind2y))
    ind2 = list(map(lambda x,y : list(set(x) & set(y)), BPMind2x,BPMind1y))

    simtmp2 = list(map(lambda x,y : len(x)*len(y), ind1,ind2))
    simtmp2 = np.divide(simtmp2, mm)
    
    sim = np.maximum(simtmp1,simtmp2)
    sim = np.reshape(sim,result_shape)

    if test:
        sim = squareform(sim)

    return sim

def cellsquareform(Y):
    n = len(Y[0])
    Z = Y[np.nonzero(np.logical_and(np.triu(np.ones((n, n), dtype=bool), 1), Y))]
    Z = np.transpose(Z)
    return Z
