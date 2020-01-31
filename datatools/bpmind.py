import pickle
from itertools import combinations
import numpy as np
import pandas as pd
from classes import bpmindclass as bpmc

def bpmind(snpPathwayFile):

    # Reading in pickle datafile
    pklin = open(snpPathwayFile,"rb")
    snpset = pickle.load(pklin)
    pklin.close()

    # Retrieving pathways list from snpset
    pathways = snpset.pathways
    snpmat = snpset.spmatrix
    snpfile = snpset.geneset

    # Finding all possible combinations of pairs for pathway names and sizes.
    comb = np.array(list(combinations(pathways, 2)))
    combnames = np.array(list(combinations(pathways.index, 2)))

    WPMind = []
    for column in snpmat:
        nonzeros = list(np.nonzero(snpmat[column].to_numpy())[0])
        WPMind.append(nonzeros)

    BPMind1, BPMind2 = [], []
    for i in range(len(snpmat.columns)):
        for j in range(i+1, len(snpmat.columns)):
            ind1, ind2 = [], []
            for k in range(len(snpmat)):
                if (snpmat.iloc[k][i]==1 and snpmat.iloc[k][j]!=1):
                    ind1.append(k)
                if (snpmat.iloc[k][i]!=1 and snpmat.iloc[k][j]==1):
                    ind2.append(k)
            BPMind1.append(ind1)
            BPMind2.append(ind2)

    # Getting between pathway sizes by multiplying combination available pairs.
    if (len(pathways) > 1):
        size = comb[:,0]*comb[:,1]
        # Orienting bpm/wpm data and converting to dataframes.
        bpmdata = {'path1names': combnames[:,0], 'ind1size': comb[:,0], 'ind1': BPMind1,
                'path2names': combnames[:,1], 'ind2size': comb[:,1], 'ind2': BPMind2,
                'size': size}

    else:
        bpmdata = {'path1names': [], 'ind1size': [],
                'path2names': [], 'ind2size': [],'size': []}

    bpm = pd.DataFrame(bpmdata)

    wpmdata = {'pathway': pathways.index, 'indsize': pathways.values, 'ind': WPMind,
                'size': (pathways.values*pathways.values)-pathways.values}
    wpm = pd.DataFrame(wpmdata)

    # Reading bpm and wpm models into bpmind class for pickle storage.
    bpmobj = bpmc.bpmindclass(bpm, wpm)

    # Saving bpmind data to pickle file.
    final = open('data/BPMind.pkl', 'wb')
    pickle.dump(bpmobj, final)
    final.close()
