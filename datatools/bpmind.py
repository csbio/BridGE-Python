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
    print(pathways)

    # Finding all possible combinations of pairs for pathway names and sizes.
    comb = np.array(list(combinations(pathways, 2)))
    combnames = np.array(list(combinations(pathways.index, 2)))

    # Getting between pathway sizes by multiplying combination pairs.
    print(comb.shape)
    size = comb[:,0]*comb[:,1]

    # Orienting bpm/wpm data and converting to dataframes.
    bpmdata = {'path1names': combnames[:,0], 'ind1size': comb[:,0],
            'path2names': combnames[:,1], 'ind2size': comb[:,1],'size': size}
    bpm = pd.DataFrame(bpmdata)

    wpmdata = {'names': pathways.index, 'indsize': pathways.values,
                'size': (pathways.values*pathways.values)-pathways.values}
    wpm = pd.DataFrame(wpmdata)

    # Reading bpm and wpm models into bpmind class for pickle storage.
    bpmobj = bpmc.bpmindclass(bpm, wpm)

    # Saving bpmind data to pickle file.
    final = open('data/BPMind.pkl', 'wb')
    pickle.dump(bpmobj, final)
    final.close()
