import pickle
from itertools import combinations
import numpy as np
import pandas as pd
from classes import bpmindclass as bpmc
import datetime
import sys

def bpmind(snpPathwayFile):

    # Reading in pickle datafile
    print('Hello from bpmind')
    print(datetime.datetime.now())
    pklin = open(snpPathwayFile,"rb")
    snpset = pickle.load(pklin)
    pklin.close()

    # Retrieving pathways list from snpset
    pathways = snpset.pathways
    snpmat = snpset.spmatrix
    snpfile = snpset.geneset
    print('Reading pathways completed')
    print(datetime.datetime.now())
    sys.stdout.flush()


    # Finding all possible combinations of pairs for pathway names and sizes.
    comb = np.array(list(combinations(pathways, 2)))
    combnames = np.array(list(combinations(pathways.index, 2)))

    WPMind = []
    for column in snpmat:
        nonzeros = list(np.nonzero(snpmat[column].values)[0])
        WPMind.append(nonzeros)

    print('Starting the for loop, # of outer loops:'+str(len(snpmat.columns)))
    print(datetime.datetime.now())
    sys.stdout.flush()
    BPMind1, BPMind2, ind1size, ind2size = [], [], [], []
    for i in range(len(snpmat.columns)):
        p1 = snpmat.iloc[:,i].to_numpy()
        p1[p1>1] = 1
        if i % 10 == 0:
            print('in loop:'+str(i))
            print(datetime.datetime.now())
            sys.stdout.flush()
        for j in range(i+1, len(snpmat.columns)):
            p2 = snpmat.iloc[:,j].to_numpy()
            p2[p2>1] = 1
            ind1, ind2 = [], []
            #for k in range(len(snpmat)):
            #    if (snpmat.iloc[k][i]==1 and snpmat.iloc[k][j]!=1):
            #        ind1.append(k)
            #    if (snpmat.iloc[k][i]!=1 and snpmat.iloc[k][j]==1):
            #        ind2.append(k)
            d1 = p1 - p2
            ind1 = np.where(d1==1)
            ind1 = ind1[0].tolist()
            d2 = p2 - p1
            ind2 = np.where(d2==1)
            ind2 = ind2[0].tolist()
            ind1size.append(len(ind1))
            ind2size.append(len(ind2))
            BPMind1.append(ind1)
            BPMind2.append(ind2)

    print(len(pathways))
    sys.stdout.flush()


    # Getting between pathway sizes by multiplying combination available pairs.
    if (len(pathways) > 1):
        size = np.array(ind1size)*np.array(ind2size)
        # Orienting bpm/wpm data and converting to dataframes.
        bpmdata = {'path1names': combnames[:,0], 'ind1size': ind1size, 'ind1': BPMind1,
                'path2names': combnames[:,1], 'ind2size': ind2size, 'ind2': BPMind2,
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
