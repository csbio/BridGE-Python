import pandas as pd
import pickle
from classes import genesetdataclass as gsc
import numpy as np
import math

# MSIGDB2PKL convert MsigDB gene set file (.gmt) to pickle file (Python pkl).
#
# INPUTS:
#   symbolsFile: MsigDB gene set file using gene symbols (.symbols.gmt).
#   entrezFile: MsigDB gene set file using gene entrez ids (.entrez.gmt).
#
# OUTPUTS:
#   <symbolsFile>.pkl -This pickle file uses a geneset class with fields:
#       geneset.entrezids - gene {symbol: entrezID} lookup dictionary
#       geneset.gpmatrix - gene pathway binary dataframe

def msigdb2pkl(symbolsFile, entrezFile):
    # Reading files into dataframes.
    sdf = pd.read_csv(symbolsFile, sep=r"\s+", header=None, engine='python')
    edf = pd.read_csv(entrezFile, sep=r"\s+", header=None, engine='python')

    # Accumulator list, symbol list, and entrezID list
    acclist, symlist, idlist = [], [], []

    # Iterating over (3 to last) columns, building pathway, symbol and ID lists.
    for i in range(2, sdf.shape[1]):
        acclist += tuple(zip(sdf[0], sdf[i]))
        symlist += list(sdf[i])
        idlist += list(edf[i])

    # Creating dictionary for easy lookup of entrezID by symbol.
    symboldict = dict(zip(symlist, pd.to_numeric(idlist)))

    # Removing extraneous values from each list.
    symboldict = filter(lambda k:   not math.isnan(k[1]), symboldict.items())
    acclist = filter(lambda k: k[1] != None, acclist)
    symlist = filter(lambda k: k != None, symlist)
    idlist = filter(lambda k: not math.isnan(k), idlist)

    # Set conversion to get unique entries, then numpy conversion for storage.
    idlist = np.array(list(set(idlist)))
    symlist = np.array(list(set(symlist)))
    pwlist = np.array(list(sdf[0]))
    symboldict = dict(set(symboldict))

    # Building gene pathway matrix (gpm) of appropriate size, and adding labels.
    gpm = pd.DataFrame(np.zeros((len(symlist), sdf.shape[0])),
                                index=symlist, columns=pwlist, dtype = int)

    # Set conversion for unique values to iterate over.
    accset = set(acclist)
    for tup in accset:
        # Marking true in gpm where symbol is in pathway.
        gpm[tup[0]][tup[1]] = 1

    # Converting data to pickle storage file with geneset class.
    geneset = gsc.genesetclass(symboldict, gpm)
    symbolsFile = symbolsFile.replace(".symbols.gmt", ".pkl")
    final = open(symbolsFile, 'wb')
    pickle.dump(geneset, final)
    final.close()
