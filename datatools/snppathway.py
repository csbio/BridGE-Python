import pickle
import numpy as np
from classes import snpsetclass as snps
import pandas as pd


# SNPPathway creates a snp-to-pathway mapping.
#
# INPUTS:
#   dataFile: Path to the file with genotype data in the Pickle format.
#   sgmFile: SNP to gene mapping file in the Pickle format.
#   genesets: Gene-set file in pickle format.
#   minPath: minimum size for a pathway to be in the mapping.
#   maxPath: maximum size for a pathway to be in the mapping.
#
# OUTPUTS:
#   snp_pathway_min<minPath>_max<maxPath>.pkl -This pickle file contains a snpset class object with following fields:
#       - pathways: List of pathway names
#       - spmatrix: Matrix of snp-pathway mapping (Numpy 2d array)
#       - geneset: Path to geneset file in .pkl format.


def snppathway(dataFile,sgmFile,genesets,minPath,maxPath):

    # Loading pickle files into objects
    pklin = open(dataFile,"rb")
    SNPdata = pickle.load(pklin)
    pklin.close()
    pklin = open(sgmFile,"rb")
    sgm = pickle.load(pklin)
    pklin.close()
    pklin = open(genesets,"rb")
    geneset = pickle.load(pklin)
    pklin.close()
    # Dropping the set difference of the genes in the gpmatrix and sgmatrix.
    genedroplist = list(set(geneset.gpmatrix.index).difference(set(sgm.columns)))
    sgpm = geneset.gpmatrix.drop(genedroplist, axis=0)
    orig_order = list(sgpm.columns)
    sgpm = sgpm.replace(to_replace=(2), value=1)
    sgpm = sgpm.sort_index(axis=1)

    
    # sort sgm rsids based on the SNPdata
    rsg = sgm.index.values
    rss = SNPdata.rsid
    tmp_sgm = np.zeros((rss.shape[0],sgm.shape[1]))
    rs_order = []
    #for rs in rss:
    for i in range(rss.shape[0]):
        rs = rss[i]
        tmp = np.where(rsg==rs)
        if tmp[0].shape[0] > 0:
            #rs_order.append(tmp[0][0])
            tmp_sgm[i,:] = sgm.iloc[tmp[0][0],:]

    #sgm = sgm.iloc[rs_order]
    sgm = pd.DataFrame(tmp_sgm,index=rss,columns=sgm.columns)

    # Only keeping the columns in both the snp-gene pathway
    testm = sgm[sgpm.index]
    # Matrix multiply dataframes to associate rsids with pathways through genes.
    final = testm.dot(sgpm)
    final = final.reindex(columns=orig_order)

    # Summing columns of pathway to known snps and filtering those not in range.
    fnlsum = final.astype(bool).sum()
    validpwysfnl = fnlsum[fnlsum.apply(lambda x: filter_val(minPath, maxPath, x))]
    spm = final[validpwysfnl.index]
    spm = spm.replace(to_replace='^[1-9][0-9]*$', value=1,regex=True)
    pathways = spm.astype(bool).sum()

    # Preparing data and filename for pickle storage.
    snpset = snps.snpsetclass(pathways, spm, genesets)
    outfilename = "data/snp_pathway_min"+str(minPath)+"_max"+str(maxPath)+".pkl"

    # Saving data to pickle file.
    final = open(outfilename, 'wb')
    pickle.dump(snpset, final)
    final.close()

    # Returning the name of the output file to be used by other modules.
    return outfilename

# Filter to keep pathways with number of snps between lower and upper range.
def filter_val(lower, upper, length):
    return ((lower <= length) and (length <= upper))
