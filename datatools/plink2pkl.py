import pandas as pd
import pickle
from datatools import imputesnp as isnp
from classes import SNPdataclass as snpc

# PLINK2PKL convert plink .raw file to pickle file format
#
# This function converts .raw plink file into .pkl file.
# It extracts all information from the .raw file and separates
# the genotype information from the rest. It saves all into
# a <outputFile>.pkl file.
#
# INPUTS:
# plinkRawFile - plink.raw file
# plinkBimFile - plink.bim file that associated with plinkRawFile
# plinkFamFile - plink.fam file that associated with plinkRawFile
# outputFile - name for output pickle file
#
# OUTPUTS:
# <outputFile>.pkl
#   The .pkl file uses an SNPdata class with the following fields:
#   - rsid: snp names
#   - data: genotype data
#   - chr: chromosome id
#   - loc: physical location
#   - pheno: sample's phenotype
#   - fid: sample's family id
#   - pid: sample id
#   - gender: sample gender

def plink2pkl(plinkRawFile,plinkBimFile,plinkFamFile,outputFile):

    # Creating headers for columns reading files into dataframes.
    bimHeader = ['chr', 'rsid', 'tmp1', 'loc', 'tmp2', 'tmp3']
    bdf = pd.read_csv(plinkBimFile, sep=r"\s+", names=bimHeader, engine='python')

    famHeader = ['fid', 'pid', 'tmp1', 'tmp2', 'gender', 'pheno']
    fdf = pd.read_csv(plinkFamFile, sep=r"\s+", names=famHeader, engine='python')

    rdf = pd.read_csv(plinkRawFile, sep=r"\s+", header=0, engine='python')

    # Imputing missing values from raw file dataframe.
    data = isnp.imputesnp(rdf[rdf.columns[6:]])

    # Structuring data to be saved into pickle format.
    SNPdata = snpc.SNPclass(data, bdf.rsid, bdf.chr, bdf.loc,
                        fdf.pheno-1, fdf.fid, fdf.pid, fdf.gender)

    # Save data to pickle file.
    final = open(outputFile, 'wb')
    pickle.dump(SNPdata, final)
    final.close()
