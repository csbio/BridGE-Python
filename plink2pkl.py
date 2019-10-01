import pandas as pd
from imputesnp import imputesnp
from SNPdataclass import SNPclass
import pickle

# FUNCTION plink2mat(plinkRawFile,plinkBimFile,plinkFamFile,outputFileName)
#
# PLINK2MAT convert plink .raw file to mat-file format
#
# This function converts .raw plink file into .mat file.
# It extracts all information from the .raw file and separates
# the genotype information from the rest. It saves all into
# a <outputFile>.mat file.
#
# INPUTS:
# plinkRawFile - plink.raw file
# plinkBimFile - plink.bim file that associated with plinkRawFile
# plinkFamFile - plink.fam file that associated with plinkRawFile
# outputFile - name for output mat-file
#
# OUTPUTS:
# <outputFile>.mat
#   The .mat file consists a structure array SNPdata with
#   the following fields:
#   - rsid:snp names
#   - data:genotype data
#   - chr: chromosome id
#   - loc: physical location
#   - pheno: sample's phenotype
#   - fid: sample's family id
#   - pid: sample id
#   - gender: sample gender

def plink2pkl(plinkRawFile,plinkBimFile,plinkFamFile,outputFile):

    rdf = pd.read_csv(plinkRawFile, sep=r"\s*", header=0)

    bimHeader = ['chr', 'rsid', 'tmp1', 'loc', 'tmp2', 'tmp3']
    bdf = pd.read_csv(plinkBimFile, sep=r"\s*", names=bimHeader)

    famHeader = ['fid', 'pid', 'tmp1', 'tmp2', 'gender', 'pheno']
    fdf = pd.read_csv(plinkFamFile, sep=r"\s*", names=famHeader)

    data = imputesnp(rdf[rdf.columns[6:]])

    SNPdata = SNPclass(data, bdf.rsid, bdf.chr, bdf.loc,
                        fdf.pheno, fdf.fid, fdf.pid, fdf.gender)

    final = open(outputFile, 'wb')
    pickle.dump(SNPdata, final)
    final.close()
