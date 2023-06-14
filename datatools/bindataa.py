import pickle
import numpy as np

# BINDATAA binarize 012 format SNP data based on dominant/recessive assumptions.
#
# INPUTS:
#   project_dir: directory of all the project files
#   dataFile - name of the data file. This .mat file consists
#       a structure array SNPdata with the following fields:
#       - rsid:snp names
#       - data:genotype data
#       - chr: chromosome id
#       - loc: physical location
#       - pheno: sample's phenotype
#       - fid: sample's family id
#       - pid: sample id
#       - gender
#   expr - flag used to designate dominant ('d'/'D') or recessive ('r'/'R')
#
# OUTPUTS:
#   a pickle file SNPdataA(D or R).pkl

def bindataa(project_dir,dataFile, expr):
    # Checking expression flag to proceed as dominant or recessive (D or R).
    if expr == 'r' or expr == 'R':
        # If recessive, set 1s to 0s, 2s to 1s, and set appropriate filename.
        set1, set2 = 0, 1
        filename = project_dir+'/intermediate/SNPdataAR.pkl'
    elif expr == 'd' or expr == 'D':
        # If dominant, set 1s to 1s, 2s to 1s, and set appropriate filename.
        set1, set2 = 1, 1
        filename = project_dir+'/intermediate/SNPdataAD.pkl'
    else:
        # Default case where expression provided was neither D or R
        print("Provide 'd'/'D' or 'r'/'R' to designate dominant/recessive.")
        return

    # Reading in pickle datafile
    pklin = open(dataFile,"rb")
    SNPdata = pickle.load(pklin)
    pklin.close()

    # Replacing 1s and 2s with corresponding values set above as D or R.
    SNPdata.data = SNPdata.data.replace(to_replace=(1), value=set1)
    SNPdata.data = SNPdata.data.replace(to_replace=(2), value=set2)

    # Saving updated SNPdata in output pickle file.
    final = open(filename, 'wb')
    pickle.dump(SNPdata, final)
    final.close()
