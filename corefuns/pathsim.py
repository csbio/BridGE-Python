import pandas as pd
import numpy as np



# pathsim() computes similarity between pathways by dividing number of shared SNPs by the size of the smaller pathway 
#
# INPUTS:
#   pathind: A dataframe contating all the pathways needed to be compared 
#
# OUTPUTS:
#   returns a matrix with entries equal to similarity between the pathway in the row and the pathway in the column
# 



def pathsim(pathind):

    sim = np.zeros((len(pathind), len(pathind)))

    for i in range(len(pathind)):
        for j in range((i+1),len(pathind)):
            sim[i][j] = len(list(set(pathind.iloc[i]) & set(pathind.iloc[j])))/min(len(pathind.iloc[i]), len(pathind.iloc[j]))


    sim = np.maximum(sim, sim.transpose())


    return sim
