import pandas as pd
import numpy as np

def pathsim(pathind):

    sim = np.zeros((len(pathind), len(pathind)))

    for i in range(len(pathind)):
        for j in range((i+1),len(pathind)):
            sim[i][j] = len(list(set(pathind.iloc[i]) & set(pathind.iloc[j])))/min(len(pathind.iloc[i]), len(pathind.iloc[j]))


    sim = np.maximum(sim, sim.transpose())


    return sim
