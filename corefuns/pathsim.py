import pandas as pd
import numpy as np

def pathsim(pathind):

    sim = np.zeros((len(pathind), len(pathind)))
    # print("=== Inside pathsim ===")
    # print("--- sim zeros ---")
    # print(sim)
    # input()

    for i in range(len(pathind)):
        for j in range((i+1),len(pathind)):
            sim[i][j] = len(list(set(pathind[i]) & set(pathind[j])))/min(len(pathind[i]), len(pathind[j]))

    # print("--- sim ratio ---")
    # print(sim)
    # input()

    sim = np.maximum(sim, sim.transpose())

    # print("--- sim final ---")
    # print(sim)
    # input()

    return sim
