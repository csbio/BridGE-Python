import pandas as pd
import numpy as np
import pickle
from corefuns import bpmsim as bpmsim
from corefuns import pathsim as pathsim
from scipy import sparse

def check_BPM_WPM_redundancy(fdrBPM,fdrWPM,fdrPATH,bpmindfile,FDRcut):
    pklin = open(bpmindfile,"rb")
    bpmind = pickle.load(pklin)
    pklin.close()

    BPM_group, BPM_nosig_noRD = [], []
    WPM_group, WPM_nosig_noRD = [], []
    PATH_group, PATH_nosig_noRD = [], []

    start = 0.05
    increment = 0.05

    # Adding increment to FDRcut for inclusive arange.
    for fdrcut in np.arange(start, FDRcut+increment, increment):
        ind = np.array(fdrBPM[fdrBPM<=fdrcut].dropna().index)
        nnz = (ind<=len(bpmind.bpm['size'])).sum()

        ind1 = ind[0:nnz]
        ind2 = ind[nnz+1:]

        ind = np.array(fdrWPM[fdrWPM<=fdrcut].dropna().index)
        nnz = (ind<=len(bpmind.wpm['size'])).sum()

        ind3 = ind[0:nnz]
        ind4 = ind[nnz+1:]
        #
        # ind1.rename(columns={"bpm2": "ind1"}, inplace = True)
        # ind2.rename(columns={"bpm2": "ind2"}, inplace = True)
        # ind3.rename(columns={"wpm2": "ind3"}, inplace = True)
        # ind4.rename(columns={"wpm2": "ind4"}, inplace = True)

        group, BPM_sim, noRD = [0]*4, [0]*4, [0]*4

        # Should ind have more than one row
        if len(ind1) > 1:
            bpm_ind1_ind1 = bpmind.bpm['ind1'].iloc[(ind1)%15]
            bpm_ind2_ind1 = bpmind.bpm['ind2'].iloc[(ind1)%15]

            # TODO: bpmsim incomplete
            BPM_sim[0] = bpmsim.bpmsim(bpm_ind1_ind1, bpm_ind2_ind1, bpm_ind1_ind1, bpm_ind2_ind1)

            TTT = (BPM_sim[0]>=0.25).astype(int)
            noRD[0], group[0] = sparse.csgraph.connected_components(TTT)
            a0=len(list(set(group[0])))

        # Has exactly one row
        elif len(ind1) == 1:
            a0=1
            group[0] = [1]
        else:
            a0=0
            group[0] = []

        # Should ind have more than one row
        if len(ind2) > 1:
            bpm_ind1_ind2 = bpmind.bpm['ind1'].iloc[(ind2)%15]
            bpm_ind2_ind2 = bpmind.bpm['ind2'].iloc[(ind2)%15]
            # print('====  bpm_ind1_ind2 ====')
            # print('====  bpm_ind2_ind2 ====')
            # print(bpm_ind1_ind2)
            # print(bpm_ind2_ind2)
            # input()

            # TODO: bpmsim incomplete
            BPM_sim[1] = bpmsim.bpmsim(bpm_ind1_ind2, bpm_ind2_ind2, bpm_ind1_ind2, bpm_ind2_ind2)

            TTT = (BPM_sim[0]>=0.25).astype(int)
            noRD[0], group[0] = sparse.csgraph.connected_components(TTT)
            b0=len(list(set(group[1])))


        # Has exactly one row
        elif len(ind2) == 1:
            b0=1
            group[1] = [1]
        else:
            b0=0
            group[1] = []

        # Should ind have more than one row
        if len(ind3) > 1:
            wpm_ind_ind3 = bpmind.wpm['ind'].iloc[(ind3)%6]
            # print('====  wpm_ind_ind3 ====')
            # print(wpm_ind_ind3)
            # input()

            # TODO: bpmsim incomplete
            BPM_sim[2] = bpmsim.bpmsim(wpm_ind_ind3, wpm_ind_ind3, wpm_ind_ind3, wpm_ind_ind3)

            TTT = (BPM_sim[0]>=0.25).astype(int)
            noRD[0], group[0] = sparse.csgraph.connected_components(TTT)
            c0=len(list(set(group[2])))


        # Has exactly one row
        elif len(ind3) == 1:
            c0=1
            group[2] = [1]
        else:
            c0=0
            group[2] = []

        # Should ind have more than one row
        if len(ind4) > 1:
            wpm_ind_ind4 = bpmind.wpm['ind'].iloc[(ind4)%6]
            # print('====  wpm_ind_ind4 ====')
            # print(wpm_ind_ind4)
            # input()

            # TODO: bpmsim incomplete
            BPM_sim[3] = bpmsim.bpmsim(wpm_ind_ind4, wpm_ind_ind4, wpm_ind_ind4, wpm_ind_ind4)

            TTT = (BPM_sim[0]>=0.25).astype(int)
            noRD[0], group[0] = sparse.csgraph.connected_components(TTT)
            d0=len(list(set(group[3])))


        # Has exactly one row
        elif len(ind4) == 1:
            d0=1
            group[3] = [1]
        else:
            d0=0
            group[3] = []

        # print("\n==== First group print ====")
        # print(group)
        # input()

        g1 = len(list(set(group[0])))
        g2 = len(list(set(group[1])))
        g3 = len(list(set(group[2])))
        g4 = len(list(set(group[3])))

        # print("== printing g1 - g4 ==")
        # print(g1)
        # print(g2)
        # print(g3)
        # print(g4)
        # input()

        BPM_group.append(list(group[0]) + list(map(lambda x : x + g1, list(group[1]))))
        BPM_nosig_noRD.append(a0 + b0)

        WPM_group.append(list(group[2]) + list(map(lambda x : x + g3, list(group[3]))))
        WPM_nosig_noRD.append(c0 + d0)

        ind = np.array(fdrPATH[fdrPATH<=fdrcut].dropna().index)
        nnz = (ind<=len(bpmind.wpm['size'])).sum()

        ind1 = ind[0:nnz]
        ind2 = ind[nnz+1:]

        group, PATH_sim, noRD = [0]*2, [0]*2, [0]*2

        # Should ind have more than one row
        if len(ind1) > 1:
            wpm_ind_ind1 = bpmind.wpm['ind'].iloc[(ind1)%6]

            PATH_sim[0] = pathsim.pathsim(wpm_ind_ind1)
            TTT = (PATH_sim[0]>=0.25).astype(int)
            noRD[0], group[0] = sparse.csgraph.connected_components(TTT)
            a0=len(list(set(group[0])))

        # Has exactly one row
        elif len(ind1) == 1:
            a0=1
            group[0] = [1]
        else:
            a0=0
            group[0] = []

        # Should ind have more than one row
        if len(ind2) > 1:
            wpm_ind_ind2 = bpmind.wpm['ind'].iloc[(ind2.index)%6]

            PATH_sim[1] = pathsim.pathsim(wpm_ind_ind2)
            TTT = (PATH_sim[1]>=0.25).astype(int)
            noRD[1], group[1] = sparse.csgraph.connected_components(TTT)
            b0=len(list(set(group[1])))

        # Has exactly one row
        elif len(ind2) == 1:
            b0=1
            group[1] = [1]
        else:
            b0=0
            group[1] = []

        g1 = len(list(set(group[0])))
        g2 = len(list(set(group[1])))

        # print("=== g1 & g2 ===")
        # print(g1)
        # print(g2)

        PATH_group.append(list(group[0]) + list(map(lambda x : x + g1, list(group[1]))))
        PATH_nosig_noRD.append(a0 + b0)

        # print('==== Group Testing ====')
        # print(BPM_group)
        # print(BPM_nosig_noRD)
        # print(WPM_group)
        # print(WPM_nosig_noRD)
        # print(PATH_group)
        # print(PATH_nosig_noRD)
        # print('=======================')

    return BPM_nosig_noRD,WPM_nosig_noRD,PATH_nosig_noRD,BPM_group,WPM_group,PATH_group
