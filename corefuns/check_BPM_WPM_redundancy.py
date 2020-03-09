import pandas as pd
import numpy as np
import pickle
from corefuns import bpmsim as bpmsim
from corefuns import pathsim as pathsim

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

        # Double check test cases to ensure slicing conversion is accurate between languages.
        ind = fdrBPM[fdrBPM<=fdrcut]
        indt = ind<=len(bpmind.bpm['size'])
        nnz = indt.sum().values[0]
        ind1 = ind[0:nnz]
        ind2 = ind[nnz+1:].dropna()

        ind = fdrWPM[fdrWPM<=fdrcut]
        indt = ind<=len(bpmind.wpm['size'])
        nnz = indt.sum().values[0]
        ind3 = ind[0:nnz]
        ind4 = ind[nnz+1:].dropna()

        ind1.rename(columns={"bpm2": "ind1"}, inplace = True)
        ind2.rename(columns={"bpm2": "ind2"}, inplace = True)
        ind3.rename(columns={"wpm2": "ind3"}, inplace = True)
        ind4.rename(columns={"wpm2": "ind4"}, inplace = True)

        group, BPM_sim = [0]*4, [0]*4

        # Should ind have more than one row
        if len(ind1.index) > 1:
            bpm_ind1_ind1 = bpmind.bpm['ind1'].iloc[(ind1.index)%15]
            bpm_ind2_ind1 = bpmind.bpm['ind2'].iloc[(ind1.index)%15]
            # print('====  ind1.shape[1] > 1  ====')
            # print(bpm_ind1_ind1)
            # print(bpm_ind2_ind1)
            # input()

            # TODO: bpmsim incomplete
            BPM_sim[0] = bpmsim.bpmsim(bpm_ind1_ind1, bpm_ind2_ind1, bpm_ind1_ind1, bpm_ind2_ind1)

            a0=len(list(set(group[0])))
            a0=0

        # Has exactly one row
        elif len(ind1.index) == 1:
            a0=1
            group[0] = [1]
        else:
            a0=0
            group[0] = []

        # Should ind have more than one row
        if len(ind2.index) > 1:
            bpm_ind1_ind2 = bpmind.bpm['ind1'].iloc[(ind2.index)%15]
            bpm_ind2_ind2 = bpmind.bpm['ind2'].iloc[(ind2.index)%15]
            # print('====  ind2.shape[1] > 1  ====')
            # print(bpm_ind1_ind2)
            # print(bpm_ind2_ind2)
            # input()

            # TODO: bpmsim incomplete
            BPM_sim[1] = bpmsim.bpmsim(bpm_ind1_ind2, bpm_ind2_ind2, bpm_ind1_ind2, bpm_ind2_ind2)
            b0=len(list(set(group[1])))
            b0=0

        # Has exactly one row
        elif len(ind2.index) == 1:
            b0=1
            group[1] = [1]
        else:
            b0=0
            group[1] = []

        # Should ind have more than one row
        if len(ind3.index) > 1:
            wpm_ind_ind3 = bpmind.wpm['ind'].iloc[(ind3.index)%6]
            # print('====  ind3.shape[1] > 1  ====')
            # print(wpm_ind_ind3)
            # input()

            # TODO: bpmsim incomplete
            BPM_sim[2] = bpmsim.bpmsim(wpm_ind_ind3, wpm_ind_ind3, wpm_ind_ind3, wpm_ind_ind3)
            c0=len(list(set(group[2])))
            c0=0

        # Has exactly one row
        elif len(ind3.index) == 1:
            c0=1
            group[2] = [1]
        else:
            c0=0
            group[2] = []

        # Should ind have more than one row
        if len(ind4.index) > 1:
            wpm_ind_ind4 = bpmind.wpm['ind'].iloc[(ind4.index)%6]
            # print('====  ind4.shape[1] > 1  ====')
            # print(wpm_ind_ind4)
            # input()

            # TODO: bpmsim incomplete
            BPM_sim[3] = bpmsim.bpmsim(wpm_ind_ind4, wpm_ind_ind4, wpm_ind_ind4, wpm_ind_ind4)
            d0=len(list(set(group[3])))
            d0=0

        # Has exactly one row
        elif len(ind4.index) == 1:
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

        BPM_group.append(group[0] + list(map(lambda x : x + g1, group[1])))
        BPM_nosig_noRD.append(a0 + b0)

        WPM_group.append(group[2] + list(map(lambda x : x + g3, group[3])))
        WPM_nosig_noRD.append(c0 + d0)

        ind = fdrPATH[fdrPATH<=fdrcut]
        indt = ind<=len(bpmind.wpm['size'])
        nnz = indt.sum().values[0]
        ind1 = ind[0:nnz]
        ind2 = ind[nnz+1:].dropna()

        group, PATH_sim = [0]*2, [0]*2

        # Should ind have more than one row
        if len(ind1.index) > 1:
            wpm_ind_ind1 = bpmind.wpm['ind'].iloc[(ind1.index)%6]
            # print('====  ind1.shape[1] > 1  ====')
            # print(wpm_ind_ind1)
            # input()

            # TODO: pathsim incomplete
            PATH_sim[0] = pathsim.pathsim(wpm_ind_ind1)
            a0=len(list(set(group[0])))
            a0=0

        # Has exactly one row
        elif len(ind1.index) == 1:
            a0=1
            group[0] = [1]
        else:
            a0=0
            group[0] = []

        # Should ind have more than one row
        if len(ind2.index) > 1:
            wpm_ind_ind2 = bpmind.wpm['ind'].iloc[(ind2.index)%6]
            # print('====  ind2.shape[1] > 1  ====')
            # print(wpm_ind_ind1)
            # input()

            # TODO: pathsim incomplete
            PATH_sim[1] = pathsim.pathsim(wpm_ind_ind2)
            b0=len(list(set(group[0])))
            b0=0

        # Has exactly one row
        elif len(ind2.index) == 1:
            b0=1
            group[1] = [1]
        else:
            b0=0
            group[1] = []

        g1 = len(list(set(group[0])))
        g2 = len(list(set(group[1])))

        PATH_group.append(group[0] + list(map(lambda x : x + g1, group[1])))
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
