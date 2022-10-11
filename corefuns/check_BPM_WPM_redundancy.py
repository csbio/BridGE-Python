import pandas as pd
import numpy as np
import pickle
from corefuns import bpmsim as bpmsim
from corefuns import pathsim as pathsim
from scipy import sparse




# check_BPM_WPM_redundancy() finds redundant BPM/WPM/PATHs and group the similar ones into redundant groups for different FDR thresholds.
#
# INPUTS:
#   fdrBPM: A dataframe contating the information(FDR,pathway names, etc) about BPMs that passed general FDR threshold 
#   fdrWPM: A dataframe contating the information(FDR,pathway names, etc) about WPMs that passed general FDR threshold 
#   fdrPATH: A dataframe contating the information(FDR,pathway names, etc) about PATHs that passed general FDR threshold 
#   bpmindfile: file containing SNP ids for BPM/WPMs in pickle format.
#   FDRcut: Maximum FDR threshold.
#
# OUTPUTS:
#   returns an array containing 6 lists
#       - BPM_nosig_noRD: List containing number of non-redundant BPMs for each FDR threshold(incremented by 0.05)
#       - WPM_nosig_noRD: List containing number of non-redundant WPMs for each FDR threshold(incremented by 0.05)
#       - PATH_nosig_noRD: List containing number of non-redundant PATHs for each FDR threshold(incremented by 0.05)
#       - BPM_group: List of multiple lists each containing redundant group membership indices for BPMs in the same order as fdrBPM.
#       - WPM_group: List of multiple lists each containing redundant group membership indices for WPMs in the same order as fdrwPM.
#       - PATH_group: List of multiple lists each containing redundant group membership indices for PATHs in the same order as fdrPATH.
# 


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
        ind2 = ind[nnz:] - len(bpmind.bpm['size'])

        ind = np.array(fdrWPM[fdrWPM<=fdrcut].dropna().index)
        nnz = (ind<=len(bpmind.wpm['size'])).sum()

        ind3 = ind[0:nnz]
        ind4 = ind[nnz:] - len(bpmind.wpm['size'])


        group, BPM_sim, noRD = [0]*4, [0]*4, [0]*4

        # Should ind have more than one row
        if len(ind1) > 1:
            bpm_ind1_ind1 = bpmind.bpm['ind1'][ind1]
            bpm_ind2_ind1 = bpmind.bpm['ind2'][ind1]

            BPM_sim[0] = bpmsim.bpmsim(bpm_ind1_ind1, bpm_ind2_ind1, bpm_ind1_ind1, bpm_ind2_ind1)

            TTT = (BPM_sim[0]>=0.25).astype(int)
            #noRD[0], group[0] = sparse.csgraph.connected_components(TTT)
            ## retrieve indices and sort based on FDR
            tmp_ind = ind1
            fdrs = fdrBPM.loc[tmp_ind]['bpm2']
            tmp_idx = fdrs.argsort(kind = 'stable').to_numpy()
            fdrs_sorted = fdrs.iloc[tmp_idx]
            ## search for group and assign
            tmp_group = np.zeros((fdrs_sorted.shape[0],))
            for x1 in range(fdrs_sorted.shape[0]):
                if x1 == 0:
                    continue
                ttt_idx1 = tmp_idx[x1]
                for x2 in range(x1+1):
                    if x2 == x1:
                        tmp_group[x1] = max(tmp_group) + 1
                    else:
                        ttt_idx2 = tmp_idx[x2]
                        if TTT[ttt_idx1,ttt_idx2] == 1:
                            tmp_group[x1] = tmp_group[x2]
                            break
            group[0] = tmp_group
            a0=len(list(set(group[0])))

        # Has exactly one row
        elif len(ind1) == 1:
            a0=1
            group[0] = [0]
        else:
            a0=0
            group[0] = []

        # Should ind have more than one row
        if len(ind2) > 1:
            bpm_ind1_ind2 = bpmind.bpm['ind1'][ind2]
            bpm_ind2_ind2 = bpmind.bpm['ind2'][ind2]

            BPM_sim[1] = bpmsim.bpmsim(bpm_ind1_ind2, bpm_ind2_ind2, bpm_ind1_ind2, bpm_ind2_ind2)

            TTT = (BPM_sim[1]>=0.25).astype(int)
            #noRD[1], group[1] = sparse.csgraph.connected_components(TTT)
            #b0=len(list(set(group[1])))
            ## retrieve indices and sort based on FDR
            tmp_ind = ind2 +  len(bpmind.bpm['size'])
            fdrs = fdrBPM.loc[tmp_ind]['bpm2']
            tmp_idx = fdrs.argsort(kind='stable').to_numpy()
            fdrs_sorted = fdrs.iloc[tmp_idx]
            ## search for group and assign
            tmp_group = np.zeros((fdrs_sorted.shape[0],))
            for x1 in range(fdrs_sorted.shape[0]):
                if x1 == 0:
                    continue
                ttt_idx1 = tmp_idx[x1]
                for x2 in range(x1+1):
                    if x2 == x1:
                        tmp_group[x1] = max(tmp_group) + 1
                    else:
                        ttt_idx2 = tmp_idx[x2]
                        if TTT[ttt_idx1,ttt_idx2] == 1:
                            tmp_group[x1] = tmp_group[x2]
                            break
            group[1] = tmp_group
            b0=len(list(set(group[1])))
            


        # Has exactly one row
        elif len(ind2) == 1:
            b0=1
            group[1] = [0]
        else:
            b0=0
            group[1] = []

        # Should ind have more than one row
        if len(ind3) > 1:
            wpm_ind_ind3 = bpmind.wpm['ind'][ind3]

            BPM_sim[2] = bpmsim.bpmsim(wpm_ind_ind3, wpm_ind_ind3, wpm_ind_ind3, wpm_ind_ind3)

            TTT = (BPM_sim[2]>=0.25).astype(int)
            #noRD[2], group[2] = sparse.csgraph.connected_components(TTT)
            ## retrieve indices and sort based on FDR
            tmp_ind = ind3
            fdrs = fdrWPM.loc[tmp_ind]['wpm2']
            tmp_idx = fdrs.argsort(kind = 'stable').to_numpy()
            fdrs_sorted = fdrs.iloc[tmp_idx]
            ## search for group and assign
            tmp_group = np.zeros((fdrs_sorted.shape[0],))
            for x1 in range(fdrs_sorted.shape[0]):
                if x1 == 0:
                    continue
                ttt_idx1 = tmp_idx[x1]
                for x2 in range(x1+1):
                    if x2 == x1:
                        tmp_group[x1] = max(tmp_group) + 1
                    else:
                        ttt_idx2 = tmp_idx[x2]
                        if TTT[ttt_idx1,ttt_idx2] == 1:
                            tmp_group[x1] = tmp_group[x2]
                            break
            group[2] = tmp_group
            c0=len(list(set(group[2])))


        # Has exactly one row
        elif len(ind3) == 1:
            c0=1
            group[2] = [0]
        else:
            c0=0
            group[2] = []

        # Should ind have more than one row
        if len(ind4) > 1:
            wpm_ind_ind4 = bpmind.wpm['ind'][ind4]

            BPM_sim[3] = bpmsim.bpmsim(wpm_ind_ind4, wpm_ind_ind4, wpm_ind_ind4, wpm_ind_ind4)

            TTT = (BPM_sim[3]>=0.25).astype(int)
            #noRD[3], group[3] = sparse.csgraph.connected_components(TTT)
            tmp_ind = ind4 + len(bpmind.wpm['size'])
            fdrs = fdrWPM.loc[tmp_ind]['wpm2']
            tmp_idx = fdrs.argsort(kind = 'stable').to_numpy()
            fdrs_sorted = fdrs.iloc[tmp_idx]
            ## search for group and assign
            tmp_group = np.zeros((fdrs_sorted.shape[0],))
            for x1 in range(fdrs_sorted.shape[0]):
                if x1 == 0:
                    continue
                ttt_idx1 = tmp_idx[x1]
                for x2 in range(x1+1):
                    if x2 == x1:
                        tmp_group[x1] = max(tmp_group) + 1
                    else:
                        ttt_idx2 = tmp_idx[x2]
                        if TTT[ttt_idx1,ttt_idx2] == 1:
                            tmp_group[x1] = tmp_group[x2]
                            break
            group[3] = tmp_group
            d0=len(list(set(group[3])))


        # Has exactly one row
        elif len(ind4) == 1:
            d0=1
            group[3] = [0]
        else:
            d0=0
            group[3] = []


        g1 = len(list(set(group[0])))
        g2 = len(list(set(group[1])))
        g3 = len(list(set(group[2])))
        g4 = len(list(set(group[3])))


        BPM_group.append(list(group[0]) + list(map(lambda x : x + g1, list(group[1]))))
        BPM_nosig_noRD.append(a0 + b0)

        WPM_group.append(list(group[2]) + list(map(lambda x : x + g3, list(group[3]))))
        WPM_nosig_noRD.append(c0 + d0)

        ind = np.array(fdrPATH[fdrPATH<=fdrcut].dropna().index)
        nnz = (ind<=len(bpmind.wpm['size'])).sum()

        ind1 = ind[0:nnz]
        ind2 = ind[nnz:] - len(bpmind.wpm['size'])

        group, PATH_sim, noRD = [0]*2, [0]*2, [0]*2

        # Should ind have more than one row
        if len(ind1) > 1:
            wpm_ind_ind1 = bpmind.wpm['ind'][ind1]

            PATH_sim[0] = pathsim.pathsim(wpm_ind_ind1)
            TTT = (PATH_sim[0]>=0.25).astype(int)
            #noRD[0], group[0] = sparse.csgraph.connected_components(TTT)
            tmp_ind = ind1
            fdrs = fdrWPM.loc[tmp_ind]['wpm2']
            tmp_idx = fdrs.argsort(kind = 'stable').to_numpy()
            fdrs_sorted = fdrs.iloc[tmp_idx]
            ## search for group and assign
            tmp_group = np.zeros((fdrs_sorted.shape[0],))
            for x1 in range(fdrs_sorted.shape[0]):
                if x1 == 0:
                    continue
                ttt_idx1 = tmp_idx[x1]
                for x2 in range(x1+1):
                    if x2 == x1:
                        tmp_group[x1] = max(tmp_group) + 1
                    else:
                        ttt_idx2 = tmp_idx[x2]
                        if TTT[ttt_idx1,ttt_idx2] == 1:
                            tmp_group[x1] = tmp_group[x2]
                            break
            group[0] = tmp_group
            a0=len(list(set(group[0])))

        # Has exactly one row
        elif len(ind1) == 1:
            a0=1
            group[0] = [0]
        else:
            a0=0
            group[0] = []

        # Should ind have more than one row
        if len(ind2) > 1:
            wpm_ind_ind2 = bpmind.wpm['ind'][ind2]

            PATH_sim[1] = pathsim.pathsim(wpm_ind_ind2)
            TTT = (PATH_sim[1]>=0.25).astype(int)
            #noRD[1], group[1] = sparse.csgraph.connected_components(TTT)
            tmp_ind = ind2 + len(bpmind.wpm['size'])
            fdrs = fdrWPM.loc[tmp_ind]['wpm2']
            tmp_idx = fdrs.argsort(kind = 'stable').to_numpy()
            fdrs_sorted = fdrs.iloc[tmp_idx]
            ## search for group and assign
            tmp_group = np.zeros((fdrs_sorted.shape[0],))
            for x1 in range(fdrs_sorted.shape[0]):
                if x1 == 0:
                    continue
                ttt_idx1 = tmp_idx[x1]
                for x2 in range(x1+1):
                    if x2 == x1:
                        tmp_group[x1] = max(tmp_group) + 1
                    else:
                        ttt_idx2 = tmp_idx[x2]
                        if TTT[ttt_idx1,ttt_idx2] == 1:
                            tmp_group[x1] = tmp_group[x2]
                            break
            group[1] = tmp_group
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


        PATH_group.append(list(group[0]) + list(map(lambda x : x + g1, list(group[1]))))
        PATH_nosig_noRD.append(a0 + b0)


    return BPM_nosig_noRD,WPM_nosig_noRD,PATH_nosig_noRD,BPM_group,WPM_group,PATH_group
