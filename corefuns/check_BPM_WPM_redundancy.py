import pandas as pd
import numpy as np
import pickle
from corefuns import bpmsim as bpmsim
from corefuns import pathsim as pathsim
from scipy import sparse
import sys
import datetime

def time_log():
    print(datetime.datetime.now())
    sys.stdout.flush()

def check_BPM_WPM_redundancy(fdrBPM,fdrWPM,fdrPATH,bpmindfile,FDRcut):
    pklin = open(bpmindfile,"rb")
    bpmind = pickle.load(pklin)
    pklin.close()

    BPM_group, BPM_nosig_noRD = [], []
    WPM_group, WPM_nosig_noRD = [], []
    PATH_group, PATH_nosig_noRD = [], []

    start = 0.05
    increment = 0.05
   
    print('in check redundancy...')

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
            print('returned from bpmsim in ind1')
            time_log()

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
            bpm_ind1_ind2 = bpmind.bpm['ind1'][ind2]
            bpm_ind2_ind2 = bpmind.bpm['ind2'][ind2]

            BPM_sim[1] = bpmsim.bpmsim(bpm_ind1_ind2, bpm_ind2_ind2, bpm_ind1_ind2, bpm_ind2_ind2)

            TTT = (BPM_sim[1]>=0.25).astype(int)
            noRD[1], group[1] = sparse.csgraph.connected_components(TTT)
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
            wpm_ind_ind3 = bpmind.wpm['ind'][ind3]

            BPM_sim[2] = bpmsim.bpmsim(wpm_ind_ind3, wpm_ind_ind3, wpm_ind_ind3, wpm_ind_ind3)

            TTT = (BPM_sim[2]>=0.25).astype(int)
            noRD[2], group[2] = sparse.csgraph.connected_components(TTT)
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
            wpm_ind_ind4 = bpmind.wpm['ind'][ind4]

            BPM_sim[3] = bpmsim.bpmsim(wpm_ind_ind4, wpm_ind_ind4, wpm_ind_ind4, wpm_ind_ind4)

            TTT = (BPM_sim[3]>=0.25).astype(int)
            noRD[3], group[3] = sparse.csgraph.connected_components(TTT)
            d0=len(list(set(group[3])))


        # Has exactly one row
        elif len(ind4) == 1:
            d0=1
            group[3] = [1]
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

        print('processing paths')
        # Should ind have more than one row
        if len(ind1) > 1:
            wpm_ind_ind1 = bpmind.wpm['ind'][ind1]

            PATH_sim[0] = pathsim.pathsim(wpm_ind_ind1)
            print('returned from pathsim in ind1')
            time_log()
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
            wpm_ind_ind2 = bpmind.wpm['ind'][ind2]

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


        PATH_group.append(list(group[0]) + list(map(lambda x : x + g1, list(group[1]))))
        PATH_nosig_noRD.append(a0 + b0)


    return BPM_nosig_noRD,WPM_nosig_noRD,PATH_nosig_noRD,BPM_group,WPM_group,PATH_group
