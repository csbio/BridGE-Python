import pandas as pd
import pickle
from analysis_tools import summarize_bpm as sbpm
# from analysis_tools import summarize_wpm as swpm

def collectresults(resultsfile,fdrcut,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,validationfile):

    pklin = open(resultsfile,"rb")
    rfd = pickle.load(pklin)
    pklin.close()

    fdrBPM, fdrWPM, fdrPATH = rfd.fdrbpm2, rfd.fdrwpm2, rfd.fdrpath2
    bpm_pv, wpm_pv, path_pv = rfd.bpm_pv, rfd.wpm_pv, rfd.path_pv
    bpm_ranksum, wpm_ranksum, path_ranksum = rfd.bpm_ranksum, rfd.wpm_ranksum, rfd.path_ranksum

    ind_bpm = fdrBPM[fdrBPM<=fdrcut].dropna()
    ind_wpm = fdrWPM[fdrWPM<=fdrcut].dropna()
    ind_path = fdrPATH[fdrPATH<=fdrcut].dropna()

    diseasemodel = resultsfile.split('_')[-2]

    print(diseasemodel)
    print(ind_bpm)
    print(ind_wpm)
    print(ind_path)

    if len(ind_bpm) > 0:
        print('bpm_check')
        fdrBPM = fdrBPM.iloc[ind_bpm.index].transpose()
        bpm_pv_discovery = bpm_pv.iloc[ind_bpm.index].transpose()
        bpm_ranksum_discovery = bpm_ranksum.iloc[ind_bpm.index].transpose()
        print(fdrBPM)
        print(bpm_pv_discovery)
        print(bpm_ranksum_discovery)

    if len(ind_wpm) > 0:
        print('wpm_check')
        fdrWPM = fdrWPM.iloc[ind_wpm.index].transpose()
        wpm_pv_discovery = wpm_pv.iloc[ind_wpm.index].transpose()
        wpm_ranksum_discovery = wpm_ranksum.iloc[ind_wpm.index].transpose()
        print(fdrWPM)
        print(wpm_pv_discovery)
        print(wpm_ranksum_discovery)

    if len(ind_path) > 0:
        print('path_check')
        fdrPATH = fdrPATH.iloc[ind_path.index].transpose()
        path_pv_discovery = path_pv.iloc[ind_path.index].transpose()
        path_ranksum_discovery = path_ranksum.iloc[ind_path.index].transpose()
        print(fdrPATH)
        print(path_pv_discovery)
        print(path_ranksum_discovery)


    if (len(ind_bpm) > 0 or len(ind_wpm) > 0 or len(ind_path) > 0):
        pklin = open(bpmindfile,"rb")
        bpm = pickle.load(pklin)
        pklin.close()
        pklin = open(snppathwayfile,"rb")
        snpset = pickle.load(pklin)
        pklin.close()

        genesetfile = "../refdata/" + snpset.geneset.split('/')[-1]
        print(genesetfile)

        if len(ind_bpm) > 0:
            path1 = bpm.bpm['path1names']
            path2 = bpm.bpm['path2names']
            path1 = path1.append(path1).reset_index().drop(columns = ['index'])
            path2 = path2.append(path2).reset_index().drop(columns = ['index'])

            path1 = path1.iloc[ind_bpm.index]
            path2 = path2.iloc[ind_bpm.index]
            print(path1)
            print(path2)

            bpm_size = bpm.bpm['size']
            bpm_size = bpm_size.append(bpm_size).reset_index().drop(columns = ['index'])
            bpm_size = bpm_size.iloc[ind_bpm.index]
            print(bpm_size)

            eff_bpm = []
            for i in range(len(ind_bpm)):
                if ind_bpm.iloc[i].name<=len(bpm.bpm.index):
                    eff_bpm.append('protective')
                else:
                    eff_bpm.append('risk')
            print(eff_bpm)

    return
