import numpy as np
import pandas as pd
import pickle
from classes import fdrresultsclass as fdrr
from classes import GenstatsOut
from classes import Stats


# fdrsampleperm() computes False Discovery Rates for BPM/WPM/PATH modules
#
# INPUTS:
#   ssmFile: Interaction networks file(path to file) in the pickle format.
#   BPMindFile: file containing SNP ids for BPM/WPMs in pickle format.
#   pcut: p-value cutoff for BPM/WPM/PATH to be considered significant and to be in FDR computing process
#   minPath: minimum size for a pathway to be considered as WPM and in BPM.
#   N: Number of random networks
#
# OUTPUTS:
#   results_<ssmFile without extension>.pkl - This pickle file contains a fdrresultclass class with fields:
#       - bpm_pv: empirical p-values for BPMs
#       - wpm_pv: empirical p-values for WPMs
#       - path_pv: empirical p-values for PATHs
#       - bpm_ranksum: -log10 ranksum p-values for BPMs
#       - wpm_ranksum: -log10 ranksum p-values for WPMs
#       - path_ranksum: -log10 ranksum p-values for PATHs
#       - fdrbpm2: FDR for BPMs
#       - fdrwpm2: FDR for WPMs
#       - fdrpath2: FDR for PATHs
# 





def fdrsampleperm(ssmFile, BPMindFile, pcut, minpath, N, genesetname=None):

    pklin = open(BPMindFile,"rb")
    bpmf = pickle.load(pklin)
    pklin.close()

    tdf = pd.DataFrame()

    ## create dataframes
    bpmdf = pd.DataFrame()
    bpm_pvdf = pd.DataFrame()
    wpmdf = pd.DataFrame()
    wpm_pvdf = pd.DataFrame()
    pathdf = pd.DataFrame()
    path_pvdf = pd.DataFrame()
    bpm_cols, bpm_pv_cols, wpm_cols, wpm_pv_cols, path_cols, path_pv_cols = [], [], [], [], [], []


    for i in range(0, N+1):
        tssmFile = ssmFile.replace("_R0", "_R" + str(i))
        tssm_tmp = tssmFile.split('/')
        tssm_tmp[-1] = 'genstats_' + tssm_tmp[-1]
        genstatsfile = '/'.join(tssm_tmp)

        ## load genstats file
        pklin = open(genstatsfile,"rb")
        gs = pickle.load(pklin)
        pklin.close()

        ## retrieve bpm and add to dataframe
        bpm1 = gs.protective_stats.bpm_local
        bpm2 = gs.risk_stats.bpm_local
        bpm_tmp = np.concatenate((bpm1,bpm2))
        bpm_colname = "bpm"+str(i)
        ## add bpm to dataframe
        bpmdf[bpm_colname] = pd.Series(bpm_tmp)

        ## bpm_pv
        bpm_pv1 = gs.protective_stats.bpm_local_pv
        bpm_pv2 = gs.risk_stats.bpm_local_pv
        bpm_pv_tmp = np.concatenate((bpm_pv1,bpm_pv2))
        bpm_pv_colname = "bpm_pv"+str(i)
        bpm_pvdf[bpm_pv_colname] = pd.Series(bpm_pv_tmp)

        ## wpm
        wpm1 = gs.protective_stats.wpm_local
        wpm2 = gs.risk_stats.wpm_local
        wpm_tmp = np.concatenate((wpm1,wpm2))
        wpm_colname = "wpm"+str(i)
        wpmdf[wpm_colname] = pd.Series(wpm_tmp)

        ## wpm_pv
        wpm_pv1 = gs.protective_stats.wpm_local_pv
        wpm_pv2 = gs.risk_stats.wpm_local_pv
        wpm_pv_tmp = np.concatenate((wpm_pv1,wpm_pv2))
        wpm_pv_colname = "wpm_pv"+str(i)
        wpm_pvdf[wpm_pv_colname] = pd.Series(wpm_pv_tmp)

        ## path
        path1 = gs.protective_stats.path_degree
        path2 = gs.risk_stats.path_degree
        path_tmp = np.concatenate((path1,path2))
        path_colname = "path"+str(i)
        pathdf[path_colname] = pd.Series(path_tmp)

        ## path_pv
        path_pv1 = gs.protective_stats.path_degree_pv
        path_pv2 = gs.risk_stats.path_degree_pv
        path_pv_tmp = np.concatenate((path_pv1,path_pv2))
        path_pv_colname = "path_pv"+str(i)
        path_pvdf[path_pv_colname] = pd.Series(path_pv_tmp)

        ## update column names
        bpm_cols.append(bpm_colname)
        bpm_pv_cols.append(bpm_pv_colname)
        wpm_cols.append(wpm_colname)
        wpm_pv_cols.append(wpm_pv_colname)
        path_cols.append(path_colname)
        path_pv_cols.append(path_pv_colname)



    fdf = pd.concat([bpmf.bpm, bpmdf, bpm_pvdf],axis=1)
    fdf.loc[((fdf['ind1size'] < minpath) | (fdf['ind2size'] < minpath)), bpm_cols] = 0
    fdf.loc[((fdf['ind1size'] < minpath) | (fdf['ind2size'] < minpath)), bpm_pv_cols] = 1
    # print(fdf)

    bpmdf = fdf[bpm_cols]
    bpm_pvdf = fdf[bpm_pv_cols]

    # calling calculate_fdr() function to compute FDRs
    fdrBPM1, fdrBPM2 = calculate_fdr(bpmdf, bpm_cols, bpm_pvdf, bpm_pv_cols, pcut, N, 'bpm')
    fdrWPM1, fdrWPM2 = calculate_fdr(wpmdf, wpm_cols, wpm_pvdf, wpm_pv_cols, 0.05, N, 'wpm')
    fdrPATH1, fdrPATH2 = calculate_fdr(pathdf, path_cols, path_pvdf, path_pv_cols, 0.05, N, 'path')

    bpm_ranksum = bpmdf[bpm_cols[0:1]]
    bpm_ranksum.columns = ['bpm_ranksum']
    wpm_ranksum = wpmdf[wpm_cols[0:1]]
    wpm_ranksum.columns = ['wpm_ranksum']
    path_ranksum = pathdf[path_cols[0:1]]
    path_ranksum.columns = ['path_ranksum']
    bpm_pv = bpm_pvdf[bpm_pv_cols[0:1]]
    bpm_pv.columns = ['bpm_pv']
    wpm_pv = wpm_pvdf[wpm_pv_cols[0:1]]
    wpm_pv.columns = ['wpm_pv']
    path_pv = path_pvdf[path_pv_cols[0:1]]
    path_pv.columns = ['path_pv']


    ssm_tmp = ssmFile.split('/')
    ssm_tmp[-1] = 'results_' + ssm_tmp[-1]
    outfilename = '/'.join(ssm_tmp)
    save_obj = fdrr.fdrrclass(bpm_pv, wpm_pv, path_pv, bpm_ranksum, wpm_ranksum,
        path_ranksum, fdrBPM1, fdrBPM2, fdrWPM1, fdrWPM2, fdrPATH1, fdrPATH2)

    final = open(outfilename, 'wb')
    pickle.dump(save_obj, final)
    final.close()

def calculate_fdr(sdf, sdf_cols, pvdf, pv_cols, pcut, N, type):
    first_pv_col = pv_cols[0:1]
    first_sdf_col = sdf_cols[0:1]
    sdf1 = sdf[first_sdf_col]
    sdf_rest = sdf[sdf_cols[1:]]
    pv1 = pvdf[first_pv_col]
    pv_rest = pvdf[pv_cols[1:]]
    vals = pv1[pv1<=pcut]

    valid_pvs, vrows, vcols, vpv1 = [], [], [], []

    valid_row = vals.first_valid_index()
    while valid_row != None:

        valid_col = vals.loc[valid_row].first_valid_index()
        valid_pv = vals.loc[valid_row][valid_col]

        vrows.append(valid_row)
        vcols.append(valid_col)
        valid_pvs.append(valid_pv)
        vpv1.append(sdf1.loc[valid_row][0])

        vals.loc[valid_row][first_pv_col] = np.nan
        valid_row = vals.first_valid_index()


    m1, m2, n1, n2 = [], [], [], []

    for i in range(len(valid_pvs)):
        pv1_val_cnts = pv1[pv1<=valid_pvs[i]].count().sum()
        pv1rest_val_cnts = pv_rest[pv_rest<=valid_pvs[i]].count().sum()/N
        if (type == 'bpm'):
            # checking BPMs with lower emp. p-vals in real network
            pvf1 = pv1[(pv1<=valid_pvs[i])].ge(0).to_numpy()
            # checking BPMs with lower ranksum p-vals in real network
            svf1 = sdf1[(sdf1>=vpv1[i])].ge(0).to_numpy()
            tfs1 = pvf1&svf1
            # checking BPMs with lower emp. p-vals in random networks
            pvf2 = pv_rest[pv_rest<=valid_pvs[i]].ge(0).to_numpy()
            # checking BPMs with lower ranksum p-vals in random networks
            svf2 = sdf_rest[sdf_rest>=vpv1[i]].ge(0).to_numpy()
            tfs2 = pvf2&svf2
            rs_val_cnts = tfs1.sum()
            rsrest_val_cnts = tfs2.sum()/N
        else:
            # checking WPMs/PATHs with lower emp. p-vals in real network
            pvf1 = pv1[(pv1<=valid_pvs[i])].ge(0).to_numpy()
            # checking WPMs/PATHs with lower ranksum p-vals in real network
            svf1 = sdf1[(sdf1>=vpv1[i])].ge(0).to_numpy()
            tfs1 = pvf1&svf1
            # checking WPMs/PATHs with lower emp. p-vals in random networks
            pvf2 = pv_rest[pv_rest<=valid_pvs[i]].ge(0).to_numpy()
            # checking WPMs/PATHs with lower ranksum p-vals in random networks
            svf2 = sdf_rest[sdf_rest>=vpv1[i]].ge(0).to_numpy()
            tfs2 = pvf2&svf2
            rs_val_cnts = tfs1.sum()
            rsrest_val_cnts = tfs2.sum()/N

        m1.append(pv1_val_cnts)
        m2.append(pv1rest_val_cnts)
        n1.append(rs_val_cnts)
        n2.append(rsrest_val_cnts)


    fdr1 = np.nan_to_num(np.array(m2)/np.array(m1))
    fdr2 = np.nan_to_num(np.array(n2)/np.array(n1))


    rfdr1 = pd.DataFrame(np.ones(pv1.shape), columns=first_pv_col)
    rfdr2 = pd.DataFrame(np.ones(pv1.shape), columns=first_pv_col)



    testM = list(map(list, zip(vrows, valid_pvs, fdr1)))

    testM.sort(key=lambda x: -x[2])

    testM = np.array(testM)
    for i in range(len(testM)):
        idx = testM[:,1] >= testM[i,1]
        testM[i,2] = np.min(testM[idx,2])

    # assign FDR to BPMs
    for i in range(len(testM)):
        rfdr1.loc[testM[i][0]] = testM[i][2]

    if type == 'bpm':
        testM = list(map(list, zip(vrows, valid_pvs, list(fdr2), vpv1)))
    else:
        testM = list(map(list, zip(vrows, valid_pvs, list(fdr2), [round(x) for x in vpv1])))
    testM.sort(key=lambda x: (-x[1]))

    # correct FDRs so BPM/WPM/PATH with lower p-vals does not have larger FDRs
    testM = np.array(testM)
    for i in range(len(testM)):
        idx = (testM[:,1] >= testM[i,1]) & (testM[:,3] <= testM[i,3])
        testM[i,2] = np.min(testM[idx,2])

    # assign FDR to BPMs
    for i in range(len(testM)):
        rfdr2.loc[testM[i][0]] = testM[i][2]

    rfdr1.columns = [(type + str(1))]
    rfdr2.columns = [(type + str(2))]

    return rfdr1, rfdr2
