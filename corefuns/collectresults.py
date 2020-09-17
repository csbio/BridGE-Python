import pandas as pd
import pickle
from corefuns import check_BPM_WPM_redundancy as cbwr
from classes import fdrresultsclass
from classes import bpmindclass
from corefuns import get_interaction_pair as gpair

def collectresults(resultsfile,fdrcut,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,validationfile=None):

    pklin = open(resultsfile,"rb")
    rfd = pickle.load(pklin)
    pklin.close()

    fdrBPM, fdrWPM, fdrPATH = rfd.fdrbpm2, rfd.fdrwpm2, rfd.fdrpath2
    bpm_pv, wpm_pv, path_pv = rfd.bpm_pv, rfd.wpm_pv, rfd.path_pv
    bpm_ranksum, wpm_ranksum, path_ranksum = rfd.bpm_ranksum, rfd.wpm_ranksum, rfd.path_ranksum


    ind_bpm = fdrBPM[fdrBPM<=fdrcut].dropna()
    ind_wpm = fdrWPM[fdrWPM<=fdrcut].dropna()
    ind_path = fdrPATH[fdrPATH<=fdrcut].dropna()

    fdrBPM_in = fdrBPM
    fdrWPM_in = fdrWPM
    fdrPATH_in = fdrPATH
    diseasemodel = resultsfile.split('_')[-2]

    # print(diseasemodel)
    # print(ind_bpm)
    # print(ind_wpm)
    # print(ind_path)

    if not ind_bpm.empty:
        # print('bpm_check')
        fdrBPM = fdrBPM.iloc[ind_bpm.index]
        bpm_pv_discovery = bpm_pv.iloc[ind_bpm.index]
        bpm_ranksum_discovery = bpm_ranksum.iloc[ind_bpm.index]
        # print(fdrBPM)
        # print(bpm_pv_discovery)
        # print(bpm_ranksum_discovery)

    if not ind_wpm.empty:
        # print('wpm_check')
        fdrWPM = fdrWPM.iloc[ind_wpm.index]
        wpm_pv_discovery = wpm_pv.iloc[ind_wpm.index]
        wpm_ranksum_discovery = wpm_ranksum.iloc[ind_wpm.index]
        # print(fdrWPM)
        # print(wpm_pv_discovery)
        # print(wpm_ranksum_discovery)

    if  not ind_path.empty:
        # print('path_check')
        fdrPATH = fdrPATH.iloc[ind_path.index]
        path_pv_discovery = path_pv.iloc[ind_path.index]
        path_ranksum_discovery = path_ranksum.iloc[ind_path.index]
        # print(fdrPATH)
        # print(path_pv_discovery)
        # print(path_ranksum_discovery)

    if not ind_bpm.empty or not ind_wpm.empty or not ind_path.empty:
        pklin = open(bpmindfile,"rb")
        bpm = pickle.load(pklin)
        pklin.close()
        pklin = open(snppathwayfile,"rb")
        snpset = pickle.load(pklin)
        pklin.close()

        genesetfile = snpset.geneset


        # Remove when refactoring for get_interaction_pair()
        header = ['snps', 'genes', 'snp_mean_gi', 'snp_mean_gi_bg', 'gi_fold', 'gi_hyge']

        if not ind_bpm.empty:
            path1 = bpm.bpm['path1names']
            path2 = bpm.bpm['path2names']
            path1 = path1.append(path1).reset_index().drop(columns = ['index'])
            path2 = path2.append(path2).reset_index().drop(columns = ['index'])

            path1 = path1.iloc[ind_bpm.index]
            path2 = path2.iloc[ind_bpm.index]
            # print(path1)
            # print(path2)

            bpm_size = bpm.bpm['size']
            bpm_size = bpm_size.append(bpm_size).reset_index().drop(columns = ['index'])
            bpm_size = bpm_size.iloc[ind_bpm.index]
            # print(bpm_size)

            eff_bpm = []
            for i in range(len(ind_bpm)):
                if ind_bpm.iloc[i].name<=len(bpm.bpm.index):
                    eff_bpm.append('protective')
                else:
                    eff_bpm.append('risk')

            eff_bpm = pd.DataFrame(eff_bpm, columns=['eff_bpm'])
            # print(eff_bpm)

            bpm_path1_drivers, bpm_path2_drivers = [], []

            for i in range(len(ind_bpm)):
                [output_path1_snp, output_path2_snp] = gpair.get_interaction_pair(path1.iloc[i][0],path2.iloc[i][0],eff_bpm.iloc[i][0],ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile)

                idx = output_path1_snp[output_path1_snp['gi_fold']>1].dropna()

                bpm_snp_tmp = []
                if not idx.empty:
                    for row in idx.head(20).iterrows():
                        temp_str = row[1]['snps']+'_'+row[1]['genes']+'_'+str(round(row[1]['gi_fold'], 2))
                        bpm_snp_tmp.append(temp_str)
                    bpm_path1_drivers.append(';'.join(bpm_snp_tmp))

                else:
                    bpm_path1_drivers.append('')

                idx = output_path2_snp[output_path2_snp['gi_fold']>1].dropna()

                bpm_snp_tmp = []

                if not idx.empty:
                    for row in idx.head(20).iterrows():
                        temp_str = row[1]['snps']+'_'+row[1]['genes']+'_'+str(round(row[1]['gi_fold'], 2))
                        bpm_snp_tmp.append(temp_str)
                        # print(temp_str)
                    bpm_path2_drivers.append(';'.join(bpm_snp_tmp))

                else:
                    bpm_path2_drivers.append('')


            bpm_path1_drivers = pd.DataFrame(bpm_path1_drivers, columns=['bpm_path1_drivers'])
            bpm_path2_drivers = pd.DataFrame(bpm_path2_drivers, columns=['bpm_path2_drivers'])

        if not ind_wpm.empty:
            pathway = bpm.wpm['pathway']
            pathway = pathway.append(pathway).reset_index().drop(columns = ['index'])
            path_wpm = pathway.iloc[ind_wpm.index]

            wpm_size = bpm.wpm['size']
            wpm_size = wpm_size.append(wpm_size).reset_index().drop(columns = ['index'])
            wpm_size = wpm_size.iloc[ind_wpm.index]

            eff_wpm = []
            for i in range(len(ind_wpm)):
                if ind_wpm.iloc[i].name<=len(bpm.wpm.index):
                    eff_wpm.append('protective')
                else:
                    eff_wpm.append('risk')
            eff_wpm = pd.DataFrame(eff_wpm, columns=['eff_wpm'])
            # print(eff_wpm)

            for i in range(len(ind_wpm)):
                [output_path_snp,x] = get_interaction_pair(path_wpm.iloc[i][0],path_wpm.iloc[i][0],eff_bpm.iloc[i][0],ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile)

                idx = output_path_snp[output_path_snp['gi_fold']>1].dropna()

                wpm_snp_tmp, wpm_path_drivers = [], []
                if idx.shape[1] > 0:
                    for row in idx.head(20).iterrows():
                        temp_str = row[1]['snps']+'_'+row[1]['genes']+'_'+str(round(row[1]['gi_fold'], 2))
                        wpm_snp_tmp.append(temp_str)
                    wpm_path_drivers.append(';'.join(bpm_snp_tmp))

                else:
                    wpm_path_drivers.append('')

            wpm_path_drivers = pd.DataFrame(wpm_path_drivers, columns=['wpm_path_drivers'])

        if not ind_path.empty:
            pathway = bpm.wpm['pathway']
            pathway = pathway.append(pathway).reset_index().drop(columns = ['index'])
            path_path = pathway.iloc[ind_path.index]

            path_size = bpm.wpm['indsize']
            path_size = path_size.append(path_size).reset_index().drop(columns = ['index'])
            path_size = path_size.iloc[ind_path.index]

            eff_path = []
            for i in range(len(ind_path)):
                if ind_path.iloc[i].name<=len(bpm.wpm.index):
                    eff_path.append('protective')
                else:
                    eff_path.append('risk')
            eff_path = pd.DataFrame(eff_path, columns=['eff_path'])
            # print(eff_path)

    if not ind_bpm.empty:
        # fdrBPM.rename(columns={"bpm2": "fdrBPM"}, inplace = True)
        # bpm_size.rename(columns={"bpm2": "bpm_size"}, inplace = True)

        path1 = path1.reset_index(drop = True)
        path2 = path2.reset_index(drop = True)
        fdrBPM = fdrBPM.reset_index(drop = True)
        bpm_size = bpm_size.reset_index(drop = True)
        bpm_pv_discovery = bpm_pv_discovery.reset_index(drop = True)
        bpm_ranksum_discovery = bpm_ranksum_discovery.reset_index(drop = True)


        fdrBPM = round(fdrBPM, 2)
        bpm_ranksum_discovery = round(bpm_ranksum_discovery, 2)
        results = [path1,path2,fdrBPM,eff_bpm,bpm_size,bpm_pv_discovery,bpm_ranksum_discovery,bpm_path1_drivers,bpm_path2_drivers]
        output_bpm_table = pd.concat(results, axis=1, sort=False)
        # print('output_bpm_table')
        # print(output_bpm_table)
        # print(output_bpm_table['size'])
        # input()

    if not ind_wpm.empty:
        # fdrWPM.rename(columns={"wpm2": "fdrWPM"}, inplace = True)
        # wpm_size.rename(columns={"wpm2": "wpm_size"}, inplace = True)

        path = path_wpm.reset_index(drop = True)
        fdrWPM = fdrWPM.reset_index(drop = True)
        wpm_size = wpm_size.reset_index(drop = True)
        wpm_pv_discovery = wpm_pv_discovery.reset_index(drop = True)
        wpm_ranksum_discovery = wpm_ranksum_discovery.reset_index(drop = True)


        fdrWPM = round(fdrWPM, 2)
        wpm_ranksum_discovery = round(wpm_ranksum_discovery, 2)
        results = [path,fdrWPM,eff_wpm,wpm_size,wpm_pv_discovery,wpm_ranksum_discovery,wpm_path_drivers]
        output_wpm_table = pd.concat(results, axis=1, sort=False)
        # print('output_wpm_table')
        # print(output_wpm_table)
        # input()

    if not ind_path.empty:
        # fdrPATH.rename(columns={"path2": "fdrPATH"}, inplace = True)
        # path_size.rename(columns={"path2": "path_size"}, inplace = True)

        path = path_path.reset_index(drop = True)
        fdrPATH = fdrPATH.reset_index(drop = True)
        path_size = path_size.reset_index(drop = True)
        path_pv_discovery = path_pv_discovery.reset_index(drop = True)
        path_ranksum_discovery = path_ranksum_discovery.reset_index(drop = True)


        fdrPATH = round(fdrPATH, 2)
        path_ranksum_discovery = round(path_ranksum_discovery, 2)
        results = [path,fdrPATH,eff_path,path_size,path_pv_discovery,path_ranksum_discovery]
        output_path_table = pd.concat(results, axis=1, sort=False)
        # print('output_path_table')
        # print(output_path_table)
        # input()


    output_discovery_summary, output_noRD_discovery_summary = [], []

    output_discovery_summary.append([fdrBPM['bpm2'].min(), (fdrBPM['bpm2']<=0.05).sum(), (fdrBPM['bpm2']<=0.1).sum(), (fdrBPM['bpm2']<=0.15).sum(), (fdrBPM['bpm2']<=0.2).sum(), (fdrBPM['bpm2']<=0.25).sum(), (fdrBPM['bpm2']<=0.3).sum(), (fdrBPM['bpm2']<=0.35).sum(), (fdrBPM['bpm2']<=0.4).sum()])

    if not fdrWPM.empty:
        output_discovery_summary.append([fdrWPM['wpm2'].min(), (fdrWPM['wpm2']<=0.05).sum(), (fdrWPM['wpm2']<=0.1).sum(), (fdrWPM['wpm2']<=0.15).sum(), (fdrWPM['wpm2']<=0.2).sum(), (fdrWPM['wpm2']<=0.25).sum(), (fdrWPM['wpm2']<=0.3).sum(), (fdrWPM['wpm2']<=0.35).sum(), (fdrWPM['wpm2']<=0.4).sum()])
        output_discovery_summary.append([fdrPATH['path2'].min(), (fdrPATH['path2']<=0.05).sum(), (fdrPATH['path2']<=0.1).sum(), (fdrPATH['path2']<=0.15).sum(), (fdrPATH['path2']<=0.2).sum(), (fdrPATH['path2']<=0.25).sum(), (fdrPATH['path2']<=0.3).sum(), (fdrPATH['path2']<=0.35).sum(), (fdrPATH['path2']<=0.4).sum()])


    output_header = ['minfdr','fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40']


    BPM_nosig_noRD,WPM_nosig_noRD,PATH_nosig_noRD,BPM_group_tmp,WPM_group_tmp,PATH_group_tmp = cbwr.check_BPM_WPM_redundancy(fdrBPM_in, fdrWPM_in, fdrPATH_in, bpmindfile, 0.4)

    output_noRD_discovery_summary.append([fdrBPM['bpm2'].min()] + BPM_nosig_noRD)

    if not fdrWPM.empty:
        output_noRD_discovery_summary.append([fdrWPM['wpm2'].min()] + WPM_nosig_noRD)
        output_noRD_discovery_summary.append([fdrPATH['path2'].min()] + PATH_nosig_noRD)

    if not fdrWPM.empty:
        output_discovery_index = ['BPM','WPM','PATH']
    else:
        output_discovery_index = ['BPM']

    output_discovery_summary = pd.DataFrame(output_discovery_summary, columns=output_header, index=output_discovery_index)
    output_noRD_discovery_summary = pd.DataFrame(output_noRD_discovery_summary, columns=output_header, index=output_discovery_index)


    final_filename = resultsfile.split('/').pop()
    final_filename = "data/output_" +  final_filename.rsplit('.',1)[0] + ".xls"


    #### DOUBLE CHECK SORTING WHEN FIXING BPMSIM AND PATHSIM. ####
    if not ind_bpm.empty:
        output_bpm_table.rename(columns={"bpm2": "fdrBPM"}, inplace = True)
        output_bpm_table.rename(columns={"indsize": "bpm_size"}, inplace = True)
        output_bpm_table.sort_values(by=['fdrBPM', 'bpm_pv', 'bpm_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_bpm_table = output_bpm_table.reset_index(drop = True)
        # print(output_bpm_table)
        # input()

    if not ind_wpm.empty:
        output_wpm_table.rename(columns={"wpm2": "fdrWPM"}, inplace = True)
        output_wpm_table.rename(columns={"indsize": "wpm_size"}, inplace = True)
        output_wpm_table.sort_values(by=['wpm2', 'wpm_pv', 'wpm_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_wpm_table = output_wpm_table.reset_index(drop = True)
        # print(output_wpm_table)
        # input()

    if  not ind_path.empty:
        output_path_table.rename(columns={"path2": "fdrPATH"}, inplace = True)
        output_path_table.rename(columns={"indsize": "path_size"}, inplace = True)
        output_path_table.sort_values(by=['fdrPATH', 'path_pv', 'path_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_path_table = output_path_table.reset_index(drop = True)
        # print(output_path_table)
        # input()

    writer = pd.ExcelWriter(final_filename)

    output_discovery_summary.to_excel(writer, sheet_name='output_discovery_summary')
    output_noRD_discovery_summary.to_excel(writer, sheet_name='output_noRD_discovery_summary')

    if not ind_bpm.empty:
        output_bpm_table.to_excel(writer, sheet_name='output_bpm_table')
    if not ind_wpm.empty:
        output_wpm_table.to_excel(writer, sheet_name='output_wpm_table')
    if not ind_path.empty:
        output_path_table.to_excel(writer, sheet_name='output_path_table')

    writer.save()

    return


