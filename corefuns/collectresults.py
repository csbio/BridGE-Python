import pandas as pd
import numpy as np
import pickle
from corefuns import check_BPM_WPM_redundancy as cbwr
from classes import fdrresultsclass
from classes import bpmindclass
from corefuns import get_interaction_pair as gpair
from corefuns import pathway_map as pmap


# collectresults() is responsible for collecting and exporting the results to an Excel file. It also calls functions for finding driver genes and non-redundunt modules.
#
# INPUTS:
#   resultsfile: Pickle file containing the FDR, empirical p-values, and ranksum scores of BPM/WPM/PATH modules.
#   fdrcut: FDR threshold for keeping BPM/WPM/PATH modules. 
#   ssmfile: Interaction network file(path to file) for real network in the pickle format.
#   bpmindfile: file containing SNP ids for BPM/WPMs in pickle format.
#   snppathwayfile: Pickle file containing mapping of SNPs to pathways.
#   snpgenemappingfile: Pickle file containing mapping of SNPs to Genes.
#
# OUTPUTS:
#   output_results_<ssmFile without extension>.xls - This Excel file has the following pages:
#       - output_discovery_summary: summary of findings(numbers of) for BPM/WPM/PATH
#       - output_noRD_discovery_summary: summary of non-redundant findings(numbers of) for BPM/WPM/PATH 
#       - output_bpm_table: Table of the all BPMs with FDR below fdrcut with all the stats and driver genes 
#       - BPM_redundant_groups: Table connecting redundunt BPM modules for different FDR thresholds
#       - output_wpm_table: Table of the all WPMs with FDR below fdrcut with all the stats and driver genes 
#       - WPM_redundant_groups: Table connecting redundunt WPM modules for different FDR thresholds
#       - output_path_table: Table of the all PATHs with FDR below fdrcut with all the stats and driver genes 
#       - PATH_redundant_groups: Table connecting redundunt PATH modules for different FDR thresholds
# 



def collectresults(resultsfile,fdrcut,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile):

    pklin = open(resultsfile,"rb")
    rfd = pickle.load(pklin)
    pklin.close()
    # Find project dir
    dir_sp = resultsfile.split('/')
    s = '/'
    project_dir = s.join(dir_sp[:-1])

    fdrBPM, fdrWPM, fdrPATH = rfd.fdrbpm2, rfd.fdrwpm2, rfd.fdrpath2
    bpm_pv, wpm_pv, path_pv = rfd.bpm_pv, rfd.wpm_pv, rfd.path_pv
    bpm_ranksum, wpm_ranksum, path_ranksum = rfd.bpm_ranksum, rfd.wpm_ranksum, rfd.path_ranksum


    ind_bpm = fdrBPM[fdrBPM<=fdrcut].dropna()
    ind_wpm = fdrWPM[fdrWPM<=fdrcut].dropna()
    ind_path = fdrPATH[fdrPATH<=fdrcut].dropna()

    diseasemodel = resultsfile.split('_')[-2]


    if not ind_bpm.empty:
        fdrBPM = fdrBPM.iloc[ind_bpm.index]
        bpm_pv_discovery = bpm_pv.iloc[ind_bpm.index]
        bpm_ranksum_discovery = bpm_ranksum.iloc[ind_bpm.index]

    if not ind_wpm.empty:
        fdrWPM = fdrWPM.iloc[ind_wpm.index]
        wpm_pv_discovery = wpm_pv.iloc[ind_wpm.index]
        wpm_ranksum_discovery = wpm_ranksum.iloc[ind_wpm.index]

    if  not ind_path.empty:
        fdrPATH = fdrPATH.iloc[ind_path.index]
        path_pv_discovery = path_pv.iloc[ind_path.index]
        path_ranksum_discovery = path_ranksum.iloc[ind_path.index]

    if not ind_bpm.empty or not ind_wpm.empty or not ind_path.empty:
        pklin = open(bpmindfile,"rb")
        bpm = pickle.load(pklin)
        pklin.close()
        pklin = open(snppathwayfile,"rb")
        snpset = pickle.load(pklin)
        pklin.close()

        genesetfile = snpset.geneset
        # load pathway ids into a map
        path_ids = {}
        pathway_size = bpm.wpm['pathway'].shape[0]
        for i in range(pathway_size):
            p = bpm.wpm['pathway'][i]
            path_ids[p] = i

        if not ind_bpm.empty:
            path1 = bpm.bpm['path1names']
            path2 = bpm.bpm['path2names']
            path1 = path1.append(path1).reset_index().drop(columns = ['index'])
            path2 = path2.append(path2).reset_index().drop(columns = ['index'])

            path1 = path1.iloc[ind_bpm.index]
            path2 = path2.iloc[ind_bpm.index]

            bpm_size = bpm.bpm['size']
            bpm_size = bpm_size.append(bpm_size).reset_index().drop(columns = ['index'])
            bpm_size = bpm_size.iloc[ind_bpm.index]

            eff_bpm = []
            for i in range(len(ind_bpm)):
                if ind_bpm.iloc[i].name<=len(bpm.bpm.index):
                    eff_bpm.append('protective')
                else:
                    eff_bpm.append('risk')

            eff_bpm = pd.DataFrame(eff_bpm, columns=['eff_bpm'])

            [bpm_path1_drivers, bpm_path2_drivers,tmp] = gpair.get_interaction_pair(len(ind_bpm),path1,path2,eff_bpm,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,path_ids)            

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

            [tmp1,tmp2,wpm_path_drivers] = gpair.get_interaction_pair(len(ind_wpm),path_wpm,path_wpm,eff_wpm,ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile,path_ids)

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

    BPM_nosig_noRD,WPM_nosig_noRD,PATH_nosig_noRD,BPM_group_tmp,WPM_group_tmp,PATH_group_tmp = cbwr.check_BPM_WPM_redundancy(fdrBPM, fdrWPM, fdrPATH, bpmindfile, 0.4)
    pmap.draw_map(project_dir,fdrcut,resultsfile,BPM_group_tmp,WPM_group_tmp,PATH_group_tmp) ## for drawing graph of BPMs/WPMs/PATHs

    if not ind_bpm.empty:

        ### resuls for redundant groups
        p1,p2,fdr_c,fdrs,g = [],[],[],[],[]
        for i in range(len(BPM_group_tmp)):
            fdr_th = (i+1) * 0.05
            bpms = fdrBPM[fdrBPM['bpm2'] <= fdr_th ]
            bpms = bpms.sort_values(kind='stable',by='bpm2')
            for x in range(bpms.shape[0]):
                fdr_c.append(fdr_th)
                fdrs.append(bpms['bpm2'].iloc[x])
                g.append(int(BPM_group_tmp[i][x]))
                p1.append(path1.loc[bpms.iloc[x].name]['path1names'])
                p2.append(path2.loc[bpms.iloc[x].name]['path2names'])

        d = {'Pathway1': p1, 'Pathway2': p2, 'FDR cut': fdr_c, 'FDR': fdrs, 'Group': g }
        output_BPM_groups_table = pd.DataFrame(data=d)


        bpms = fdrBPM[fdrBPM['bpm2'] <= fdrcut ]
        bpms = bpms.sort_values(kind='stable',by='bpm2')
        bpm_groups = np.array(BPM_group_tmp[-1],dtype=np.int64)
        path1 = path1.reindex(index=bpms.index)
        path1 = path1.reset_index(drop = True)
        path2 = path2.reindex(index=bpms.index)
        path2 = path2.reset_index(drop = True)
        fdrBPM = fdrBPM.reindex(index=bpms.index)
        fdrBPM = fdrBPM.reset_index(drop = True)
        bpm_size = bpm_size.reindex(index=bpms.index)
        bpm_size = bpm_size.reset_index(drop = True)
        bpm_pv_discovery = bpm_pv_discovery.reindex(index=bpms.index)
        bpm_pv_discovery = bpm_pv_discovery.reset_index(drop = True)
        bpm_ranksum_discovery= bpm_ranksum_discovery.reindex(index=bpms.index)
        bpm_ranksum_discovery = bpm_ranksum_discovery.reset_index(drop = True)
        group_numbers = pd.DataFrame(data=bpm_groups,columns=['group'])


        fdrBPM = round(fdrBPM, 2)
        bpm_ranksum_discovery = round(bpm_ranksum_discovery, 2)
        results = [path1,path2,group_numbers,fdrBPM,eff_bpm,bpm_size,bpm_pv_discovery,bpm_ranksum_discovery,bpm_path1_drivers,bpm_path2_drivers]
        output_bpm_table = pd.concat(results, axis=1, sort=False)


    if not ind_wpm.empty:
        ### resuls for redundant groups
        p1,fdr_c,fdrs,g = [],[],[],[]
        for i in range(len(WPM_group_tmp)):
            fdr_th = (i+1) * 0.05
            wpms = fdrWPM[fdrWPM['wpm2'] <= fdr_th ]
            wpms = wpms.sort_values(kind='stable',by='wpm2')
            for x in range(wpms.shape[0]):
                fdr_c.append(fdr_th)
                fdrs.append(wpms['wpm2'].iloc[x])
                g.append(int(WPM_group_tmp[i][x]))
                p1.append(path_wpm.loc[wpms.iloc[x].name]['pathway'])
                
        d = {'Pathway': p1, 'FDR cut': fdr_c, 'FDR': fdrs, 'Group': g }
        output_WPM_groups_table = pd.DataFrame(data=d)

        wpms = fdrWPM[fdrWPM['wpm2'] <= fdrcut ]
        wpms = wpms.sort_values(kind='stable',by='wpm2')
        wpm_groups = np.array(WPM_group_tmp[-1],dtype=np.int64)
        path_wpm = path_wpm.reindex(index=wpms.index)
        path = path_wpm.reset_index(drop = True)
        fdrWPM = fdrWPM.reindex(index=wpms.index)
        fdrWPM = fdrWPM.reset_index(drop = True)
        wpm_size = wpm_size.reindex(index=wpms.index)
        wpm_size = wpm_size.reset_index(drop = True)
        wpm_pv_discovery = wpm_pv_discovery.reindex(index=wpms.index)
        wpm_pv_discovery = wpm_pv_discovery.reset_index(drop = True)
        wpm_ranksum_discovery = wpm_ranksum_discovery.reindex(index=wpms.index)
        wpm_ranksum_discovery = wpm_ranksum_discovery.reset_index(drop = True)
        group_numbers = pd.DataFrame(data=wpm_groups,columns=['group'])


        fdrWPM = round(fdrWPM, 2)
        wpm_ranksum_discovery = round(wpm_ranksum_discovery, 2)
        results = [path,group_numbers,fdrWPM,eff_wpm,wpm_size,wpm_pv_discovery,wpm_ranksum_discovery,wpm_path_drivers]
        output_wpm_table = pd.concat(results, axis=1, sort=False)

    if not ind_path.empty:

        ### resuls for redundant groups
        p1,fdr_c,fdrs,g = [],[],[],[]
        for i in range(len(PATH_group_tmp)):
            fdr_th = (i+1) * 0.05
            path_results = fdrPATH[fdrPATH['path2'] <= fdr_th ]
            path_results = path_results.sort_values(kind='stable',by='path2')
            for x in range(path_results.shape[0]):
                fdr_c.append(fdr_th)
                fdrs.append(path_results['path2'].iloc[x])
                g.append(int(PATH_group_tmp[i][x]))
                p1.append(path_path.loc[path_results.iloc[x].name]['pathway'])
                

        d = {'Pathway': p1, 'FDR cut': fdr_c, 'FDR': fdrs, 'Group': g }
        output_PATH_groups_table = pd.DataFrame(data=d)

        paths = fdrPATH[fdrPATH['path2'] <= fdrcut ]
        paths = paths.sort_values(kind='stable',by='path2')
        path_groups = np.array(PATH_group_tmp[-1],dtype=np.int64)
        path_path = path_path.reindex(index=paths.index)
        path = path_path.reset_index(drop = True)
        fdrPATH = fdrPATH.reindex(index=paths.index)
        fdrPATH = fdrPATH.reset_index(drop = True)
        path_size = path_size.reindex(index=paths.index)
        path_size = path_size.reset_index(drop = True)
        path_pv_discovery = path_pv_discovery.reindex(index=paths.index)
        path_pv_discovery = path_pv_discovery.reset_index(drop = True)
        path_ranksum_discovery = path_ranksum_discovery.reindex(index=paths.index)
        path_ranksum_discovery = path_ranksum_discovery.reset_index(drop = True)
        group_numbers = pd.DataFrame(data=path_groups,columns=['group'])


        fdrPATH = round(fdrPATH, 2)
        path_ranksum_discovery = round(path_ranksum_discovery, 2)
        results = [path,group_numbers,fdrPATH,eff_path,path_size,path_pv_discovery,path_ranksum_discovery]
        output_path_table = pd.concat(results, axis=1, sort=False)


    output_discovery_summary, output_noRD_discovery_summary = [], []

    output_discovery_summary.append([fdrBPM['bpm2'].min(), (fdrBPM['bpm2']<=0.05).sum(), (fdrBPM['bpm2']<=0.1).sum(), (fdrBPM['bpm2']<=0.15).sum(), (fdrBPM['bpm2']<=0.2).sum(), (fdrBPM['bpm2']<=0.25).sum(), (fdrBPM['bpm2']<=0.3).sum(), (fdrBPM['bpm2']<=0.35).sum(), (fdrBPM['bpm2']<=0.4).sum()])

    if not fdrWPM.empty:
        output_discovery_summary.append([fdrWPM['wpm2'].min(), (fdrWPM['wpm2']<=0.05).sum(), (fdrWPM['wpm2']<=0.1).sum(), (fdrWPM['wpm2']<=0.15).sum(), (fdrWPM['wpm2']<=0.2).sum(), (fdrWPM['wpm2']<=0.25).sum(), (fdrWPM['wpm2']<=0.3).sum(), (fdrWPM['wpm2']<=0.35).sum(), (fdrWPM['wpm2']<=0.4).sum()])
        output_discovery_summary.append([fdrPATH['path2'].min(), (fdrPATH['path2']<=0.05).sum(), (fdrPATH['path2']<=0.1).sum(), (fdrPATH['path2']<=0.15).sum(), (fdrPATH['path2']<=0.2).sum(), (fdrPATH['path2']<=0.25).sum(), (fdrPATH['path2']<=0.3).sum(), (fdrPATH['path2']<=0.35).sum(), (fdrPATH['path2']<=0.4).sum()])


    output_header = ['minfdr','fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40']


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
    final_filename = project_dir + "/output_" +  final_filename.rsplit('.',1)[0] + ".xls"


    #### DOUBLE CHECK SORTING WHEN FIXING BPMSIM AND PATHSIM. ####
    if not ind_bpm.empty:
        output_bpm_table.rename(columns={"bpm2": "fdrBPM"}, inplace = True)
        output_bpm_table.rename(columns={"indsize": "bpm_size"}, inplace = True)
        output_bpm_table.sort_values(by=['fdrBPM', 'bpm_pv', 'bpm_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_bpm_table = output_bpm_table.reset_index(drop = True)

    if not ind_wpm.empty:
        output_wpm_table.rename(columns={"wpm2": "fdrWPM"}, inplace = True)
        output_wpm_table.rename(columns={"indsize": "wpm_size"}, inplace = True)
        output_wpm_table.sort_values(by=['fdrWPM', 'wpm_pv', 'wpm_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_wpm_table = output_wpm_table.reset_index(drop = True)

    if  not ind_path.empty:
        output_path_table.rename(columns={"path2": "fdrPATH"}, inplace = True)
        output_path_table.rename(columns={"indsize": "path_size"}, inplace = True)
        output_path_table.sort_values(by=['fdrPATH', 'path_pv', 'path_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_path_table = output_path_table.reset_index(drop = True)

    writer = pd.ExcelWriter(final_filename)

    output_discovery_summary.to_excel(writer, sheet_name='output_discovery_summary')
    output_noRD_discovery_summary.to_excel(writer, sheet_name='output_noRD_discovery_summary')

    if not ind_bpm.empty:
        output_bpm_table.to_excel(writer, sheet_name='output_bpm_table')
        #output_BPM_groups_table.to_excel(writer, sheet_name='BPM_redundant_groups')
    if not ind_wpm.empty:
        output_wpm_table.to_excel(writer, sheet_name='output_wpm_table')
        #output_WPM_groups_table.to_excel(writer, sheet_name='WPM_redundant_groups')
    if not ind_path.empty:
        output_path_table.to_excel(writer, sheet_name='output_path_table')
        #output_PATH_groups_table.to_excel(writer, sheet_name='PATH_redundant_groups')

    writer.save()

    return


