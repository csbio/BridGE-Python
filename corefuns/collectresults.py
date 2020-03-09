import pandas as pd
import pickle
from corefuns import check_BPM_WPM_redundancy as cbwr

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

    if not ind_bpm.empty:
        print('bpm_check')
        fdrBPM = fdrBPM.iloc[ind_bpm.index]
        bpm_pv_discovery = bpm_pv.iloc[ind_bpm.index]
        bpm_ranksum_discovery = bpm_ranksum.iloc[ind_bpm.index]
        print(fdrBPM)
        print(bpm_pv_discovery)
        print(bpm_ranksum_discovery)

    if not ind_wpm.empty:
        print('wpm_check')
        fdrWPM = fdrWPM.iloc[ind_wpm.index]
        wpm_pv_discovery = wpm_pv.iloc[ind_wpm.index]
        wpm_ranksum_discovery = wpm_ranksum.iloc[ind_wpm.index]
        print(fdrWPM)
        print(wpm_pv_discovery)
        print(wpm_ranksum_discovery)

    if  not ind_path.empty:
        print('path_check')
        fdrPATH = fdrPATH.iloc[ind_path.index]
        path_pv_discovery = path_pv.iloc[ind_path.index]
        path_ranksum_discovery = path_ranksum.iloc[ind_path.index]
        print(fdrPATH)
        print(path_pv_discovery)
        print(path_ranksum_discovery)

    if not ind_bpm.empty or not ind_wpm.empty or not ind_path.empty:
        pklin = open(bpmindfile,"rb")
        bpm = pickle.load(pklin)
        pklin.close()
        pklin = open(snppathwayfile,"rb")
        snpset = pickle.load(pklin)
        pklin.close()

        genesetfile = "../refdata/" + snpset.geneset.split('/')[-1]
        print(genesetfile)


        # Remove when refactoring for get_interaction_pair()
        header = ['snp', 'genes', 'snp_mean_gi', 'snp_mean_gi_bg', 'gi_fold', 'gi_hyge']

        if not ind_bpm.empty:
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

            eff_bpm = pd.DataFrame(eff_bpm, columns=['eff_bpm'])
            print(eff_bpm)

            for i in range(len(ind_bpm)):
                # [output_path1_snp, output_path2_snp] = get_interaction_pair(path1.iloc[i],path2.iloc[i],eff_bpm.iloc[i],ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile)

                output_path1_snp_lol = [['rs13006682',	'FHL2',	0.178571428571429,	0.0582524271844660,	3.06547619047619,	1.80044249891492],
                ['rs9818',	'RPA1',	0.214285714285714,	0.0841423948220065,	2.54670329670330,	1.68837784881871],
                ['rs11651877',	'RPA1',	0.392857142857143,	0.220064724919094,	1.78518907563025,	1.63164194881376],
                ['rs4851774',	'FHL2',	0.214285714285714,	0.0873786407766990,	2.45238095238095,	1.60783957534689],
                ['rs2048982',	'FHL2',	0.357142857142857,	0.200647249190939,	1.77995391705069,	1.48560941293369],
                ['rs4724443',	'IGFBP3',	0.107142857142857,	0.0323624595469256,	3.31071428571429,	1.28403118126272],
                ['rs6738207',	'FHL2',	0.178571428571429,	0.0776699029126214,	2.29910714285714,	1.27604005230902],
                ['rs4405588',	'RPA1',	0.142857142857143,	0.0614886731391586,	2.32330827067669,	1.09385874837223],
                ['rs11101282',	'MAPK8',	0.107142857142857,	0.0420711974110032,	2.54670329670330,	0.987752130042409],
                ['rs1048483',	'HIC1',	0.214285714285714,	0.126213592233010,	1.69780219780220,	0.910347974215635],
                ['rs2099339',	'NFKBIB',	0.178571428571429,	0.110032362459547,	1.62289915966387,	0.744923832292046],
                ['rs1914748',	'FHL2',	0.178571428571429,	0.116504854368932,	1.53273809523810,	0.669813755595502],
                ['rs1131636',	'RPA1',	0.500000000000000,	0.414239482200647,	1.20703125000000,	0.655189098945817],
                ['rs1202169',	'ABCB1'	,0.0714285714285714,	0.0323624595469256,	2.20714285714286,	0.644751264382561],
                ['rs2453839',	'IGFBP3',	0.0357142857142857,	0.0129449838187702,	2.75892857142857,	0.498345915720932],
                ['rs3803804',	'HIC1',	0.142857142857143,	0.103559870550162,	1.37946428571429,	0.486338452552602],
                ['rs8067702',	'RPA1'	,0.0714285714285714,	0.0420711974110032,	1.69780219780220,	0.477564371718690],
                ['rs7208143',	'RPA1',	0.178571428571429,	0.142394822006473,	1.25405844155844,	0.436112987224146],
                ['rs1001413',	'NFKBIB',	0.0357142857142857,	0.0161812297734628,	2.20714285714286,	0.420106445628633],
                ['rs11621560',	'HSP90AA1',	0.0357142857142857,	0.0291262135922330,	1.22619047619048,	0.236818616694704]]

                output_path1_snp = pd.DataFrame(output_path1_snp_lol, columns = header)

                output_path2_snp_lol = [['rs8098694',	'BCL2',	0.242424242424242,	0.0938511326860841,	2.58307210031348,	2.21240545634186],
                ['rs1564483',	'BCL2',	0.272727272727273,	0.126213592233010,	2.16083916083916,	1.89157716580341],
                ['rs1541296',	'BCL2',	0.545454545454545,	0.365695792880259,	1.49155269509252,	1.69327220197517],
                ['rs12457893',	'BCL2',	0.424242424242424,	0.268608414239482,	1.57940854326397,	1.51790067614466],
                ['rs11171710',	'CDK2',	0.454545454545455,	0.304207119741100,	1.49419729206963,	1.39965046264644],
                ['rs10860361',	'APAF1',	0.212121212121212,	0.110032362459547,	1.92780748663102,	1.27171986253408],
                ['rs4941185',	'BCL2',	0.151515151515152,	0.0744336569579288,	2.03557312252964,	1.07869510863049],
                ['rs4456611',	'BCL2',	0.0909090909090909,	0.0355987055016181,	2.55371900826446,	0.995554309813982],
                ['rs7297662',	'CDK2',	0.0606060606060606,	0.0226537216828479,	2.67532467532468,	0.781805139376270],
                ['rs4805475',	'CCNE1',	0.0909090909090909,	0.0485436893203884,	1.87272727272727,	0.683253708802447],
                ['rs1381548',	'BCL2',	0.212121212121212,	0.152103559870550,	1.39458413926499,	0.661980486630221],
                ['rs1048691',	'CDK4',	0.303030303030303,	0.239482200647249,	1.26535626535627,	0.618733470606919],
                ['rs10860363',	'APAF1',	0.0909090909090909,	0.0550161812297735,	1.65240641711230,	0.571970140447586],
                ['rs5754248',	'TIMP3',	0.0606060606060606,	0.0323624595469256,	1.87272727272727,	0.537553673625681]]

                output_path2_snp = pd.DataFrame(output_path2_snp_lol, columns = header)

                idx = output_path1_snp[output_path1_snp['gi_fold']>1].dropna()

                bpm_snp_tmp, bpm_path1_drivers, bpm_path2_drivers = [], [], []
                if not idx.empty:
                    for row in idx.head(20).iterrows():
                        temp_str = row[1]['snp']+'_'+row[1]['genes']+'_'+str(round(row[1]['gi_fold'], 2))
                        bpm_snp_tmp.append(temp_str)
                    bpm_path1_drivers.append(';'.join(bpm_snp_tmp))

                else:
                    bpm_path1_drivers.append('')

                idx = output_path2_snp[output_path2_snp['gi_fold']>1].dropna()
                bpm_snp_tmp = []

                if not idx.empty:
                    for row in idx.head(20).iterrows():
                        temp_str = row[1]['snp']+'_'+row[1]['genes']+'_'+str(round(row[1]['gi_fold'], 2))
                        bpm_snp_tmp.append(temp_str)
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
            print(eff_wpm)

            for i in range(len(ind_wpm)):
                # [output_path_snp] = get_interaction_pair(path_wpm.iloc[i],path_wpm.iloc[i],eff_bpm.iloc[i],ssmfile,bpmindfile,snppathwayfile,snpgenemappingfile)
                output_path_snp = pd.DataFrame([[]], columns = header)

                idx = output_path_snp[output_path_snp['gi_fold']>1].dropna()

                bpm_snp_tmp, bpm_path1_drivers, bpm_path2_drivers = [], [], []
                if idx.shape[1] > 0:
                    for row in idx.head(20).iterrows():
                        temp_str = row[1]['snp']+'_'+row[1]['genes']+'_'+str(round(row[1]['gi_fold'], 2))
                        bpm_snp_tmp.append(temp_str)
                    bpm_path1_drivers.append(';'.join(bpm_snp_tmp))

                else:
                    bpm_path1_drivers.append('')

            bpm_path1_drivers = pd.DataFrame(bpm_path1_drivers, columns=['bpm_path1_drivers'])

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
            print(eff_path)

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
        print('output_bpm_table')
        print(output_bpm_table)
        input()

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
        print('output_wpm_table')
        print(output_wpm_table)
        input()

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
        print('output_path_table')
        print(output_path_table)
        input()


    output_discovery_summary, output_noRD_discovery_summary = [], []

    output_discovery_summary.append([fdrBPM['bpm2'].min(), (fdrBPM['bpm2']<=0.05).sum(), (fdrBPM['bpm2']<=0.1).sum(), (fdrBPM['bpm2']<=0.15).sum(), (fdrBPM['bpm2']<=0.2).sum(), (fdrBPM['bpm2']<=0.25).sum(), (fdrBPM['bpm2']<=0.3).sum(), (fdrBPM['bpm2']<=0.35).sum(), (fdrBPM['bpm2']<=0.4).sum()])

    if not fdrWPM.empty:
        output_discovery_summary.append([fdrWPM['wpm2'].min(), (fdrWPM['wpm2']<=0.05).sum(), (fdrWPM['wpm2']<=0.1).sum(), (fdrWPM['wpm2']<=0.15).sum(), (fdrWPM['wpm2']<=0.2).sum(), (fdrWPM['wpm2']<=0.25).sum(), (fdrWPM['wpm2']<=0.3).sum(), (fdrWPM['wpm2']<=0.35).sum(), (fdrWPM['wpm2']<=0.4).sum()])
        output_discovery_summary.append([fdrPATH['path2'].min(), (fdrPATH['path2']<=0.05).sum(), (fdrPATH['path2']<=0.1).sum(), (fdrPATH['path2']<=0.15).sum(), (fdrPATH['path2']<=0.2).sum(), (fdrPATH['path2']<=0.25).sum(), (fdrPATH['path2']<=0.3).sum(), (fdrPATH['path2']<=0.35).sum(), (fdrPATH['path2']<=0.4).sum()])


    # print('== output_discovery_summary ==')
    # for row in output_discovery_summary:
    #     print(row)
    # print('==============================')

    output_header = ['minfdr','fdr05','fdr10','fdr15','fdr20','fdr25','fdr30','fdr35','fdr40']



    BPM_nosig_noRD,WPM_nosig_noRD,PATH_nosig_noRD,BPM_group_tmp,WPM_group_tmp,PATH_group_tmp = cbwr.check_BPM_WPM_redundancy(fdrBPM, fdrWPM, fdrPATH, bpmindfile, 0.4)

    print('==== Group Testing ====')
    print(BPM_group_tmp)
    print(BPM_nosig_noRD)
    print(WPM_group_tmp)
    print(WPM_nosig_noRD)
    print(PATH_group_tmp)
    print(PATH_nosig_noRD)
    print('=======================')
    input()

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

    # print(output_discovery_summary)
    # print(output_noRD_discovery_summary)
    # input()

    final_filename = resultsfile.split('/').pop()
    final_filename = "data/output_" +  final_filename.rsplit('.',1)[0] + ".xls"

    if not ind_bpm.empty:
        output_bpm_table.rename(columns={"bpm2": "fdrBPM"}, inplace = True)
        output_bpm_table.rename(columns={"indsize": "bpm_size"}, inplace = True)
        output_bpm_table.sort_values(by=['fdrBPM', 'bpm_pv', 'bpm_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_bpm_table = output_bpm_table.reset_index(drop = True)
        print(output_bpm_table)
        input()

    if not ind_wpm.empty:
        output_wpm_table.rename(columns={"wpm2": "fdrWPM"}, inplace = True)
        output_wpm_table.rename(columns={"indsize": "wpm_size"}, inplace = True)
        output_wpm_table.sort_values(by=['wpm2', 'wpm_pv', 'wpm_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_wpm_table = output_wpm_table.reset_index(drop = True)
        print(output_wpm_table)
        input()

    if  not ind_path.empty:
        output_path_table.rename(columns={"path2": "fdrPATH"}, inplace = True)
        output_path_table.rename(columns={"indsize": "path_size"}, inplace = True)
        output_path_table.sort_values(by=['fdrPATH', 'path_pv', 'path_ranksum'], ascending=[1, 1, 0], inplace=True)
        output_path_table = output_path_table.reset_index(drop = True)
        print(output_path_table)
        input()

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
