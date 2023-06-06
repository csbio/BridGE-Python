import numpy as np
import math
import pandas as pd
import pickle
from scipy.io import loadmat
from scipy.stats import hypergeom
from scipy.io import loadmat,savemat
from scipy.stats import chi2_contingency
from scipy.stats import mannwhitneyu
from scipy.stats import rankdata
from scipy.stats import norm
from classes import bpmindclass as bpmc
from classes import snpsetclass as snpp
from classes import SNPdataclass
from classes import InteractionNetwork
from classes import genesetdataclass
from corefuns import hygetest as ht
np.seterr(divide='ignore', invalid='ignore')


# get_interaction_pair() finds driver SNPs and genes for BPMs, WPMs. Figures out BPM or WPM by comparing path1 and path2 names.
#
# INPUTS:
#   n: Number of BPMs/WPMs provided to the function.
#	path1: Names of pathway1 in BPMs/WPMs
#	path2: Names of pathway2 in BPMs (same as path1 for WPMs)
#	effects: array of 'risk'/'protective' indicating the interaction type
#   bpmfile: file containing SNP ids for BPM/WPMs in pickle format.
#   snp2pathwayfile: Pickle file containing mapping of SNPs to pathways.
#   snp2genefile: Pickle file containing mapping of SNPs to Genes.
#	path_ids: a map of pathway names to their id in the dataset, made in collectresults function
#
# OUTPUTS:
#   returns an array containing 3 data frames
#       - bpm_path1_drivers: driver SNPs(abd their mapped gene) for pathway1 in BPM
#       - bpm_path2_drivers: driver SNPs(abd their mapped gene) for pathway2 in BPM
#       - wpm_path_drivers: driver SNPs(abd their mapped gene) for pathway with WPM
# 



def hygetest_caller(input_row):
	return ht.hygetest(input_row[0],input_row[1],input_row[2],input_row[3])


def get_interaction_pair(n,path1,path2,effects,ssmfile,bpmfile,snp2pathwayfile,snp2genefile,path_ids,densitycutoff=None):
	pklin = open(snp2genefile,'rb')
	snp2gene = pickle.load(pklin)
	pklin = open(snp2pathwayfile,'rb')
	snp2path = pickle.load(pklin)



	

	pklin = open(bpmfile,'rb')
	bpm = pickle.load(pklin)

	# Find project dir
	dir_sp = snp2pathwayfile.split('/')
	s = '/'
	project_dir = s.join(dir_sp[:-1])


	geneset_file = snp2path.geneset
	pklin = open(geneset_file,'rb')
	geneset = pickle.load(pklin)


	## load snp data
	pklin = open(project_dir+'/SNPdataAD.pkl','rb')
	snpdataAD = pickle.load(pklin)

	pklin = open(project_dir+'/SNPdataAR.pkl','rb')
	snpdataAR = pickle.load(pklin)

	#load interaction network
	pklin = open(ssmfile,'rb')
	int_network = pickle.load(pklin)

	# find score cutoffs if densitycutoff provided
	if densitycutoff == None:
		pos_cutoff = 0.2
		neg_cutoff = 0.2
	else:
		if densitycutoff <= 0 or densitycutoff >=1:
			densitycutoff = 0.1
		pos_cutoff = np.quantile(int_network.protective,1-densitycutoff)
		neg_cutoff = np.quantile(int_network.risk,1-densitycutoff)


	bpm_path1_drivers, bpm_path2_drivers = [], []
	wpm_path_drivers = []
	path_index = []

	for i_path in range(n):
		pathname1 = path1.iloc[i_path][0]
		pathname2 = path2.iloc[i_path][0]
		effect = effects.iloc[i_path][0]


		if effect == 'protective':
			ssm = int_network.protective
			max_id = int_network.protective_max_id
			score_cutoff = pos_cutoff
		else:
			ssm = int_network.risk
			max_id = int_network.risk_max_id
			score_cutoff = neg_cutoff

		ssm_pr = ssmfile.split('.')
		model = ssm_pr[0].split('_')[-2]




		## find pathway index in wpm and find snp ids
		p_id1 = 0
		p_id2 = 0 
		pathway_size = bpm.wpm['pathway'].shape[0]
		p_id1 = path_ids[pathname1]
		p_id2 = path_ids[pathname2]
		if p_id2 < p_id1:
			tmp = p_id1
			p_id1 = p_id2
			p_id2 = tmp

		## find pair bpm id
		n_pairs = int(pathway_size * (pathway_size-1) / 2)
		bpm_id = 0

		if p_id1 != p_id2:
			#find bpm_id
			bpm_id = int(n_pairs - (pathway_size-p_id1)*(pathway_size-p_id1-1)/2 + p_id2 - p_id1 - 1)
			if effect == 'protective':
				rem = 0
			else:
				rem = bpm.bpm.shape[0]
			path_index.append(bpm_id+rem)
			## get snp ids
			bpm_ind1 = bpm.bpm['ind1']
			bpm_ind2 = bpm.bpm['ind2']

			ind1 = bpm_ind1[bpm_id]
			ind2 = bpm_ind2[bpm_id]
		else:
			wpm_ind1 = bpm.wpm['ind']
			ind1 = wpm_ind1[p_id1]
			ind2 = ind1
			if effect == 'protective':
				rem = 0
			else:
				rem = bpm.wpm.shape[0]
			path_index.append(p_id1+rem)

		## get snp rsids
		ind1_snp = snpdataAR.rsid[ind1].values
		ind2_snp = snpdataAR.rsid[ind2].values

		## find all genes for the 2 pathways from the geneset(index)
		tmp = geneset.gpmatrix[pathname1]
		ind1_gp = tmp[tmp==1].index.values
		tmp = geneset.gpmatrix[pathname2]
		ind2_gp = tmp[tmp==1].index.values

		## find genes common between snp and pathways(excluding multiple gene membership for pathways)
		## first remove genes in pathways which are not in the snp2gene matrix
		all_genes = snp2gene.columns.values
		ind1_gp_updated = np.intersect1d(ind1_gp,all_genes)
		ind2_gp_updated = np.intersect1d(ind2_gp,all_genes)

		## remove snps in pathways which are not in the snp2gene
		all_snps = snp2gene.index.values
		ind1_snp_updated = np.intersect1d(ind1_snp,all_snps)
		ind2_snp_updated = np.intersect1d(ind2_snp,all_snps)
		tmp = snp2gene.loc[ind1_snp_updated,ind1_gp_updated]
		ind1_gene = []
		for snp in ind1_snp:
			if not snp in ind1_snp_updated:
				ind1_gene.append('/')
				continue
			t1 = tmp.loc[snp]
			ind1_gene_tmp = t1[t1==1].index.values
			if len(ind1_gene_tmp) > 0:
				ind1_gene.append('/'.join(ind1_gene_tmp))
			else:
				ind1_gene.append('/')
		if p_id1 == p_id2:
			ind2_gene = ind1_gene
		else:
			tmp = snp2gene.loc[ind2_snp_updated,ind2_gp_updated]
			ind2_gene = []
			for snp in ind2_snp:
				if not snp in ind2_snp_updated:
					ind2_gene.append('/')
					continue
				t1 = tmp.loc[snp]
				ind2_gene_tmp = t1[t1==1].index.values
				if len(ind2_gene_tmp) > 0:
					ind2_gene.append('/'.join(ind2_gene_tmp))
				else:
					ind2_gene.append('/')

		ind1_gene = np.array(ind1_gene)
		ind2_gene = np.array(ind2_gene)
		ssm_dis = ssm[ind1,:][:,ind2]
		if model == 'combined': 
			m_id = max_id[ind1,:][:,ind2]
		if p_id1 == p_id2:
			tmp_dis = np.tril(ssm_dis)
		else:
			tmp_dis = ssm_dis
		bin_int_index = np.argwhere(tmp_dis>score_cutoff)
		i = bin_int_index[:,0]
		j = bin_int_index[:,1]
		snps1 = ind1_snp[i]
		snps2 = ind2_snp[j]
		genes1 = ind1_gene[i]
		genes2 = ind2_gene[j]

		interaction_pairs = i.shape[0]
		GI = np.zeros(interaction_pairs)
		GT_type = []
		freq_case = np.zeros(interaction_pairs)
		freq_control = np.zeros(interaction_pairs)
		pheno = snpdataAD.pheno
		pheno_size = pheno.shape[0]
		snpdataAD1 = snpdataAD.data.iloc[:,ind1]
		snpdataAD2 = snpdataAD.data.iloc[:,ind2]
		snpdataAR1 = snpdataAR.data.iloc[:,ind1]
		snpdataAR2 = snpdataAR.data.iloc[:,ind2]
		snp_pair = np.zeros((pheno_size,interaction_pairs))

	
		pheno_nnz = np.sum(pheno)
		pheno_nz = pheno_size - pheno_nnz
		I = np.ones(pheno.shape)
		res_pheno = I - pheno
	
		for k in range(interaction_pairs):
			GI[k] = ssm_dis[i[k],j[k]]
			tmp_model = model
			if model == 'combined':
				dm = m_id[i[k],j[k]] # disease model with maximum interaction in combined model
				if dm == 1: ## RR	
					tmp_model = 'RR'
				elif dm == 2: # DD
					tmp_model = 'DD'
				elif dm == 3: # RD
					tmp_model = 'RD'
			if tmp_model == 'RR':
				GT_type.append('recessive')
				## case frequency on 1-1
				p1 = np.multiply(pheno,snpdataAR1.iloc[:,i[k]])
				p2 = np.multiply(pheno,snpdataAR2.iloc[:,j[k]])
				s = np.multiply(p1,p2)
				s = np.sum(s)
				freq_case[k] = s/pheno_nnz
				## control frequency on 1-1
				p1 = np.multiply(res_pheno,snpdataAR1.iloc[:,i[k]])
				p2 = np.multiply(res_pheno,snpdataAR2.iloc[:,j[k]])
				s = np.multiply(p1,p2)
				s = np.sum(s)
				freq_control[k] = s/pheno_nz
				snp_pair[:,k] = np.multiply(snpdataAR1.iloc[:,i[k]],snpdataAR2.iloc[:,j[k]])
			elif tmp_model == 'DD':
				GT_type.append('dominant')
				## case frequency on 1-1
				p1 = np.multiply(pheno,snpdataAD1.iloc[:,i[k]])
				p2 = np.multiply(pheno,snpdataAD2.iloc[:,j[k]])
				s = np.multiply(p1,p2)
				s = np.sum(s)
				freq_case[k] = s/pheno_nnz
				## control frequency on 1-1
				p1 = np.multiply(res_pheno,snpdataAD1.iloc[:,i[k]])
				p2 = np.multiply(res_pheno,snpdataAD2.iloc[:,j[k]])
				s = np.multiply(p1,p2)
				s = np.sum(s)
				freq_control[k] = s/pheno_nz
				snp_pair[:,k] = np.multiply(snpdataAD1.iloc[:,i[k]],snpdataAD2.iloc[:,j[k]])
			elif tmp_model == 'RD':
				GT_type.append('recessive_dominant')
				## case frequency on 1-1 RD
				p1 = np.multiply(pheno,snpdataAR1.iloc[:,i[k]])
				p2 = np.multiply(pheno,snpdataAD2.iloc[:,j[k]])
				s = np.multiply(p1,p2)
				s = np.sum(s)
				freq_case_1 = s/pheno_nnz
				## control frequency on 1-1 RD
				p1 = np.multiply(res_pheno,snpdataAR1.iloc[:,i[k]])
				p2 = np.multiply(res_pheno,snpdataAD2.iloc[:,j[k]])
				s = np.multiply(p1,p2)
				s = np.sum(s)
				freq_control_1 = s/pheno_nz
				## case frequency on 1-1 DR
				p1 = np.multiply(pheno,snpdataAD1.iloc[:,i[k]])
				p2 = np.multiply(pheno,snpdataAR2.iloc[:,j[k]])
				s = np.multiply(p1,p2)
				s = np.sum(s)
				freq_case_2 = s/pheno_nnz
				## control frequency on 1-1 DR
				p1 = np.multiply(res_pheno,snpdataAD1.iloc[:,i[k]])
				p2 = np.multiply(res_pheno,snpdataAR2.iloc[:,j[k]])
				s = np.multiply(p1,p2)
				s = np.sum(s)
				freq_control_2 = s/pheno_nz

				if freq_control_1 == 0:
					freq_control_1 = 1
				if freq_control_2 == 0:
					freq_control_2 = 1
				freq_R1 = freq_case_1/freq_control_1
				freq_R2 = freq_case_2/freq_control_2
				if effect == 'protective':
					if freq_R1 > freq_R2:
						freq_case[k] = freq_case_2
						freq_control[k] = freq_control_2
						snp_pair[:,k] = np.multiply(snpdataAD1.iloc[:,i[k]],snpdataAR2.iloc[:,j[k]])
					else:
						freq_case[k] = freq_case_1
						freq_control[k] = freq_control_1
						snp_pair[:,k] = np.multiply(snpdataAR1.iloc[:,i[k]],snpdataAD2.iloc[:,j[k]])
				else:
					if freq_R1 > freq_R2:
						freq_case[k] = freq_case_1
						freq_control[k] = freq_control_1
						snp_pair[:,k] = np.multiply(snpdataAR1.iloc[:,i[k]],snpdataAD2.iloc[:,j[k]])
					else:
						freq_case[k] = freq_case_2
						freq_control[k] = freq_control_2
						snp_pair[:,k] = np.multiply(snpdataAD1.iloc[:,i[k]],snpdataAR2.iloc[:,j[k]])


		if model == 'RR' or model == 'RD' or model == 'DD' or model == 'combined':
			t = {'snp1': snps1, 'gene1': genes1, 'snp2': snps2 , 'gene2': genes2,'GI type': GT_type,'case frequency': freq_case,'control frequency': freq_control,'GI': GI}
			output_pair = pd.DataFrame(data=t,columns=['snp1', 'gene1', 'snp2', 'gene2', 'GI type', 'case frequency', 'control frequency', 'GI' ])
			output_pair.sort_values('GI', ascending=False ,inplace=True)

		ind1 = np.array(ind1)
		ind2 = np.array(ind2)

		## preparing output for pathway 1
		ind1_mean_GI  = np.sum((ssm_dis>score_cutoff),axis=1) / ind2.shape[0]
		ind1_mean_GI_bg = np.sum((ssm[ind1,:]>score_cutoff) ,axis=1) / ssm.shape[1]
		snps = ind1_snp
		genes = ind1_gene
		snp_mean_gi = ind1_mean_GI
		snp_mean_gi_bg = ind1_mean_GI_bg
		gi_fold = np.divide(snp_mean_gi,snp_mean_gi_bg)
		gi_fold[np.isnan(gi_fold)] = 0
		in_int = np.sum((ssm_dis>score_cutoff),axis=1)
		in_int = np.reshape(in_int,(in_int.shape[0],1))
		all_int = np.sum((ssm[ind1,:]>score_cutoff) ,axis=1)
		all_int = np.reshape(all_int,(all_int.shape[0],1))
		N = np.broadcast_to(ssm.shape[1],in_int.shape)
		D = np.broadcast_to(ind2.shape[0],in_int.shape)
		hype_in = np.concatenate((N,D,in_int,all_int),axis=1)
		gi_hyge = np.apply_along_axis(hygetest_caller,1,hype_in)
		gi_hyge_log	= -1 * np.log10(gi_hyge)


		t = {'snps': snps,'genes': genes,'snp_mean_gi': snp_mean_gi, 'snp_mean_gi_bg': snp_mean_gi_bg, 'gi_fold': gi_fold, 'gi_hyge': gi_hyge,'gi_hyge_log': gi_hyge_log}
		output_path1_snp = pd.DataFrame(data=t,columns=['snps','genes', 'snp_mean_gi','snp_mean_gi_bg','gi_fold','gi_hyge','gi_hyge_log'])
		output_path1_snp = output_path1_snp[(output_path1_snp['gi_hyge']<=0.05)&(output_path1_snp['gi_fold'] > 1.2)]
		output_path1_snp.sort_values('gi_fold', ascending=False ,inplace=True)
		#if output_path1_snp.shape[0] > 20:
		#	output_path1_snp = output_path1_snp.iloc[0:20,:]
		#output_path1_snp =  output_path1_snp.loc[output_path1_snp['gi_fold']>1.2]
		#output_path1_snp.sort_values('gi_hyge', ascending=False ,inplace=True)
	
		if p_id1 == p_id2:
			output_path2_snp = None
		else:
			## preparing output for pathway 2
			ind2_mean_GI  = np.sum((ssm_dis>score_cutoff),axis=0) / ind1.shape[0]
			ind2_mean_GI_bg = np.sum((ssm[ind2,:]>score_cutoff) ,axis=1) / ssm.shape[0]
			snps = ind2_snp
			genes = ind2_gene
			snp_mean_gi = ind2_mean_GI
			snp_mean_gi_bg = ind2_mean_GI_bg
			gi_fold = np.divide(snp_mean_gi,snp_mean_gi_bg)
			gi_fold[np.isnan(gi_fold)] = 0
			in_int = np.sum((ssm_dis>score_cutoff),axis=0)
			in_int = np.reshape(in_int,(in_int.shape[0],1))
			all_int = np.sum((ssm[:,ind2]>score_cutoff) ,axis=0)
			all_int = np.reshape(all_int,(all_int.shape[0],1))
			N = np.broadcast_to(ssm.shape[1],in_int.shape)
			D = np.broadcast_to(ind1.shape[0],in_int.shape)
			hype_in = np.concatenate((N,D,in_int,all_int),axis=1)
			gi_hyge = np.apply_along_axis(hygetest_caller,1,hype_in)
			gi_hyge_log = -1 * np.log10(gi_hyge)

			t = {'snps': snps,'genes': genes,'snp_mean_gi': snp_mean_gi, 'snp_mean_gi_bg': snp_mean_gi_bg, 'gi_fold': gi_fold, 'gi_hyge': gi_hyge,'gi_hyge_log': gi_hyge_log}
			output_path2_snp = pd.DataFrame(data=t,columns=['snps','genes', 'snp_mean_gi','snp_mean_gi_bg','gi_fold','gi_hyge','gi_hyge_log'])
			output_path2_snp = output_path2_snp[(output_path2_snp['gi_hyge'] <= 0.05) & (output_path2_snp['gi_fold'] > 1.2)]
			output_path2_snp.sort_values('gi_fold', ascending=False ,inplace=True)
			#if output_path2_snp.shape[0] > 20:
			#	output_path2_snp = output_path2_snp.iloc[0:20,:]
			#output_path2_snp =  output_path2_snp.loc[output_path2_snp['gi_fold']>1.2]
			#output_path2_snp.sort_values('gi_hyge', ascending=False ,inplace=True)
		if p_id1 != p_id2:
			idx = output_path1_snp[output_path1_snp['gi_fold']>1].dropna()
			bpm_snp_tmp = []
			if not idx.empty:
				for row in idx.head(20).iterrows():
					temp_str = row[1]['snps']+'_'+row[1]['genes']+'_fold'+str(round(row[1]['gi_fold'], 2))+'_hyge'+str(round(row[1]['gi_hyge_log'], 2))
					bpm_snp_tmp.append(temp_str)
				bpm_path1_drivers.append(';'.join(bpm_snp_tmp))
			else:
				bpm_path1_drivers.append('')
			idx = output_path2_snp[output_path2_snp['gi_fold']>1].dropna()
			bpm_snp_tmp = []
			if not idx.empty:
				for row in idx.head(20).iterrows():
					temp_str = row[1]['snps']+'_'+row[1]['genes']+'_fold'+str(round(row[1]['gi_fold'], 2))+'_hyge'+str(round(row[1]['gi_hyge_log'], 2))
					bpm_snp_tmp.append(temp_str)
				bpm_path2_drivers.append(';'.join(bpm_snp_tmp))
			else:
				bpm_path2_drivers.append('')
		else:
			idx = output_path1_snp[output_path1_snp['gi_fold']>1].dropna()
			wpm_snp_tmp = []
			if not idx.empty:
				for row in idx.head(20).iterrows():
					temp_str = row[1]['snps']+'_'+row[1]['genes']+'_fold'+str(round(row[1]['gi_fold'], 2))+'_hyge'+str(round(row[1]['gi_hyge_log'], 2))
					wpm_snp_tmp.append(temp_str)
				wpm_path_drivers.append(';'.join(wpm_snp_tmp))
			else:
				wpm_path_drivers.append('')



	bpm_path1_drivers = pd.DataFrame(bpm_path1_drivers, columns=['bpm_path1_drivers'],index = path_index)
	bpm_path2_drivers = pd.DataFrame(bpm_path2_drivers, columns=['bpm_path2_drivers'],index = path_index)
	wpm_path_drivers = pd.DataFrame(wpm_path_drivers, columns=['wpm_path_drivers'],index = path_index)

	out_triple = [bpm_path1_drivers,bpm_path2_drivers,wpm_path_drivers]
	return out_triple





	
	

	


	


