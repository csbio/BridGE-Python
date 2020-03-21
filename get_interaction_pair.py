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

def hygetest_caller(input_row):
	return ht.hygetest(input_row[0],input_row[1],input_row[2],input_row[3])


def get_interaction_pair(pathname1,pathname2,bpmfile,snp2genefile,snp2pathwayfile,effect,ssmfile):
	pklin = open(snp2genefile,'rb')
	snp2gene = pickle.load(pklin)
	#print(snp2gene)
	pklin = open(snp2pathwayfile,'rb')
	snp2path = pickle.load(pklin)
	#print(snp2path.spmatrix)

	## only for testing 
	geneset_file = 'data/toy.genesets.pkl'
	pklin = open(geneset_file,'rb')
	geneset = pickle.load(pklin)

	pklin = open(bpmfile,'rb')
	bpm = pickle.load(pklin)

	## load snp data
	pklin = open('data/SNPdataAD.pkl','rb')
	snpdataAD = pickle.load(pklin)

	pklin = open('data/SNPdataAR.pkl','rb')
	snpdataAR = pickle.load(pklin)

	#load interaction network
	pklin = open(ssmfile,'rb')
	int_network = pickle.load(pklin)

	if effect == 'protective':
		ssm = int_network.protective
	else:
		ssm = int_network.risk

	ssm_pr = ssmfile.split('.')
	model = ssm_pr[0].split('_')[-2]



	#print(bpm.wpm['pathway'][0][0].shape[0])
	## find pathway index in wpm and find snp ids
	p_id1 = 0
	p_id2 = 0 
	pathway_size = bpm.wpm['pathway'][0][0].shape[0]
	for i in range(pathway_size):
		p = bpm.wpm['pathway'][0][0][i][0]
		if p == pathname1:
			p_id1 = i
		if p == pathname2:
			p_id2 = i
		if p_id1 != 0 and p_id2 != 0:
			break

	## find pair bpm id
	n_pairs = int(pathway_size * (pathway_size-1) / 2)
	#print(n_pairs)
	bpm_id = 0
	bpm_path1 = bpm.bpm['path1idx'][0][0]
	bpm_path2 = bpm.bpm['path2idx'][0][0]
	bpm_path1 = bpm_path1 - 1
	bpm_path2 = bpm_path2 - 1

	if p_id1 != p_id2:
		for i in range(n_pairs):
			p1 = bpm_path1[i]
			if p1 == p_id1:
				p2 = bpm_path2[i]
				if p2 == p_id2:
					bpm_id = i
					break

		## get snp ids
		#print(bpm.bpm['ind1'][0][0][1].shape)
		bpm_ind1 = bpm.bpm['ind1'][0][0]
		bpm_ind2 = bpm.bpm['ind2'][0][0]

		ind1 = np.reshape(bpm_ind1[bpm_id],bpm_ind1[bpm_id].shape[0])
		ind2 = np.reshape(bpm_ind2[bpm_id],bpm_ind2[bpm_id].shape[0])
		ind1 = ind1 - 1
		ind2 = ind2 - 1
	else:
		wpm_ind1 = bpm.wpm['ind'][0][0]
		ind1 = np.reshape(wpm_ind1[p_id1],wpm_ind1[p_id1].shape[0])
		ind1 = ind1 - 1
		ind2 = ind1

	## get snp rsids
	#print(snpdataAR.rsid)
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
	#tmp.loc['rs11607758','ACAT1'] = 0
	#tmp = tmp.sum(axis=0)
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
	bin_int_index = np.argwhere(ssm_dis>0.2)
	i = bin_int_index[:,0]
	j = bin_int_index[:,1]
	#i = np.unique(i)
	#j = np.unique(j)
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
		if model == 'RR':
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
		elif model == 'DD':
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
		elif model == 'RD':
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
	t = {'snp1': snps1, 'gene1': genes1, 'snp2': snps2 , 'gene2': genes2,'GI type': GT_type,'case frequency': freq_case,'control frequency': freq_control,'GI': GI}
	output_pair = pd.DataFrame(data=t,columns=['snp1', 'gene1', 'snp2', 'gene2', 'GI type', 'case frequency', 'control frequency', 'GI' ])
	output_pair.sort_values('GI', ascending=False ,inplace=True)

	## preparing output for pathway 1
	ind1_mean_GI  = np.sum((ssm_dis>0.2),axis=1) / ind2.shape[0]
	ind1_mean_GI_bg = np.sum((ssm[ind1,:]>0.2) ,axis=1) / ssm.shape[1]
	snps = ind1_snp
	genes = ind1_gene
	snp_mean_gi = ind1_mean_GI
	snp_mean_gi_bg = ind1_mean_GI_bg
	gi_fold = np.divide(snp_mean_gi,snp_mean_gi_bg)
	gi_fold[np.isnan(gi_fold)] = 0
	in_int = np.sum((ssm_dis>0.2),axis=1)
	in_int = np.reshape(in_int,(in_int.shape[0],1))
	all_int = np.sum((ssm[ind1,:]>0.2) ,axis=1)
	all_int = np.reshape(all_int,(all_int.shape[0],1))
	N = np.broadcast_to(ssm.shape[1],in_int.shape)
	D = np.broadcast_to(ind2.shape[0],in_int.shape)
	hype_in = np.concatenate((N,D,in_int,all_int),axis=1)
	gi_hyge = np.apply_along_axis(hygetest_caller,1,hype_in)

	t = {'snps': snps,'genes': genes,'snp_mean_gi': snp_mean_gi, 'snp_mean_gi_bg': snp_mean_gi_bg, 'gi_fold': gi_fold, 'gi_hyge': gi_hyge}
	output_path1_snp = pd.DataFrame(data=t,columns=['snps','genes', 'snp_mean_gi','snp_mean_gi_bg','gi_fold','gi_hyge'])
	output_path1_snp.sort_values('gi_fold', ascending=False ,inplace=True)
	output_path1_snp =  output_path1_snp.loc[output_path1_snp['gi_fold']>1.2]
	output_path1_snp.sort_values('gi_hyge', ascending=False ,inplace=True)
	
	if p_id1 == p_id2:
		output_path2_snp = NaN
	else:
		## preparing output for pathway 2
		ind2_mean_GI  = np.sum((ssm_dis>0.2),axis=0) / ind1.shape[0]
		ind2_mean_GI_bg = np.sum((ssm[:,ind2]>0.2) ,axis=0) / ssm.shape[1]
		snps = ind2_snp
		genes = ind2_gene
		snp_mean_gi = ind2_mean_GI
		snp_mean_gi_bg = ind2_mean_GI_bg
		gi_fold = np.divide(snp_mean_gi,snp_mean_gi_bg)
		gi_fold[np.isnan(gi_fold)] = 0
		in_int = np.sum((ssm_dis>0.2),axis=0)
		in_int = np.reshape(in_int,(in_int.shape[0],1))
		all_int = np.sum((ssm[:,ind2]>0.2) ,axis=0)
		all_int = np.reshape(all_int,(all_int.shape[0],1))
		N = np.broadcast_to(ssm.shape[1],in_int.shape)
		D = np.broadcast_to(ind1.shape[0],in_int.shape)
		hype_in = np.concatenate((N,D,in_int,all_int),axis=1)
		gi_hyge = np.apply_along_axis(hygetest_caller,1,hype_in)

		t = {'snps': snps,'genes': genes,'snp_mean_gi': snp_mean_gi, 'snp_mean_gi_bg': snp_mean_gi_bg, 'gi_fold': gi_fold, 'gi_hyge': gi_hyge}
		output_path2_snp = pd.DataFrame(data=t,columns=['snps','genes', 'snp_mean_gi','snp_mean_gi_bg','gi_fold','gi_hyge'])
		output_path2_snp.sort_values('gi_fold', ascending=False ,inplace=True)
		output_path2_snp =  output_path2_snp.loc[output_path2_snp['gi_fold']>1.2]
		output_path2_snp.sort_values('gi_hyge', ascending=False ,inplace=True)
	out_pair = [output_path1_snp,output_path2_snp]
	return out_pair













	
	

	


	


snp2genefile = 'data/snpgenemapping_50kb.pkl'
snp2pathwayfile = 'data/snp_pathway_min10_max300.pkl'
bpmfile = 'data/bpmind.pkl'
pathname1 = 'KEGG_TERPENOID_BACKBONE_BIOSYNTHESIS'
pathname2 = 'BIOCARTA_G2_PATHWAY'
ssmfile = 'data/ssM_hygessi_RD_R0.pkl'
effect = 'protective'
x = get_interaction_pair(pathname1,pathname2,bpmfile,snp2genefile,snp2pathwayfile,effect,ssmfile)
print(x[1])