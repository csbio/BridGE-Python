import numpy as np
import pandas as pd
import pickle
from classes import snpsetclass as snps
import sys
import copy



def withinclassrand(seed,plinkCluster,datafile):
	# load data file
	pklin = open(datafile,'rb')
	snpdata = pickle.load(pklin)
	pklin.close()
	# load cluster file
	cluster_data = pd.read_csv(plinkCluster,header=None,sep='\s+',engine='python',dtype=str)
	# match with snp data
	fid = snpdata.fid.astype(str)
	pid = snpdata.pid.astype(str)
	cluster_data = cluster_data[cluster_data[0].isin(fid)]
	# find unique clusters
	cluster_numbers = np.unique(cluster_data[2].to_numpy())
	# process clusters
	clusters = []
	for cl in cluster_numbers:
		ind = cluster_data[cluster_data[2]==cl].index
		if not ind.shape[0]==2:
			sys.exit('error in cluster file') 
		clusters.append([cluster_data.iloc[ind[0],0],cluster_data.iloc[ind[0],1],cl,cluster_data.iloc[ind[1],0],cluster_data.iloc[ind[1],1]])

	clusters = np.array(clusters)
	# randomize
	np.random.seed(seed)
	randidx = np.random.permutation(clusters.shape[0])

	n = int(np.floor(clusters.shape[0]/2))
	f1 = clusters[:,0]
	p1 = clusters[:,1]
	f2 = clusters[:,3]
	p2 = clusters[:,4]
	phenonew = copy.deepcopy(snpdata.pheno)

	tmp1 = fid[fid.isin(f1[randidx[0:n]])].index.to_numpy()
	tmp2 = pid[snpdata.pid.isin(p1[randidx[0:n]])].index.to_numpy()
	ind = np.intersect1d(tmp1,tmp2)
	phenonew[ind] = 1 - snpdata.pheno[ind]
	tmp1 = fid[fid.isin(f2[randidx[0:n]])].index.to_numpy()
	tmp2 = pid[pid.isin(p2[randidx[0:n]])].index.to_numpy()
	ind = np.intersect1d(tmp1,tmp2)
	phenonew[ind] = 1 - snpdata.pheno[ind]

	return phenonew


