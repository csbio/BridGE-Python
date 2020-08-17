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
from classes import InteractionNetwork
from classes import GenstatsOut, Stats


def call_chi2(table): ## input format (in each row): ## f11(bpm interactions) - f10(non-bpm interactions) - f01(bpm non-interations) - f00 (non-bpm non-interactions)
	tests = table.shape[0]
	results = np.zeros(tests)
	for i in range(tests):
		contingency_table = table[i,:]
		obs = np.reshape(contingency_table,(2,2))
		g, p, dof, expctd = chi2_contingency(obs)
		results[i] = p
	return results

def ranksum(x,y): ## custom ranksukm similar to the Matlab version
	## x is in the bpm
	## y is out of bpm
	tmp = np.concatenate((y,x))
	rnk = rankdata(tmp)
	# tie correction
	tmp_sorted = np.sort(tmp)
	tmp1 = tmp_sorted[0:-1]
	tmp2 = tmp_sorted[1:]
	tie_index = np.where(tmp1 == tmp2)
	if tie_index[0].shape[0] == 0:
		tie_correction = 0
	else:
		tie_index = np.union1d(tie_index[0],tie_index[0]+1)
		tmp1 = np.asarray(tie_index[0:-1])
		tmp2 = np.asarray(tie_index[1:])
		tie_start = tie_index[1:][(tmp2 - tmp1)>1]
		tie_start = np.insert(tie_start,0,tie_index[0])
		tie_break = np.where(tmp2 - tmp1 > 1)
		tie_break = tie_break[0]
		tie_end = tie_index[tie_break]
		tie_end = np.append(tie_end,tie_index[-1])
		n_tied = tie_end - tie_start
		tie_correction = np.sum(np.multiply((np.multiply(n_tied-1,n_tied)),n_tied+1))/2
	
	nx = x.shape[0]
	ny = y.shape[0]
	ry = np.sum(rnk[0:ny])
	rx = np.sum(rnk[ny:])
	wx = rx
	wmean = nx * (nx + ny + 1) / 2
	tie_var = 2 * tie_correction / ((nx + ny )*(nx + ny - 1))
	wvar = nx * ny * (nx + ny + 1 - tie_var) / 12
	z = (wx - wmean - 0.5)/ np.sqrt(wvar)
	p = norm.cdf(-z)
	return p











def rungenstats(input_network,bpm,wpm,minPath,binary_flag,snpPerms): ## all inputs are classes,loaded from the corresponding pickle file
	## inputs:
	## - input_network: numpy matrix of interaction network
	## - bpm: bpm dataframe
	## - wpm: wpm dataframe
	## - minPath: minimum number of snps in a pathway
	## - binary_flag: flag to make the interaction network binary

	mm = input_network
	s = mm.shape[0]
	bpm_size = bpm['size'].values.shape[0] 
	bpmsize = bpm['size'].values
	ind1 = bpm['ind1'].values
	ind2 = bpm['ind2'].values
	bpmind1size = bpm['ind1size'].values
	bpmind2size = bpm['ind2size'].values

	wpm_size = wpm['size'].values.shape[0]
	wpmsize = wpm['size'].values
	wpmindsize = wpm['indsize'].values
	ind = wpm['ind'].values
	## ?Binary
	if not binary_flag:
		mm_stored = np.copy(mm)
		mm[mm>0.2] = 1
		mm[mm<1] = 0
	sumMM = np.sum(mm,1)
	### BPM binary chi2
	# bpm genetic interaction counts
	bpmgi = np.zeros(bpm_size)
	for i in range(bpm_size):
		id1 = np.array(ind1[i])
		id2 = np.array(ind2[i]) 
		bpmgi[i] = np.sum(mm[id1,:][:,id2])
		# bpm background interactions
	path1bggi = np.zeros(bpm_size)
	path2bggi = np.zeros(bpm_size)
	for i in range(bpm_size):
		id1 = np.array(ind1[i])
		id2 = np.array(ind2[i])
		sumgi1 = np.sum(sumMM[id1])
		sumgi2 = np.sum(sumMM[id2])
		path1bggi[i] = sumgi1 - bpmgi[i]
		path2bggi[i] = sumgi2 - bpmgi[i]
	# bpm non interaction
	bpmnotgi = bpmsize - bpmgi
	# non-bpm non-interation
	path1bgsize = bpmind1size * s
	path2bgsize = bpmind2size * s

	path1notgi = path1bgsize - path1bggi - bpmsize
	path2notgi = path2bgsize - path2bggi - bpmsize

	# call chi2
	## build the tables
	table1 = np.stack((bpmgi,path1bggi,bpmnotgi,path1notgi))
	table1 = table1.transpose()
	table2 = np.stack((bpmgi,path2bggi,bpmnotgi,path2notgi))
	table2 = table2.transpose()
	## call chi2
	chi2_bpm_1 = call_chi2(table1)
	chi2_bpm_2 = call_chi2(table2)
	chi2_bpm_1 = np.log10(chi2_bpm_1) * -1
	chi2_bpm_2 = np.log10(chi2_bpm_2) * -1

	## consider under-enriched chi2s
	chi2_bpm_1[bpmgi/(bpmgi + bpmnotgi) < path1bggi / (path1bggi+path1notgi)] = -1 * chi2_bpm_1[bpmgi/(bpmgi + bpmnotgi) < path1bggi / (path1bggi+path1notgi)]
	chi2_bpm_2[bpmgi/(bpmgi + bpmnotgi) < path2bggi / (path2bggi+path2notgi)] = -1 * chi2_bpm_2[bpmgi/(bpmgi + bpmnotgi) < path2bggi / (path2bggi+path2notgi)]

	## compute densitites
	density_bpm_local_1 = (bpmgi+path1bggi) / (path1notgi+path1bggi + bpmsize)
	density_bpm_local_2 = (bpmgi+path2bggi) / (path2notgi+path2bggi + bpmsize)


	## choose the denser (or lower chi2 value)
	dense_index = np.zeros(bpm_size)
	dense_index[chi2_bpm_1<chi2_bpm_2] = 1
	dense_index[chi2_bpm_1>chi2_bpm_2] = 2
	dense_index[(dense_index==0) & (density_bpm_local_1 > density_bpm_local_2)] = 1
	dense_index[(dense_index==0) & (density_bpm_local_1 < density_bpm_local_2)] = 2

	## finalize bpm local
	chi2_bpm_local = np.zeros(chi2_bpm_1.shape[0])
	chi2_bpm_local[dense_index==1] = chi2_bpm_1[dense_index==1]
	chi2_bpm_local[dense_index==2] = chi2_bpm_2[dense_index==2]

	## keeping track of significant bpms
	ind2keep_bpm = (chi2_bpm_local >= -1 * np.log10(0.05)) & (bpmind1size > minPath) & (bpmind2size > minPath)

		
	## keeping denser pathway in ind1_new 
	ind1_new = np.copy(ind1)
	ind2_new = np.copy(ind2)
	ind1_new[dense_index==1] = ind1[dense_index==1]
	ind1_new[dense_index==2] = ind2[dense_index==2]
	ind2_new[dense_index==1] = ind2[dense_index==1]
	ind2_new[dense_index==2] = ind1[dense_index==2]

	## pairs to keep
	bpmind1 = ind1_new[ind2keep_bpm]
	bpmind2 = ind2_new[ind2keep_bpm]

	###WPM Chi2
	wpmgi = np.zeros(wpm_size)
	for i in range(wpm_size):
		id1 = np.array(ind[i])
		wpmgi[i] = np.sum(mm[id1,:][:,id1])
	wpmnotgi = wpmsize - wpmgi
	density_wpm = wpmgi / wpmsize

	## WPM background size and interactions
	pathbggi = np.zeros(wpm_size)
	for i in range(wpm_size):
		id1 = np.array(ind[i])
		pathbggi[i] = np.sum(sumMM[id1])
	pathbggi = pathbggi - wpmgi
	pathbgsize = wpmindsize * s
	pathbgnotgi = pathbgsize - pathbggi - wpmsize

	wpm_table = np.stack((wpmgi,pathbggi,wpmnotgi,pathbgnotgi))
	wpm_table = wpm_table.transpose()

	## call chi2
	chi2_wpm = call_chi2(wpm_table)
	chi2_wpm = np.log10(chi2_wpm) * -1

	## consider under-enriched chi2s
	chi2_wpm[wpmgi/(wpmgi + wpmnotgi) < pathbggi / (pathbggi+pathbgnotgi)] = -1 * chi2_wpm[wpmgi/(wpmgi + wpmnotgi) < pathbggi / (pathbggi+pathbgnotgi)]
	ind2keep_wpm = (chi2_wpm >= -1 * np.log10(0.05))
	wpmind = ind[ind2keep_wpm]





	##### mutual binary - bon-binary ends here

	if binary_flag:
		## compute bpm interaction count and density for the remaining

		bpmsum = np.zeros(ind1_new.shape[0])
		bpmsum_tmp = np.zeros(bpmind1.shape[0])
		density_bpm = np.zeros(ind1_new.shape[0])

		for i in range(bpmind1.shape[0]):
			id1 = np.array(bpmind1[i])
			id2 = np.array(bpmind2[i])
			bpmsum_tmp[i] = np.sum(mm[id1,:][:,id2])

		density_bpm_tmp = bpmsum_tmp / bpmsize[ind2keep_bpm]
		density_bpm[ind2keep_bpm] = density_bpm_tmp
		bpmsum[ind2keep_bpm] = bpmsum_tmp
		bpm_local = chi2_bpm_local ## output

		### WPM density 
		wpm_local = chi2_wpm
		wpmsum = np.zeros(wpm_size)
		density_wpm = np.zeros(wpm_size)
		wpmsum_tmp = np.zeros(wpmind.shape[0])
		for i in range(wpmind.shape[0]):
			id1 = np.array(wpmind[i])
			wpmsum_tmp[i] = np.sum(mm[id1,:][:,id1])
		density_wpm[ind2keep_wpm] = wpmsum_tmp/wpmsize[ind2keep_wpm]
		wpmsum[ind2keep_wpm] = wpmsum_tmp

	else:
		## restore non-binary mm
		mm = mm_stored
		## ranksum test
		bpm_local_tmp = np.zeros(bpmind1.shape[0])
		bpmsum_tmp = np.zeros(bpmind1.shape[0])
		bpmsum = np.zeros(ind1_new.shape[0])
		density_bpm = np.zeros(ind1_new.shape[0])
		for i in range(bpmind1.shape[0]):
			id1 = np.array(bpmind1[i])
			id2 = np.array(bpmind2[i])
			bpmsum_tmp[i] = np.sum(mm[id1,:][:,id2])
			mmtmp = mm[id1,:]
			mm_in = mmtmp[:,id2]
			mm_out = np.delete(mmtmp,id2,1)
			mm_in = mm_in.flatten()
			mm_out = mm_out.flatten()
			#s,p = wilcoxon(mm_in,y=mm_out,correction=True, alternative='greater')
			p = ranksum(mm_in,mm_out)
			bpm_local_tmp[i] = p
		density_bpm[ind2keep_bpm] = bpmsum_tmp / bpmsize[ind2keep_bpm]
		bpm_local = np.zeros(ind1_new.shape[0])
		bpm_local[ind2keep_bpm] = -1 * np.log10(bpm_local_tmp)
		bpmsum[ind2keep_bpm] = bpmsum_tmp
		## update ind2keep_bpm
		ind2keep_bpm = (bpm_local >= -1 * np.log10(0.05))

		### wpm ranksum
		denisty_wpm = np.zeros(wpm_size)
		wpmsum_tmp = np.zeros(wpmind.shape[0])
		wpmsum = np.zeros(wpm_size)
		wpm_local_tmp = np.zeros(wpmind.shape[0])
		for i in range(wpmind.shape[0]):
			id1 = np.array(wpmind[i])
			mmtmp = mm[id1,:]
			mm_in = mmtmp[:,id1]
			mm_out = np.delete(mmtmp,id1,1)
			mm_in = mm_in.flatten()
			mm_out = mm_out.flatten()
			p = ranksum(mm_in,mm_out)
			wpm_local_tmp[i] = p
			wpmsum_tmp[i] = np.sum(mm[id1,:][:,id1])
		density_wpm[ind2keep_wpm] = wpmsum_tmp/wpmsize[ind2keep_wpm]
		wpm_local = np.zeros(wpm_size)
		wpm_local[ind2keep_wpm] = -1 * np.log10(wpm_local_tmp)
		wpmsum[ind2keep_wpm] = wpmsum_tmp
		ind2keep_wpm = (wpm_local >= -1 * np.log10(0.05)) 






	## compute expected bpm density
	sumMM = np.sum(mm,1)
	bpmind1 = ind1_new[ind2keep_bpm]
	bpmind2 = ind2_new[ind2keep_bpm]
	density_bpm_expected = np.zeros(ind1_new.shape[0])
	for i in range(ind1_new.shape[0]):
		id1 = np.array(ind1_new[i])
		density_bpm_expected[i] = np.sum(sumMM[id1]) / (s * len(ind1_new[i]) )

	## compute expected wpm density
	density_wpm_expected = np.zeros(wpm_size)
	for i in range(wpm_size):
		id1 = np.array(ind[i])
		density_wpm_expected[i] = np.sum(sumMM[id1]) / (s * len(ind[i]) )

	## path degree
	path_degree = np.zeros(wpm_size)
	for i in range(wpm_size):
		id1 = np.array(ind[i])
		dist_in = sumMM[id1]
		dist_out = np.delete(sumMM,id1)
		dist_in = dist_in.flatten()
		dist_out = dist_out.flatten()
		p = ranksum(dist_in,dist_out)
		path_degree[i] =  -1 * np.log10(p)
	ind2keep_path = (path_degree >= -1 * np.log10(0.1) )



	
	## random snp permutation to compute emirical p-value for the significant bpms
	np.random.seed(349898398)
	bpm_local_pv = np.ones(ind1_new.shape[0])
	wpm_local_pv = np.ones(wpm_size)
	path_degree_pv = np.ones(wpm_size)
	arr = np.arange(s)
	wpm_local_pv
	wpmind = ind[ind2keep_wpm]
	pathind = ind[ind2keep_path]
	count_bpm = np.zeros(bpmind1.shape[0])
	count_wpm = np.zeros(wpmind.shape[0])
	count_path = np.zeros(np.sum(ind2keep_path))
	for perm in range(snpPerms):
		np.random.shuffle(arr)
		mmtmp = mm[:,arr]
		sumMMtmp = np.sum(mmtmp,axis=0)
		bpmsum_tmp = np.zeros(bpmind1.shape[0])
		for i in range(bpmind1.shape[0]):
			id1 = np.array(bpmind1[i])
			id2 = np.array(bpmind2[i])
			bpmsum_tmp[i] = np.sum(mmtmp[id1,:][:,id2])
		count_bpm = count_bpm + (bpmsum_tmp > bpmsum[ind2keep_bpm])
		## wpm permutation
		wpmsum_tmp = np.zeros(wpmind.shape[0])
		for i in range(wpmind.shape[0]):
			id1 = np.array(wpmind[i])
			wpmsum_tmp[i] = np.sum(mmtmp[id1,:][:,id1])
		count_wpm = count_wpm + (wpmsum_tmp > wpmsum[ind2keep_wpm])
		## path degree permutation
		path_degree_tmp = np.zeros(pathind.shape[0])
		for i in range(pathind.shape[0]):
			id1 = np.array(pathind[i])
			dist_in = sumMMtmp[id1]
			dist_out = np.delete(sumMMtmp,id1)
			dist_in = dist_in.flatten()
			dist_out = dist_out.flatten()
			p = ranksum(dist_in,dist_out)
			path_degree_tmp[i] =  -1 * np.log10(p)
		count_path = count_path + (path_degree_tmp > path_degree[ind2keep_path])


	bpm_local_pv[ind2keep_bpm] = (count_bpm + 1) / snpPerms
	wpm_local_pv[ind2keep_wpm] = (count_wpm + 1) / snpPerms
	path_degree_pv[ind2keep_path] = (count_path + 1) / snpPerms

	return bpm_local,bpm_local_pv,density_bpm,density_bpm_expected,dense_index,wpm_local,wpm_local_pv,density_wpm,density_wpm_expected,path_degree,path_degree_pv



def genstats(ssmfile,bpmfile,binary_flag,snpPerms,minPath,netDensity=0):
	### load bpmfile
	pklin = open(bpmfile,'rb')
	bpm_obj = pickle.load(pklin)
	bpm = bpm_obj.bpm
	wpm = bpm_obj.wpm
	### load interaction network
	pklin = open(ssmfile,'rb')
	network = pickle.load(pklin)
	p_network = network.protective
	r_network = network.risk
	if binary_flag:
		if netDensity==0:
			p_network[p_network>0] = 1
			r_network[r_network>0] = 1
		else:
			p_cutoff = np.quantile(p_network,1-netDensity)
			r_cutoff = np.quantile(r_network,1-netDensity)
			p_network[p_network<p_cutoff] = 0
			p_network[p_network>=p_cutoff] = 1
			r_network[r_network<r_cutoff] = 0
			r_network[r_network>=r_cutoff] = 1
	p_bpm_local,p_bpm_local_pv,p_density_bpm,p_density_bpm_expected,p_dense_index,p_wpm_local,p_wpm_local_pv,p_density_wpm,p_density_wpm_expected,p_path_degree,p_path_degree_pv = rungenstats(p_network,bpm,wpm,minPath, binary_flag,snpPerms)
	p_stats = Stats.Stats(p_bpm_local,p_bpm_local_pv,p_density_bpm,p_density_bpm_expected,p_dense_index,p_wpm_local,p_wpm_local_pv,p_density_wpm,p_density_wpm_expected,p_path_degree,p_path_degree_pv)
	r_bpm_local,r_bpm_local_pv,r_density_bpm,r_density_bpm_expected,r_dense_index,r_wpm_local,r_wpm_local_pv,r_density_wpm,r_density_wpm_expected,r_path_degree,r_path_degree_pv = rungenstats(r_network,bpm,wpm,minPath, binary_flag,snpPerms)
	r_stats = Stats.Stats(r_bpm_local,r_bpm_local_pv,r_density_bpm,r_density_bpm_expected,r_dense_index,r_wpm_local,r_wpm_local_pv,r_density_wpm,r_density_wpm_expected,r_path_degree,r_path_degree_pv)
	out_obj = GenstatsOut.GenstatsOut(p_stats,r_stats)
	tmp = ssmfile.split('/')
	tmp[-1] = 'genstats_' + tmp[-1]
	outputfile ='/'.join(tmp)
	final = open(outputfile, 'wb')
	pickle.dump(out_obj, final)
	final.close()
	return







