import numpy as np
import math
import pandas as pd
import pickle
#from scipy.io import loadmat
from scipy.stats import hypergeom
#from scipy.io import loadmat,savemat
from scipy.stats import chi2_contingency
from scipy.stats import mannwhitneyu
from scipy.stats import rankdata
from scipy.stats import norm
from classes import bpmindclass as bpmc
from classes import InteractionNetwork
from classes import GenstatsOut, Stats
from corefuns import cyadd
import sys
import gc
import time
import multiprocessing as mp
from multiprocessing import sharedctypes
import os, psutil
import copy
import datetime
sys.path.append('classes')
import InteractionNetwork # for cassi pickle class
np.seterr(divide='ignore', invalid='ignore')

# genstats() computes BPM/WPM/PATH statistics. Can be run parallel.
#
# INPUTS:
#   ssmFile: Interaction networks file in the pickle format.
#   bpmfile: files containing SNP ids for BPM/WPMs in pickle format.
#   binary_flag: If True, interaction scores are binarized for computing BPM/WPM/PATH significances
#   snpPerms: Number of snp permutations used for computing empirical p-values
#   minPath: minimum size for a pathway to be considered as WPM and in BPM.
#	n_workers: number of parallel cpu cores the program shoud use
#
# OUTPUTS:
#   genstats_<ssmFile without extension>.pkl - This pickle file contains a GenstasOut class, which itself contains 2 Stats class oject
#       - protective_stats: Statistics for protective network including ranksum scores,empirical p-values, expected density for BPM/WPMs
#       - risk_stats: Statistics for risk network including ranksum scores,empirical p-values, expected density for BPM/WPMs



class perm_args:
	def __init__(self,id,share,bpmind1,bpmind2,bpmsum,ind2keep_bpm,ind2keep_wpm,ind2keep_path,wpmind,wpmsum,pathind,path_degree,idx1,idx2):
		self.bpmind1 = copy.deepcopy(bpmind1)
		self.bpmind2 = copy.deepcopy(bpmind2)
		self.bpmsum = copy.deepcopy(bpmsum)
		self.ind2keep_bpm = copy.deepcopy(ind2keep_bpm)
		self.ind2keep_wpm = copy.deepcopy(ind2keep_wpm)
		self.ind2keep_path = copy.deepcopy(ind2keep_path)
		self.wpmind = copy.deepcopy(wpmind)
		self.wpmsum = copy.deepcopy(wpmsum)
		self.pathind = copy.deepcopy(pathind)
		self.path_degree = copy.deepcopy(path_degree)
		self.share = share
		self.id = id
		self.idx1 = idx1
		self.idx2 = idx2

class par_rank_args:
	def __init__(self,id,start,end):
		self.id = id
		self.bpmind1 = None
		self.bpmind2 = None
		self.start = start
		self.end = end




def init_worker(mm_in):
	global shared_mm_c
	shared_mm_c = mm_in


def init_worker_perm(mm_in,xs1,xs2):
	global shared_mm_c
	shared_mm_c = mm_in
	
	global shared_xs1
	shared_xs1 = xs1

	global shared_xs2
	shared_xs2 = xs2
		
	
	


def call_chi2(table): ## input format (in each row): ## f11(bpm interactions) - f10(non-bpm interactions) - f01(bpm non-interations) - f00 (non-bpm non-interactions)
	tests = table.shape[0]
	results = np.zeros(tests)
	for i in range(tests):
		contingency_table = table[i,:]
		obs = np.reshape(contingency_table,(2,2))
		if not obs[0,0] == 0:
			g, p, dof, expctd = chi2_contingency(obs,correction=False)
		else:
			p =1
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



def parallel_ranksum(job_arg):
	bpmind1 = job_arg.bpmind1
	bpmind2 = job_arg.bpmind2
	bpm_local_tmp = []
	bpmsum_tmp = []
	tr = []
	mm = np.ctypeslib.as_array(shared_mm_c)
	mm = np.copy(mm)
	for i in range(job_arg.start,job_arg.end):
		id1 = np.array(bpmind1[i])
		id2 = np.array(bpmind2[i])
		if id1.shape[0]< 5 or id2.shape[0] < 5:
			tr.append(i)
			bpmsum_tmp.append(0)
			bpm_local_tmp.append(1)
			continue
		bpmsum_tmp.append(np.sum(mm[id1,:][:,id2]))
		mmtmp = mm[id1,:]
		mm_in = mmtmp[:,id2]
		mm_out = np.delete(mmtmp,id2,1)
		mm_in = mm_in.flatten()
		mm_out = mm_out.flatten()
		#s,p = wilcoxon(mm_in,y=mm_out,correction=True, alternative='greater')
## This section is for tracking memory. uncomment from if to sys.stdout.flush() for getting memory usage info
#		if i % 10000 == 0:
#			print('calling ranksum:'+str(i))
#			t = time.localtime()
#			current_time = time.strftime("%H:%M:%S", t)
#			print(current_time)
#			print('Memory:')
#			process = psutil.Process(os.getpid())
#			print(process.memory_info().rss)
#			sys.stdout.flush()
		p2 = mannwhitneyu(mm_in,mm_out,use_continuity=True,alternative='greater')
		p = p2[1]
		bpm_local_tmp.append(p)
	return bpmsum_tmp,bpm_local_tmp,tr


def snp_permutation_parallel(perm_args):
	# retrieve local arrays
	bpmind1 = perm_args.bpmind1
	bpmind2 = perm_args.bpmind2
	bpmsum = perm_args.bpmsum
	ind2keep_bpm = perm_args.ind2keep_bpm
	ind2keep_wpm = perm_args.ind2keep_wpm
	ind2keep_path = perm_args.ind2keep_path
	wpmind = perm_args.wpmind
	wpmsum = perm_args.wpmsum
	pathind = perm_args.pathind
	path_degree = perm_args.path_degree
	share = perm_args.share
	idx1 = perm_args.idx1
	idx2 = perm_args.idx2

	mmtmp_c = np.ctypeslib.as_array(shared_mm_c)
	mmtmp = np.copy(mmtmp_c)

	xs1 = np.ctypeslib.as_array(shared_xs1)
	xs2 = np.ctypeslib.as_array(shared_xs2)

	# set seed based on id
	np.random.seed(349898398+perm_args.id*1000)

	# permutation
	count_bpm = np.zeros(bpmind1.shape[0])
	count_wpm = np.zeros(wpmind.shape[0])
	count_path = np.zeros(np.sum(ind2keep_path))


	gc.collect()
	for perm in range(share):
		mmtmp = mmtmp[:,np.random.permutation(mmtmp.shape[1])]
		sumMMtmp = np.sum(mmtmp,axis=0)
		bpmsum_tmp = np.zeros(bpmind1.shape[0])
		d1 = 0
		d2 = 0
		for i in range(bpmind1.shape[0]):
			t1 = datetime.datetime.now()
			id1 = np.array(bpmind1[i])
			id2 = np.array(bpmind2[i])
			p1 = idx1[i]
			p2 = idx2[i]
			t2 = datetime.datetime.now()
			dt = t2 - t1
			d1 = d1 + dt.total_seconds() * 1000
			if id1.shape[0] == 0 or id2.shape[0] == 0:
				bpmsum_tmp[i] = 0
			else:
				t1 = datetime.datetime.now()
				bpmsum_tmp[i] = cyadd.csum(mmtmp,id1,id2)
				t2 = datetime.datetime.now()
				dt = t2 - t1
				d2 = d2 + dt.total_seconds() * 1000
		count_bpm = count_bpm + (bpmsum_tmp > bpmsum[ind2keep_bpm])
		## wpm permutation
		wpmsum_tmp = np.zeros(wpmind.shape[0])
		for i in range(wpmind.shape[0]):
			id1 = np.array(wpmind[i])
			#wpmsum_tmp[i] = np.sum(mmtmp[id1,:][:,id1])
			wpmsum_tmp[i] = cyadd.csum(mmtmp,id1,id1)
		count_wpm = count_wpm + (wpmsum_tmp > wpmsum[ind2keep_wpm])
		## path degree permutation
		path_degree_tmp = np.zeros(pathind.shape[0])
		for i in range(pathind.shape[0]):
			id1 = np.array(pathind[i])
			dist_in = sumMMtmp[id1]
			dist_out = np.delete(sumMMtmp,id1)
			dist_in = dist_in.flatten()
			dist_out = dist_out.flatten()
			p2 = mannwhitneyu(dist_in,dist_out,use_continuity=True,alternative='greater')
			p = p2[1] 
			path_degree_tmp[i] =  -1 * np.log10(p)
		count_path = count_path + (path_degree_tmp > path_degree[ind2keep_path])
	return count_bpm,count_wpm,count_path








def rungenstats(input_network,bpm,wpm,minPath,binary_flag,snpPerms,n_workers): ## all inputs are classes,loaded from the corresponding pickle file
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
		mm[mm>=0.2] = 1
		mm[mm<1] = 0
	sumMM = np.sum(mm,1)
	### BPM binary chi2
	# bpm genetic interaction counts
	tr = []
	bpmgi = np.zeros(bpm_size)
	for i in range(bpm_size):
		id1 = np.array(ind1[i])
		id2 = np.array(ind2[i])
		if id1.shape[0] < 5 or id2.shape[0] < 5:
			tr.append(i)
			continue
		bpmgi[i] = np.sum(mm[id1,:][:,id2])
	# bpm background interactions
	path1bggi = np.zeros(bpm_size)
	path2bggi = np.zeros(bpm_size)
	for i in range(bpm_size):
		if i in tr:
			continue
		id1 = np.array(ind1[i])
		id2 = np.array(ind2[i])
		if id1.shape[0] == 0:
			sumgi1 = 1
		else:
			sumgi1 = np.sum(sumMM[id1])
		if id2.shape[0] == 0:
			sumgi2 = 1
		else:
			sumgi2 = np.sum(sumMM[id2])
		path1bggi[i] = sumgi1 - bpmgi[i]
		path2bggi[i] = sumgi2 - bpmgi[i]
	# bpm non interaction
	bpmnotgi = bpmsize - bpmgi
	bpmnotgi[bpmnotgi<0] = 0
	# non-bpm non-interation
	path1bgsize = bpmind1size * s
	path2bgsize = bpmind2size * s

	path1notgi = path1bgsize - path1bggi - bpmsize
	path2notgi = path2bgsize - path2bggi - bpmsize

	# call chi2
	## build the tables
	table1 = np.stack((bpmgi,path1bggi,bpmnotgi,path1notgi))
	table1 = table1.transpose()
	table1[tr,:] = 5
	sys.stdout.flush()
	table2 = np.stack((bpmgi,path2bggi,bpmnotgi,path2notgi))
	table2 = table2.transpose()
	table2[tr,:] = 5
	## call chi2
	chi2_bpm_1 = call_chi2(table1)
	chi2_bpm_2 = call_chi2(table2)
	chi2_bpm_1 = np.log10(chi2_bpm_1) * -1
	chi2_bpm_2 = np.log10(chi2_bpm_2) * -1
	chi2_bpm_1[tr] = 0
	chi2_bpm_2[tr] = 0

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
	ind2keep_bpm = (chi2_bpm_local >= -1 * np.log10(0.1)) & (bpmind1size >= minPath) & (bpmind2size >= minPath)

		
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
	ind2keep_wpm = (chi2_wpm >= -1 * np.log10(0.1))
	wpmind = ind[ind2keep_wpm]
	gc.collect()
	



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
		bpmsum = np.zeros(ind1_new.shape[0])
		density_bpm = np.zeros(ind1_new.shape[0])
		# parallel run for computing ranksum
		## compute share for each process
		share = math.floor(bpmind1.shape[0]/n_workers)
		job_args = []
		## init args
		for i in range(n_workers):
			start = i * share
			end = (i+1) * share
			if i == n_workers-1:
				end = bpmind1.shape[0]
			j = par_rank_args(i,start,end)
			j.bpmind1 = copy.deepcopy(bpmind1)
			j.bpmind2 = copy.deepcopy(bpmind2)
			job_args.append(j)
		## share mm
		mmtmp  = np.ctypeslib.as_ctypes(mm)
		shared_mm = sharedctypes.RawArray(mmtmp._type_, mmtmp)
		## call processes
		pool = mp.Pool(processes=n_workers,initializer=init_worker,initargs=(shared_mm,))
		results = pool.map(parallel_ranksum, job_args)
		## retrieve the results
		bpm_tmp = []
		bpmsum_tmp = []
		tr = []
		for i in range(n_workers):
			bpmsum_tmp = bpmsum_tmp + results[i][0]
			bpm_tmp = bpm_tmp + results[i][1]
			tr = tr + results[i][2]
		bpmsum_tmp = np.array(bpmsum_tmp)
		bpm_local_tmp = np.array(bpm_tmp)

		
		bpm_local_tmp[tr] = 1
		bpmsum_tmp[tr] = 0
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
			p2 = mannwhitneyu(mm_in,mm_out,use_continuity=True,alternative='greater')
			p = p2[1] 
			wpm_local_tmp[i] = p
			wpmsum_tmp[i] = np.sum(mm[id1,:][:,id1])
		density_wpm[ind2keep_wpm] = wpmsum_tmp/wpmsize[ind2keep_wpm]
		wpm_local = np.zeros(wpm_size)
		wpm_local[ind2keep_wpm] = -1 * np.log10(wpm_local_tmp)
		wpmsum[ind2keep_wpm] = wpmsum_tmp
		ind2keep_wpm = (wpm_local >= -1 * np.log10(0.05)) 

	gc.collect()
	## compute expected bpm density
	sumMM = np.sum(mm,1)
	bpmind1 = ind1_new[ind2keep_bpm]
	bpmind2 = ind2_new[ind2keep_bpm]
	density_bpm_expected = np.zeros(ind1_new.shape[0])
	for i in range(ind1_new.shape[0]):
		id1 = np.array(ind1_new[i])
		if id1.shape[0] == 0:
			continue 
		try:
			density_bpm_expected[i] = np.sum(sumMM[id1]) / (s * len(ind1_new[i]) )
		except Exception as e:
			print(i)
			print(id1.shape)
			print(id1)

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
		p2 = mannwhitneyu(dist_in,dist_out,use_continuity=True,alternative='greater')
		p = p2[1] 
		path_degree[i] =  -1 * np.log10(p)
	ind2keep_path = (path_degree >= -1 * np.log10(0.1) )



	gc.collect()
	## random snp permutation to compute emirical p-value for the significant bpms
	np.random.seed(349898398)
	bpm_local_pv = np.ones(ind1_new.shape[0])
	wpm_local_pv = np.ones(wpm_size)
	path_degree_pv = np.ones(wpm_size)
	arr = np.arange(s)
	wpm_local_pv
	wpmind = ind[ind2keep_wpm]
	pathind = ind[ind2keep_path]

	# for speeding up bpm sum
	p_range = np.arange(0,wpm_size,1)
	xx,yy = np.meshgrid(p_range,p_range)
	i = np.triu_indices(wpm_size,1)
	idx1 = yy[i]
	idx2 = xx[i]
	idx1 = idx1[ind2keep_bpm]
	idx2 = idx2[ind2keep_bpm]

	idx1_unique = np.unique(idx1)
	idx2_unique = np.unique(idx2)

	pre_comp2 = np.zeros((s,wpm_size))
	for ptw in idx2_unique:
		loc = np.where(idx2 == ptw)
		loc = loc[0][0]
		id2 = np.array(bpmind2[loc])
		pre_comp2[id2,ptw] = 1

	pre_comp1 = np.zeros((s,wpm_size))
	for ptw in idx1_unique:
		loc = np.where(idx1 == ptw)
		loc = loc[0][0]
		id1 = np.array(bpmind1[loc])
		pre_comp1[id1,ptw] = 1

	xs1 = np.ctypeslib.as_ctypes(pre_comp1)
	xs1 = sharedctypes.RawArray(xs1._type_, xs1)

	xs2 = np.ctypeslib.as_ctypes(pre_comp1)
	xs2 = sharedctypes.RawArray(xs2._type_, xs2)





	count_bpm = np.zeros(bpmind1.shape[0])
	count_wpm = np.zeros(wpmind.shape[0])
	count_path = np.zeros(np.sum(ind2keep_path))
	mmtmp = mm
	mmtmp  = np.ctypeslib.as_ctypes(mm)
	shared_mm = sharedctypes.RawArray(mmtmp._type_, mmtmp)
	# assign paralel job args
	share = math.floor(snpPerms/n_workers)
	last_share = snpPerms - (n_workers-1)*share
	job_args = []
	for proc in range(n_workers):
		if proc == 0:
			proc_share = last_share
		else:
			proc_share = share
		j_arg = perm_args(proc,proc_share,bpmind1,bpmind2,bpmsum,ind2keep_bpm,ind2keep_wpm,ind2keep_path,wpmind,wpmsum,pathind,path_degree,idx1,idx2)
		j_arg.bpmind1 = copy.deepcopy(bpmind1)
		j_arg.bpmind2 =	copy.deepcopy(bpmind2)
		j_arg.wpmind =	copy.deepcopy(wpmind)	
		job_args.append(j_arg)
		

	pool = mp.Pool(processes=n_workers,initializer=init_worker_perm,initargs=(shared_mm,xs1,xs2))
	results = pool.map(snp_permutation_parallel, job_args)
	# combine results
	for proc in range(n_workers):
		count_bpm = count_bpm + results[proc][0]
		count_wpm = count_wpm + results[proc][1]
		count_path = count_path + results[proc][2]



	bpm_local_pv[ind2keep_bpm] = (count_bpm + 1) / snpPerms
	wpm_local_pv[ind2keep_wpm] = (count_wpm + 1) / snpPerms
	path_degree_pv[ind2keep_path] = (count_path + 1) / snpPerms

	return bpm_local,bpm_local_pv,density_bpm,density_bpm_expected,dense_index,wpm_local,wpm_local_pv,density_wpm,density_wpm_expected,path_degree,path_degree_pv



def genstats(ssmfile,bpmfile,binary_flag,snpPerms,minPath,n_workers,netDensity=None):
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
		if netDensity == None:
			p_network[p_network>0] = 1
			r_network[r_network>0] = 1
		else:
			p_cutoff = np.quantile(p_network,1-netDensity)
			r_cutoff = np.quantile(r_network,1-netDensity)
			p_network[p_network<p_cutoff] = 0
			p_network[p_network>=p_cutoff] = 1
			r_network[r_network<r_cutoff] = 0
			r_network[r_network>=r_cutoff] = 1
	p_bpm_local,p_bpm_local_pv,p_density_bpm,p_density_bpm_expected,p_dense_index,p_wpm_local,p_wpm_local_pv,p_density_wpm,p_density_wpm_expected,p_path_degree,p_path_degree_pv = rungenstats(p_network,bpm,wpm,minPath, binary_flag,snpPerms,n_workers)
	p_stats = Stats.Stats(p_bpm_local,p_bpm_local_pv,p_density_bpm,p_density_bpm_expected,p_dense_index,p_wpm_local,p_wpm_local_pv,p_density_wpm,p_density_wpm_expected,p_path_degree,p_path_degree_pv)
	r_bpm_local,r_bpm_local_pv,r_density_bpm,r_density_bpm_expected,r_dense_index,r_wpm_local,r_wpm_local_pv,r_density_wpm,r_density_wpm_expected,r_path_degree,r_path_degree_pv = rungenstats(r_network,bpm,wpm,minPath, binary_flag,snpPerms,n_workers)
	gc.collect()	
	r_stats = Stats.Stats(r_bpm_local,r_bpm_local_pv,r_density_bpm,r_density_bpm_expected,r_dense_index,r_wpm_local,r_wpm_local_pv,r_density_wpm,r_density_wpm_expected,r_path_degree,r_path_degree_pv)
	out_obj = GenstatsOut.GenstatsOut(p_stats,r_stats)
	tmp = ssmfile.split('/')
	tmp[-1] = 'genstats_' + tmp[-1]
	outputfile ='/'.join(tmp)
	final = open(outputfile, 'wb')
	pickle.dump(out_obj, final)
	final.close()
	return







