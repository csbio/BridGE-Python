import numpy
import math
import pandas as pd
from scipy.io import loadmat
from scipy.stats import hypergeom
from scipy.io import loadmat,savemat
from scipy.stats import chi2_contingency
from scipy.stats import mannwhitneyu


def call_chi2(table):
	tests = table.shape[0]
	results = numpy.zeros(tests)
	for i in range(tests):
		contingency_table = table[i,:]
		obs = numpy.reshape(contingency_table,(2,2))
		g, p, dof, expctd = chi2_contingency(obs)
		results[i] = p
	return results
def call_ranksum(table,ids):
	data2 = table[:,ids]
	data1 = numpy.delete(table,ids,1)
	data2 = numpy.reshape(data2,data2.shape[0]*data2.shape[1])
	data1 = numpy.reshape(data1,data1.shape[0]*data1.shape[1])
	w, p = mannwhitneyu(data1,data2)
	return p

# temporary arguments
input_network = 'network.mat'
path_info = 'bpm.mat'
binary = 0
network_type = 0 ## 0 for positive, 1 for negative
minpath = 5
snp_perms = 10

data = loadmat(input_network)
p_network = data.get('p_network')
n_network = data.get('n_network')

if network_type == 0:
	network = p_network
else:
	network = n_network

d2 = loadmat(path_info)
ind1s = d2.get('ind1s')
ind2s = d2.get('ind2s')
wpmind = d2.get('wpmind')
bpmsize = d2.get('bpmsize')


MM=network ## network


orig_MM = MM
MM = numpy.zeros(orig_MM.shape)
if binary == 0:
	MM[orig_MM > 0.2] = 1


num_of_pairs = ind1s.shape[1]
num_of_wpm = wpmind.shape[1]
snp_size = MM.shape[0]

sumMM = numpy.sum(MM,1)


if binary == 1 or binary==0:
	table1 = numpy.zeros((num_of_pairs,4)) ## f11(bpm interactions) - f10(non-bpm interactions) - f01(bpm non-interations) - f00 (non-bpm non-interactions)
	table2 = numpy.zeros((num_of_pairs,4))
	bpm1size = numpy.zeros(num_of_pairs)
	bpm2size = numpy.zeros(num_of_pairs)
	id1_list = []
	id2_list = []

	## create contingency table for chi-squared local
	for i in range(num_of_pairs):
		id1 = ind1s[0,i]
		id1 = id1 - 1 
		id2 = ind2s[0,i]
		id2 = id2 -  1
		id1 = numpy.reshape(id1,id1.shape[0])
		id2 = numpy.reshape(id2,id2.shape[0])
		id1_list.append(id1)
		id2_list.append(id2)
		## f11
		table1[i,0] = numpy.sum(MM[id1,:][:,id2])
		table2[i,0] = numpy.sum(MM[id1,:][:,id2])
		## f01
		p1_size = id1.shape[0]
		p2_size = id2.shape[0]
		bpm1size[i] = p1_size
		bpm2size[i] = p2_size
		full_possible = p1_size * p2_size
		bpmngi1 = full_possible - table1[i,0]
		bpmngi2 = full_possible - table2[i,0]
		table1[i,2] = bpmngi1
		table2[i,2] = bpmngi2
		## f10
		global_gi1 = numpy.sum(sumMM[id1])
		global_gi2 = numpy.sum(sumMM[id2])
		sub_gi1 = global_gi1 - table1[i,0]
		sub_gi2 = global_gi2 - table2[i,0]
		table1[i,1] = sub_gi1
		table2[i,1] = sub_gi2
		## f00
		background_size1 = p1_size * snp_size
		background_size2 = p2_size * snp_size
		sub_ngi1 = background_size1 - sub_gi1 - table1[i,0]
		sub_ngi2 = background_size2 - sub_gi2 - table2[i,0]
		table1[i,3] = sub_ngi1
		table2[i,3] = sub_ngi2

bpmgi = table1[:,0]
bpmgi = numpy.reshape(bpmgi,num_of_pairs)
bpmnotgi = table1[:,2]
bpmnotgi = numpy.reshape(bpmnotgi,num_of_pairs)
path1bggi = table1[:,1]
path1bggi = numpy.reshape(path1bggi,num_of_pairs)
path1bgnotgi = table1[:,3]
path1bgnotgi = numpy.reshape(path1bgnotgi,num_of_pairs)
path2bggi = table2[:,1]
path2bggi = numpy.reshape(path2bggi,num_of_pairs)
path2bgnotgi = table2[:,3]
path2bgnotgi = numpy.reshape(path2bgnotgi,num_of_pairs)


chi_results1 = call_chi2(table1)
chi_1 = numpy.log10(chi_results1) * -1
chi_1[(bpmgi/(bpmgi+bpmnotgi)) < (path1bggi/(path1bggi+path1bgnotgi)) ] = -1 * chi_1[(bpmgi/(bpmgi+bpmnotgi)) < (path1bggi/(path1bggi+path1bgnotgi)) ]
chi_1 = numpy.reshape(chi_1,(1,num_of_pairs))

chi_results2 = call_chi2(table2)
chi_2 = numpy.log10(chi_results2) * -1
chi_2[(bpmgi/(bpmgi+bpmnotgi)) < (path2bggi/(path2bggi+path2bgnotgi)) ] = -1 * chi_2[(bpmgi/(bpmgi+bpmnotgi)) < (path2bggi/(path2bggi+path2bgnotgi)) ]
chi_2 = numpy.reshape(chi_2,(1,num_of_pairs))

density_bpm_local1 = (bpmgi+path1bggi) / (path1bgnotgi+path1bggi+bpmsize)
density_bpm_local2 = (bpmgi+path2bggi) / (path2bgnotgi+path2bggi+bpmsize)

dense_idx = numpy.zeros((1,num_of_pairs))

# scenario 1: chi2 and density agree
max_ind1 = (density_bpm_local1>density_bpm_local2) & (chi_1<chi_2)
max_ind2 = (density_bpm_local2>density_bpm_local1) & (chi_2<chi_1) 


dense_idx[max_ind1] = 1
dense_idx[max_ind2] = 2

ind_bool = (dense_idx == 0) ## remaining

# scenario 2: only chi2s are equal
max_ind1 = (density_bpm_local1>density_bpm_local2) & (chi_1==chi_2) & ind_bool
max_ind2 = (density_bpm_local2>density_bpm_local1) & (chi_1==chi_2) & ind_bool

dense_idx[max_ind1] = 1
dense_idx[max_ind2] = 2

ind_bool = (dense_idx == 0) ## remaining


# scenario 3: only densities are equal
max_ind1 = (density_bpm_local1==density_bpm_local2) & (chi_1<chi_2) & ind_bool
max_ind2 = (density_bpm_local1==density_bpm_local2) & (chi_2<chi_1) & ind_bool

dense_idx[max_ind1] = 1
dense_idx[max_ind2] = 2

ind_bool = (dense_idx == 0) ## remaining

# scenario 4: chi2 - density don't agree, use chi2 to choose denser

max_ind1 = (chi_1<chi_2) & ind_bool
max_ind2 = (chi_2<chi_1) & ind_bool

dense_idx[max_ind1] = 1
dense_idx[max_ind2] = 2

chi2_local = numpy.zeros((1,num_of_pairs))
chi2_local[dense_idx==1] = chi_1[dense_idx==1]
chi2_local[dense_idx==2] = chi_2[dense_idx==2]

ind2keep = (chi2_local >= numpy.log10(10) ) & (bpm1size >= minpath) & (bpm2size >= minpath)

bpmind1 = ind1s
bpmind2 = ind2s

bpmind1[dense_idx==2] = ind2s[dense_idx==2]
bpmind2[dense_idx==2] = ind2s[dense_idx==2]

## chi2 for wpms
table_wpm = numpy.zeros((num_of_wpm,4))
for i in range(num_of_wpm):
	# f11: number of interactions within pathway
	id1 = wpmind[0,i]
	id1 = id1 - 1 
	id1 = numpy.reshape(id1,id1.shape[0])
	table_wpm[i,0] = numpy.sum(MM[id1,:][:,id1])
	## f01: non interactions within pathway
	p_size = id1.shape[0]
	full_possible = p_size * p_size
	table_wpm[i,2] = full_possible - table_wpm[i,0]
	## f10: interactions out of the pathway
	global_gi = numpy.sum(sumMM[id1])
	global_gi = global_gi - table_wpm[i,0]
	table_wpm[i,1] = global_gi
	## f00: non-interactions out of the pathway
	background_size = snp_size * p_size
	sub_ngi = background_size - table_wpm[i,1] - table1[i,0]
	table_wpm[i,3] = sub_ngi


wpmgi = table_wpm[:,0]
wpmgi = numpy.reshape(wpmgi,num_of_wpm)
wpmnotgi = table_wpm[:,2]
wpmnotgi = numpy.reshape(wpmnotgi,num_of_wpm)
pathbggi = table_wpm[:,1]
pathbggi = numpy.reshape(pathbggi,num_of_wpm)
pathbgnotgi = table_wpm[:,3]
pathbgnotgi = numpy.reshape(pathbgnotgi,num_of_wpm)


chi_results = call_chi2(table_wpm)
chi_wpm = numpy.log10(chi_results) * -1
chi_wpm[(wpmgi/(wpmgi+wpmnotgi)) < (pathbggi/(pathbggi+pathbgnotgi)) ] = -1 * chi_wpm[(wpmgi/(wpmgi+wpmnotgi)) < (pathbggi/(pathbggi+pathbgnotgi)) ]
chi_wpm = numpy.reshape(chi_wpm,(1,num_of_wpm))
ind2keep_wpm = (chi_wpm >= numpy.log10(10))

bpm_local = numpy.zeros((1,num_of_pairs))
bpm_density = numpy.zeros((1,num_of_pairs))
## density of denser
bpmind1_tmp = bpmind1[ind2keep]
bpmind2_tmp = bpmind2[ind2keep]
bpm_density_tmp = bpm_density[ind2keep]
bpm_local_tmp = bpm_local[ind2keep]
num2keep = bpmind1_tmp.shape[0]

for i in range(num2keep):
	id1 = bpmind1_tmp[i]
	id1 = id1 - 1 
	id2 = bpmind2_tmp[i]
	id2 = id2 -  1
	id1 = numpy.reshape(id1,id1.shape[0])
	id2 = numpy.reshape(id2,id2.shape[0])
	total = numpy.sum(orig_MM[id1,:][:,id2])
	s = bpmsize[ind2keep][i]
	bpm_density_tmp[i] = total / s

	## ranksum for non-binary
	if binary == 0:
		m = orig_MM[id1,:]
		p = call_ranksum(m,id2)
		bpm_local_tmp[i] = numpy.log10(p) * -1
	else:
		bpm_local_tmp[i] = chi2_local[ind2keep][i]

bpm_local[ind2keep] = bpm_local_tmp
bpm_density[ind2keep] = bpm_density_tmp

## wpm 
wpm_local = numpy.zeros((1,num_of_wpm))
wpm_density = numpy.zeros((1,num_of_wpm))
wpmind_tmp = wpmind[ind2keep_wpm]
wpm_density_tmp = wpm_density[ind2keep_wpm]
wpm_local_tmp = wpm_local[ind2keep_wpm]
num2keep = wpmind_tmp.shape[0]

for i in range(num2keep):
	id1 = wpmind_tmp[i]
	id1 = id1 - 1
	id1 = numpy.reshape(id1,id1.shape[0])
	total = numpy.sum(orig_MM[id1,:][:,id1]) / 2
	s = id1.shape[0]
	wpm_density_tmp[i] = total/s

	if binary == 0:
		## ranksum
		m = orig_MM[id1,:]
		p = call_ranksum(m,id1)
		wpm_local_tmp[i] = numpy.log10(p) * -1
	else:
		wpm_local_tmp[i] = chi_wpm[ind2keep_wpm][i]

wpm_local[ind2keep_wpm] = wpm_local_tmp
wpm_density[ind2keep_wpm] = wpm_density_tmp

## claculate expected density
expected_bpm_density = numpy.zeros(num_of_pairs)
for i in range(num_of_pairs):
	id1 = bpmind1[0,i]
	id1 = id1 - 1 	
	id1 = numpy.reshape(id1,id1.shape[0])
	p1_size = id1.shape[0]
	total = numpy.sum(orig_MM[id1,:])
	s = p1_size * snp_size
	expected_bpm_density[i] = total / s

expected_wpm_density = numpy.zeros(num_of_wpm)
for i in range(num_of_wpm):
	id1 = wpmind[0,i]
	id1 = id1 - 1
	id1 = numpy.reshape(id1,id1.shape[0])
	p1_size = id1.shape[0]
	total = numpy.sum(orig_MM[id1,:])
	s = p1_size * snp_size
	expected_wpm_density[i] = total / s

## finding path degree
all_snps = numpy.arange(snp_size)
path_degree = numpy.zeros(num_of_wpm)
sumMM = numpy.sum(orig_MM,axis=1)
for i in range(num_of_wpm):
	id1 = wpmind[0,i]
	id1 = id1 - 1
	id1 = numpy.reshape(id1,id1.shape[0])
	data1 = sumMM[id1]
	id_other = numpy.setdiff1d(all_snps,id1)
	data2 = sumMM[id_other]
	w, p = mannwhitneyu(data1,data2)
	path_degree[i] = numpy.log10(p) * -1
ind2keep_path = (path_degree > numpy.log10(10))


## snp permutation to find emprical p-value
bpmsum = table1[:,0]
bpmsum = numpy.reshape(bpmsum,(1,num_of_pairs))
bpm_emp_p = numpy.zeros(num_of_pairs)
numpy.random.seed(66754)
for r in range(snp_perms):
	permuted_idx = numpy.random.permutation(snp_size)
	tempMM = orig_MM[:,permuted_idx]
	## computing sum of bpm interactions in the permutated network
	tmp_sum = numpy.zeros(num_of_pairs)
	for i in range(num_of_pairs):
		id1 = id1_list[i]
		id2 = id2_list[i]
		tmp_sum[i] = numpy.sum(tempMM[id1,:][:,id2])
	bpm_emp_p = bpm_emp_p + (tmp_sum > bpmsum)

print(ind2keep.shape)
print(bpm_emp_p.shape)
bpm_emp_p[ind2keep] = 0
bpm_emp_p = bpm_emp_p / snp_perms

print(bpm_emp_p)












