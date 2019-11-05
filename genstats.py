import numpy
import math
import pandas as pd
from scipy.io import loadmat
from scipy.stats import hypergeom
from scipy.io import loadmat,savemat
from scipy.stats import chi2_contingency


def call_chi2(table):
	tests = table.shape[0]
	results = numpy.zeros(tests)
	for i in range(tests):
		contingency_table = tests[i,:]
		obs = numpy.reshape(contingency_table,(2,2))
		g, p, dof, expctd = chi2_contingency(obs)
		results[i] = p
	return results

# temporary arguments
input_network = 'network.mat'
path_info = 'bpm.mat'
binary = 1

data = loadmat(input_network)
p_network = data.get('p_network')
p_network = data.get('n_network')

d2 = loadmat(path_info)
ind1s = d2.get('ind1s')
ind2s = d2.get('ind2s')


MM=[] ## network
BPMind1=[] ## indices of pathways
BPMind2=[]

num_of_pairs = ind1s.shape[1]
snp_size = MM.shape[0]

sumMM = numpy.sum(MM,1)


if binary == 1:
	table1 = numpy.zeros(num_of_pairs,4) ## f11(bpm interactions) - f10(non-bpm interactions) - f01(bpm non-interations) - f00 (non-bpm non-interactions)
	table2 = numpy.zeros(num_of_pairs,4)
	bpmgi1 = numpy.zeros(num_of_pairs) ## bpm interactions
	bpmgi2 = numpy.zeros(num_of_pairs)
	## create contingency table for chi-squared local
	for i in range(num_of_pairs):
		id1 = ind1s[0,i]
		id2 = ind1s[0,i]
		## f11
		cid = numpy.intersect1d(id1,id2)
		bpmi1[i] = numpy.sum(sumMM[cid])
		bpmi2[i] = numpy.sum(sumMM[cid])
		table1[i,0] = bpmi1[i]
		table2[i,0] = bpmi2[i]
		## f01
		p1_size = id1.shape[0]
		p2_size = id2.shape[0]
		full_possible = p1_size * p2_size
		bpmngi1 = full_possible - bpmi1[i]
		bpmngi2 = full_possible - bpmi2[i]
		table1[i,2] = bpmngi1
		table2[i,2] = bpmngi2
		## f10
		global_gi1 = numpy.sum(sumMM[id1])
		global_gi2 = numpy.sum(sumMM[id2])
		sub_gi1 = global_gi1 - bpmgi1[i]
		sub_gi2 = global_gi2 - bpmgi2[i]
		table1[i,1] = sub_gi1
		table2[i,1] = sub_gi2
		## f00
		background_size1 = p1_size * snp_size
		background_size2 = p2_size * snp_size
		sub_ngi1 = background_size1 - sub_gi1 - bpmi1[i]
		sub_ngi2 = background_size2 - sub_gi2 - bpmi2[i]
		table1[i,3] = sub_ngi1
		table2[i,3] = sub_ngi2

chi_results1 = call_chi2(table1)
chi_results2 = call_chi2(table2)









