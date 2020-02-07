import numpy as np
import math
import pandas as pd
import pickle
from scipy.io import loadmat
from scipy.stats import hypergeom
from scipy.io import loadmat,savemat
from scipy.stats import chi2_contingency
from scipy.stats import mannwhitneyu
from classes import bpmindclass as bpmc
from classes import InteractionNetwork


def call_chi2(table): ## input format (in each row): ## f11(bpm interactions) - f10(non-bpm interactions) - f01(bpm non-interations) - f00 (non-bpm non-interactions)
	tests = table.shape[0]
	results = numpy.zeros(tests)
	for i in range(tests):
		contingency_table = table[i,:]
		obs = numpy.reshape(contingency_table,(2,2))
		g, p, dof, expctd = chi2_contingency(obs)
		results[i] = p
	return results


def rungenstats(input_network,bpm,wpm,minPath,binary_flag): ## all inputs are classes,loaded from the corresponding pickle file
	## inputs:
	## - input_network: numpy matrix of interaction network
	## - bpm: bpm dataframe
	## - wpm: wpm dataframe
	## - minPath: minimum number of snps in a pathway
	## - binary_flag: flag to make the interaction network binary

	mm = input_network
	s = mm.shape[0]
	bpm_size = bpm['size'][0].shape[1] ## change based on the input format
	bpmsize = np.reshape(bpm['size'][0],bpm_size) ## change based on the input format
	ind1 = bpm['ind1'][0] ## change based on the input format
	ind1 = np.reshape(ind1,bpm_size) ## change based on the input format
	ind2 = bpm['ind2'][0] ## change based on the input format
	ind2 = np.reshape(ind2,bpm_size) ## change based on the input format
	bpmind1size = np.reshape(bpm['ind1size'][0],bpm_size)
	bpmind2size = np.reshape(bpm['ind2size'][0],bpm_size)
	## ?Binary
	if binary_flag:
		# make the network binary
		mm[mm>0.2] = 1
		mm[mm<1] = 0
		sumMM = np.sum(mm,1)
	### BPM binary chi2
		# bpm genetic interaction counts
		bpmgi = np.zeros(bpm_size)
		for i in range(bpm_size):
			id1 = np.reshape(ind1[i],ind1[i].shape[i]) ## change based on the input format
			id2 = np.reshape(ind2[i],ind2[i].shape[i]) ## change based on the input format
			bpmgi[i] = np.sum(mm[id1,:][:,id2])
		# bpm background interactions
		path1bggi = np.zeros(bpm_size)
		path2bggi = np.zeros(bpm_size)
		for i in range(bpm_size):
			id1 = np.reshape(ind1[i],ind1[i].shape[i]) ## change based on the input format
			id2 = np.reshape(ind2[i],ind2[i].shape[i]) ## change based on the input format
			sumgi1 = np.sum(sumMM[id1])
			sumgi2 = np.sum(sumMM[id2])
			path1bggi[i] = sumgi1 - bpmgi[i]
			path2bggi[i] = sumgi2 - bpmgi[i]
		# bpm non interaction
		bpmnotgi = bpmsize - bpmgi
		# non-bpm non-interation
		path1bgsize = bpmind1size * s
		path2bgsize = bpmind2size * s

		path1notgi = path1bgsize - path1bggi - bpmgi
		path2notgi = path2bgsize - path2bggi - bpmgi

		# call chi2
		## build the tables
		table1 = np.stack((bpmgi,path1bggi,bpmnotgi,path1notgi))
		table1 = table1.transpose()
		table2 = np.stack((bpmgi,path2bggi,bpmnotgi,path2notgi))
		table2 = table2.transpose()
		## call chi2
		chi2_bpm_1 = call_chi2(table1)
		chi2_bpm_2 = call_chi2(table2)



	### WPM chi2

	### BPM perm

	### WPM perm

	return


## test parameters for loading the inputs
interaction_input = 'data/ssM_hygessi_RR_R0.pkl'
pklin = open(interaction_input,'rb')
network = pickle.load(pklin)
p_network = network.protective
print(p_network.shape)


#bmpind_input = 'data/BPMind.pkl'
bmpind_input = '../bpmind.pkl'
pklin = open(bmpind_input,'rb')
bb = pickle.load(pklin)
bpm = bb.bpm
wpm = bb.wpm

minPath = 5
binary = True
print(bpm['ind1size'][0][0,1])
#rungenstats(p_network,bpm,wpm,minPath, binary)