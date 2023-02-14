import numpy as np
import pandas as pd
import pickle
import sys
sys.path.append('classes')
import InteractionNetwork

# Pyhthon script for converting CSSI output to Python interaction networks.
#inputs:
#    cassifile: CASSI output file
#    gitype: genetic interaction type
#            'lin': logistic regression test
#            'lr': logistic regression test 
#            'je': joint effects
#            'awu': adjusted Wu method
#            'afe': adjusted fast epistasis
#            'wz': Wellek Ziegler method
#    outputfile: output file name

if len(sys.argv) < 5:
	sys.exit('not enough parameters')
plinkbimfile = sys.argv[1]
cassifile = sys.argv[2]
gittype = sys.argv[3]
outputfile = sys.argv[4]

# load CASSI interaction file
ssm_cassi = pd.read_csv(cassifile,sep='\s+')

# load bim file
bim = pd.read_csv(plinkbimfile,sep='\s+',engine='python',header=None)
# number of SNPs
p = bim.shape[0]

x = ssm_cassi['SNP1'] - 1
y = ssm_cassi['SNP2'] - 1
ssm_pro = np.zeros((p,p))
ssm_risk = np.zeros((p,p))


# if else, now for lr
if gittype == 'lr':
	# logistic regression (cassi)
	z = -np.log10(ssm_cassi['LR_P'])
	log_or = ssm_cassi['LR_LOG_OR']
	for i in range(x.shape[0]):
		s1 = x[i]
		s2 = y[i]
		if log_or[i] < 0: # protective
			ssm_pro[s1,s2] = z[i]
		if log_or[i] > 0:
			ssm_risk[s1,s2] = z[i]
elif gittype == 'lin':
	# Linear Regression
	z = -np.log10(ssm_cassi['LIN_P'])
	lin_beta = ssm_cassi['LIN_BETA']
	for i in range(x.shape[0]):
		s1 = x[i]
		s2 = y[i]
		if lin_beta[i] < 0: # protective
			ssm_pro[s1,s2] = z[i]
		if lin_beta[i] >= 0:
			ssm_risk[s1,s2] = z[i]
elif gittype == 'je':
	# joint effect (cassi)
	z = -np.log10(ssm_cassi['JE_CC_P'])
	case_log_or = ssm_cassi['JE_CASE_LOG_OR']
	ctrl_log_or = ssm_cassi['JE_CTRL_LOG_OR']
	for i in range(x.shape[0]):
		s1 = x[i]
		s2 = y[i]
		if ctrl_log_or[i] > case_log_or[i]: # protective
			ssm_pro[s1,s2] = z[i]
		if ctrl_log_or[i] <= case_log_or[i]:
			ssm_risk[s1,s2] = z[i]
elif gittype == 'awu':
	# Adjusted Wu
	z = -np.log10(ssm_cassi['AWU_CC_P'])
	case_log_or = ssm_cassi['AWU_CASE_LOG_OR']
	ctrl_log_or = ssm_cassi['AWU_CTRL_LOG_OR']
	for i in range(x.shape[0]):
		s1 = x[i]
		s2 = y[i]
		if ctrl_log_or[i] > case_log_or[i]: # protective
			ssm_pro[s1,s2] = z[i]
		if ctrl_log_or[i] <= case_log_or[i]:
			ssm_risk[s1,s2] = z[i]
elif gittype == 'afe':
	# Adjusted Fast Epistasis
	z = -np.log10(ssm_cassi['AFE_CC_P'])
	case_log_or = ssm_cassi['AFE_CASE_LOG_OR']
	ctrl_log_or = ssm_cassi['AFE_CTRL_LOG_OR']
	for i in range(x.shape[0]):
		s1 = x[i]
		s2 = y[i]
		if ctrl_log_or[i] > case_log_or[i]: # protective
			ssm_pro[s1,s2] = z[i]
		if ctrl_log_or[i] <= case_log_or[i]:
			ssm_risk[s1,s2] = z[i]
elif gittype == 'wz':
	# Wellek Ziegler
	z = -np.log10(ssm_cassi['WZ_CC_P'])
	case_r = ssm_cassi['WZ_CASE_R']
	ctrl_r = ssm_cassi['WZ_CTRL_R']
	for i in range(x.shape[0]):
		s1 = x[i]
		s2 = y[i]
		if ctrl_r[i] > case_r[i]: # protective
			ssm_pro[s1,s2] = z[i]
		if ctrl_r[i] <= case_r[i]:
			ssm_risk[s1,s2] = z[i]


# save to file
ssm_risk = np.maximum(ssm_risk,ssm_risk.T)
ssm_pro = np.maximum(ssm_pro,ssm_pro.T)
network = InteractionNetwork.InteractionNetwork(ssm_risk,ssm_pro,None,None)
final = open(outputfile, 'wb')
pickle.dump(network, final,protocol=4)
final.close()



