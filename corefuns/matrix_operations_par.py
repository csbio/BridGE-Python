import numpy
import math
import pandas as pd
from scipy.io import loadmat
from scipy.stats import hypergeom
from scipy.io import loadmat,savemat
from corefuns import HygeCache as hc
import multiprocessing as mp
import sys
import pickle
from classes import snpsetclass as snps
from classes import InteractionNetwork


# input arguments: model,alpha1,alpha2,plink-cluster file, nworkers,R




class job_quota:
	def __init__(self,population_size,case_size,control_size):
		self.population_size = population_size
		self.case_size = case_size
		self.control_size = control_size
		self.symmetric = True
		self.sx = []
		self.sy = []
		self.i1 = 0
		self.i2 = 0
		self.pheno = []
		self.case_flag = True
		self.alpha1 = 0
		self.alpha2 = 0

class job_result:
	def __init__(self):
		#self.p11 = []
		#self.p10 = []
		#self.p01 = []
		#self.p00 = []
		self.risk = None
		self.protective = None




def matrix_to_array(m,i,symmetric): ## takes a 2D matrix as an input, returns the lower triangle as an 1D array
	x = m.shape[0]
	y = m.shape[1]
	i2 = i + x
	if symmetric:
		idx_low = numpy.tril_indices(x, i-1,y)
		a = m[idx_low]
	else:
		a = numpy.reshape(m,(x*y))
	return a

def array_to_matrix(a,n,m,i,symmetric): ## takes a 1D array as an input(lower triangle), returns the 2D array
	if symmetric:
		idx = numpy.tril_indices(n+i,-1)
		t = i * (i-1) / 2
		t = int(t)
		d0 = idx[0][t:]
		d1 = idx[1][t:]
		idx = (d0,d1)
		matrix = numpy.zeros((m,m))
		matrix[idx] = a
		matrix = matrix[i:i+n,:]
		#matrix = numpy.maximum(matrix,matrix.transpose())
	else:
		matrix = numpy.reshape(a,(n,m))
	return matrix

def parallel_run(job_arg):
	# sub-matrix sx, sub-matrix sy, pheno
	pheno = job_arg.pheno
	sx = job_arg.sx
	sy = job_arg.sy
	i1 = job_arg.i1
	i2 = job_arg.i2
	symmetric_flag = job_arg.symmetric
	pheno_res = numpy.ones(pheno.shape) - pheno
	s = sy.shape[1]


	## reshaping pheno
	pheno = numpy.reshape(pheno,(pheno.shape[0],1))
	pheno_res = numpy.reshape(pheno_res,(pheno_res.shape[0],1))

	# pairwise
	### P11
	#### Risk
	tempx = sx[:] * pheno[:]
	x11 = numpy.matmul(tempx.transpose(),sy)

	#### Protective
	tempx_r = sx[:] * pheno_res[:]
	xp11 = numpy.matmul(tempx_r.transpose(),sy)
	
	### Genotype size for 1-1
	g11 = numpy.matmul(sx.transpose(),sy)
	### 1-snps
	Ix = numpy.ones(sx.shape)
	Iy = numpy.ones(sy.shape)

	sx_res = numpy.subtract(Ix,sx)
	sy_res = numpy.subtract(Iy,sy)



	### P00
	#### Risk
	temp = sx_res[:] * pheno[:]
	x00 = numpy.matmul(temp.transpose(),sy_res)
	#### Protective
	temp_r = sx_res[:] * pheno_res[:]
	xp00 = numpy.matmul(temp_r.transpose(),sy_res)

	### Genotype size for 0-0
	g00 = numpy.matmul(sx_res.transpose(),sy_res)


	### P10
	#### Risk
	x10 = numpy.matmul(tempx.transpose(),sy_res)
	#### Protective
	xp10  = numpy.matmul(tempx_r.transpose(),sy_res)

	### Genotype size for 1-0
	g10 = numpy.matmul(sx.transpose(),sy_res)


	### P01
	#### Risk
	x01 = numpy.matmul(temp.transpose(),sy)
	#### Protective
	xp01 = numpy.matmul(temp_r.transpose(),sy)

	### Genotype size for 0-1
	g01 = numpy.matmul(sx_res.transpose(),sy)

	g11 = matrix_to_array(g11,i1,symmetric_flag)
	g10 = matrix_to_array(g10,i1,symmetric_flag)
	g01 = matrix_to_array(g01,i1,symmetric_flag)
	g00 = matrix_to_array(g00,i1,symmetric_flag)
	x00 = matrix_to_array(x00,i1,symmetric_flag)
	x11 = matrix_to_array(x11,i1,symmetric_flag)
	x01 = matrix_to_array(x01,i1,symmetric_flag)
	x10 = matrix_to_array(x10,i1,symmetric_flag)
	xp00 = matrix_to_array(xp00,i1,symmetric_flag)
	xp11 = matrix_to_array(xp11,i1,symmetric_flag)
	xp01 = matrix_to_array(xp01,i1,symmetric_flag)
	xp10 = matrix_to_array(xp10,i1,symmetric_flag)


	cache = hc.HygeCache(job_arg.population_size,job_arg.case_size,job_arg.control_size)

	## for risk associated
	p11 = cache.apply_hyge(g11,x11,True)
	p01 = cache.apply_hyge(g01,x01,True)
	p10 = cache.apply_hyge(g10,x10,True)
	p00 = cache.apply_hyge(g00,x00,True)
	q = numpy.stack((p01,p10,p00))
	q_min = numpy.amin(q,0)

	raw_out = numpy.divide(p11,q_min)
	log_out = numpy.log10(raw_out)
	log_out = numpy.multiply(log_out,-1)
	# index for those not passing alpha1 filter
	id1 = p11 > job_arg.alpha1

	# index for those not passing alpha2 filter
	id10 = p10 <= job_arg.alpha2
	id01 = p01 <= job_arg.alpha2
	id00 = p00 <= job_arg.alpha2

	idx = id1 | id10 | id01 | id00

	log_out[idx] = 0
	log_out[p11 == 0] = 0
	log_out[q_min == 0] = 0
	log_out[log_out<0] = 0
	## convert 1D result to 2D matrix
	matrix = array_to_matrix(log_out,(i2-i1),s,i1,symmetric_flag)

	result = job_result()
	result.risk = matrix

	## for protective
	p11 = cache.apply_hyge(g11,xp11,False)
	p01 = cache.apply_hyge(g01,xp01,False)
	p10 = cache.apply_hyge(g10,xp10,False)
	p00 = cache.apply_hyge(g00,xp00,False)
	
	q = numpy.stack((p01,p10,p00))
	q_min = numpy.amin(q,0)
	raw_out = numpy.divide(p11,q_min)
	log_out = numpy.log10(raw_out)
	log_out = numpy.multiply(log_out,-1)

	# index for those not passing alpha1 filter
	id1 = p11 > job_arg.alpha1

	# index for those not passing alpha2 filter
	id10 = p10 <= job_arg.alpha2
	id01 = p01 <= job_arg.alpha2
	id00 = p00 <= job_arg.alpha2

	idx = id1 | id10 | id01 | id00

	log_out[idx] = 0
	log_out[p11 == 0] = 0
	log_out[q_min == 0] = 0
	log_out[log_out<0] = 0
	## convert 1D result to 2D matrix
	matrix = array_to_matrix(log_out,(i2-i1),s,i1,symmetric_flag)
	result.protective = matrix

	return result


def run(model,alpha1,alpha2,n_workers,R):

	## output name
	output_name = 'data/ssM_hygessi_'+ model+'_R'+str(R) + '.pkl'
	
	pkl_d = open('data/SNPdataAD.pkl','rb')
	pkl_r = open('data/SNPdataAR.pkl','rb')
	snpdata_d = pickle.load(pkl_d)
	snpdata_r = pickle.load(pkl_r)
	pkl_d.close()
	pkl_r.close()


	pheno = snpdata_r.pheno
	dataR = snpdata_r.data
	dataD = snpdata_d.data



	if model == 'RR':
		datai = dataR
		dataj = dataR
	elif model == 'DD':
		datai = dataD
		dataj = dataD
	else:
		datai = dataR
		dataj = dataD


	population_size = pheno.shape[0]

	## shuffle phenotypes if R != 0
	if R > 0:
		numpy.random.seed(66754)
		for i in range(R):
			permuted_idx = numpy.random.permutation(population_size)
		pheno = pheno[permuted_idx]

	case_size = numpy.count_nonzero(pheno)
	control_size = population_size - case_size
	pheno_res = numpy.ones(pheno.shape) - pheno

	sx = datai
	sy = dataj
	s = sx.shape[1]


	## dividing sx for parallel computing
	idx = [0]
	if model == 'RR' or model == 'DD':
		share = s*s / n_workers
		for i in range(n_workers):
			if i == n_workers -1 :
				idx.append(s)
			else:
				idx.append(math.floor(math.sqrt(idx[i]*idx[i]+ share)))
	else:
		share = math.floor(s / n_workers)
		for i in range(n_workers):
			if i == n_workers - 1:
				idx.append(s)
			else:
				idx.append(idx[i] + share)

	# assign job arguments
	job_args = []
	for i in range(n_workers):
		job_arg = job_quota(population_size,case_size,control_size)
		job_arg.alpha1 = alpha1
		job_arg.alpha2 = alpha2
		if model == 'RR' or model == 'DD':
			job_arg.symmetric = True
		else:
			job_arg.symmetric = False
		job_arg.i1 = idx[i]
		job_arg.i2 = idx[i+1]
		print(sx.shape)
		job_arg.sx = sx.values[:,job_arg.i1:job_arg.i2]
		job_arg.sy = sy.values
		job_arg.pheno = pheno
		job_args.append(job_arg)


	## creating parallel pool
	pool = mp.Pool(processes=n_workers)
	results = pool.map(parallel_run, job_args)

	result_risk = numpy.zeros((s,s))
	result_protective = numpy.zeros((s,s))
	for i in range(n_workers):
		m_r = results[i].risk
		m_p = results[i].protective
		i1 = idx[i]
		i2 = idx[i+1]
		result_risk[i1:i2,:] = m_r
		result_protective[i1:i2,:] = m_p

	## copy over diagonal for DD - RR
	if model == 'RR' or model == 'DD':
		i_upper = numpy.triu_indices(s, 1)
		result_risk[i_upper] = result_risk.T[i_upper]
		result_protective[i_upper] = result_protective.T[i_upper]
	else: # choose maximum, make diagonals 0
		d = numpy.diag_indices(s)
		i_upper = numpy.triu_indices(s, 1)
		#i_down = numpy.tril_indices(s,-1)
		up_risk = result_risk[i_upper]
		down_risk = result_risk.T[i_upper]
		up_protective = result_protective[i_upper]
		down_protective = result_protective.T[i_upper]
		r = numpy.stack((up_risk,down_risk))
		p = numpy.stack((up_protective,down_protective))
		r = numpy.amax(r,0)
		p = numpy.amax(p,0)
		result_risk[i_upper] = r
		result_risk.T[i_upper] = r
		result_risk[d] = 0
		result_protective[i_upper] = p
		result_protective.T[i_upper] = p
		result_protective[d] = 0
	network = InteractionNetwork.InteractionNetwork(result_risk,result_protective,None,None)

	# Save data to pickle file.
	final = open(output_name, 'wb')
	pickle.dump(network, final)
	final.close()


def combine(alpha1,alpha2,n_workers,R):
	output_name = 'data/ssM_hygessi_'+ 'combined_R'+str(R) + '.pkl'
	## run for 3 models
	run('RR',alpha1,alpha2,n_workers,R)
	run('RD',alpha1,alpha2,n_workers,R)
	run('DD',alpha1,alpha2,n_workers,R)

	## load results for 3 models
	rr_file = 'data/ssM_hygessi_'+ 'RR_R'+str(R) + '.pkl'
	rd_file = 'data/ssM_hygessi_'+ 'RD_R'+str(R) + '.pkl'
	dd_file = 'data/ssM_hygessi_'+ 'DD_R'+str(R) + '.pkl'
	pklin = open(rr_file,'rb')
	rr_network = pickle.load(pklin)
	pklin.close()
	pklin = open(rd_file,'rb')
	rd_network = pickle.load(pklin)
	pklin.close()
	pklin = open(dd_file,'rb')
	dd_network = pickle.load(pklin)
	pklin.close()

	s = rr_network.risk.shape[0]

	## compare netowrks, take the element-wise max
	risk_max_id = numpy.zeros((s,s))
	risk_max_id[rr_network.risk > dd_network.risk ] = 1
	risk_max_id[risk_max_id==0] = 2
	risk_max_temp = numpy.maximum(rr_network.risk,dd_network.risk)
	risk_max_id[risk_max_temp<rd_network.risk] = 3
	risk_max = numpy.maximum(risk_max_temp,rd_network.risk)
	risk_max_id[risk_max==0] = 0

	protective_max_id = numpy.zeros((s,s))
	protective_max_id[rr_network.protective > dd_network.protective] = 1
	protective_max_id[protective_max_id==0] = 2
	protective_max_temp = numpy.maximum(rr_network.protective,dd_network.protective)
	protective_max_id[protective_max_temp<rd_network.protective] = 3
	protective_max = numpy.maximum(protective_max_temp,rd_network.protective)
	protective_max_id[protective_max==0] = 0

	## save in the pickle format
	network = InteractionNetwork.InteractionNetwork(risk_max,protective_max,risk_max_id,protective_max_id)
	final = open(output_name, 'wb')
	pickle.dump(network, final)
	final.close()
