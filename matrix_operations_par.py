import numpy
import math
import pandas as pd
from scipy.io import loadmat
from scipy.stats import hypergeom
from scipy.io import loadmat,savemat
from hygecache import HygeCache
import multiprocessing as mp
import sys


# input arguments: input_file,model,output_name, nworkers




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
		self.log_out = None




def matrix_to_array(m,i,symmetric): ## takes a 2D matrix as an input, returns the lower triangle as an 1D array
	x = m.shape[0]
	y = m.shape[1]
	if symmetric:
		idx_low = numpy.tril_indices(x, i,y)
		a = m[idx_low]
	else:
		a = numpy.reshape(m,(x*y))
	return a

def array_to_matrix(a,n,m,i,symmetric): ## takes a 1D array as an input(lower triangle), returns the 2D array
	if symmetric:
		idx = numpy.tril_indices(n, i, m)
		matrix = numpy.zeros((n,m))
		matrix[idx] = a
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
	xp10 = p11 = numpy.matmul(tempx_r.transpose(),sy_res)

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


	cache = HygeCache(job_arg.population_size,job_arg.case_size,job_arg.control_size)

	#ps
	p11 = cache.apply_hyge(g11,x11,job_arg.case_flag)
	p01 = cache.apply_hyge(g01,x01,job_arg.case_flag)
	p10 = cache.apply_hyge(g10,x10,job_arg.case_flag)
	p00 = cache.apply_hyge(g00,x00,job_arg.case_flag)
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
	print(i1)
	print(i2)
	matrix = array_to_matrix(log_out,(i2-i1),s,i1,symmetric_flag)

	result = job_result()
	result.log_out = matrix
	return result



#model = 'RD'
#n_workers = 6

data_filename = sys.argv[1]
model = sys.argv[2]
output_name = sys.argv[3]
n_workers = int(sys.argv[4])

## load data
#data_filename = 'SNPdata-py.mat'
data = loadmat(data_filename)
dataR = data.get('dataR')
dataD = data.get('dataD')
pheno = data.get('pheno')

pheno = pheno.astype(int)
dataR = dataR.astype(int)
dataD = dataD.astype(int)


alpha1 = 0.05
alpha2 = 0.05

if model == 'RR':
	datai = dataR
	dataj = dataR
elif model == 'DD':
	datai = dataD
	dataj = dataD
else:
	datai = dataR
	dataj = dataD


population_size = pheno.size
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
	job_arg.sx = sx[:,job_arg.i1:job_arg.i2]
	job_arg.sy = sy
	job_arg.pheno = pheno
	job_args.append(job_arg)


## creating parallel pool
pool = mp.Pool(processes=n_workers)
results = pool.map(parallel_run, job_args)

result = numpy.zeros((s,s))
for i in range(n_workers):
	m = log_out = results[i].log_out
	i1 = idx[i]
	i2 = idx[i+1]
	result[i1:i2,:] = m

## copy over diagonal for DD - RR
if model == 'RR' or model == 'DD':
	i_upper = numpy.triu_indices(s, -1)
	result[i_upper] = result.T[i_upper]




savemat(output_name, {'result': result})












