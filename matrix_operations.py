import numpy
import math
import pandas as pd
from scipy.io import loadmat
from scipy.stats import hypergeom
from scipy.io import loadmat,savemat

model = 'RD'

## load data
data_filename = 'SNPdata-py.mat'
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


## individuals

### P1~
#### Risk
sx = datai
p1_ = numpy.matmul(sx.transpose(),pheno) 

#### Protective
pp1_ = numpy.matmul(sx.transpose(),pheno_res)


### P~1
#### Risk
sy = dataj
p_1 = numpy.matmul(sy.transpose(),pheno)


#### Protective
pp_1 = numpy.matmul(sy.transpose(),pheno_res)

## pairwise

### P11
#### Risk
tempx = sx[:] * pheno[:]
p11 = numpy.matmul(tempx.transpose(),sy)

#### Protective
tempx_r = sx[:] * pheno_res[:]
p11 = numpy.matmul(tempx_r.transpose(),sy)

### 1-snps
I = numpy.ones(sx.shape)
sx_res = numpy.subtract(I,sx)
sy_res = numpy.subtract(I,sy)

### P00
#### Risk
temp = sx_res[:] * pheno[:]
p00 = numpy.matmul(temp.transpose(),sy_res)
#### Protective
temp_r = sx_res[:] * pheno_res[:]
p00 = numpy.matmul(temp_r.transpose(),sy_res)

### P10
#### Risk
p10 = numpy.matmul(tempx.transpose(),sy_res)
#### Protective
pp10 = p11 = numpy.matmul(tempx_r.transpose(),sy_res)

### P01
#### Risk
p00 = numpy.matmul(temp.transpose(),sy)
#### Protective
p00 = numpy.matmul(temp_r.transpose(),sy)






