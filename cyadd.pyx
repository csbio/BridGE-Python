import numpy as np

cimport numpy as np

DTYPE = np.float64

DTYPE2 = np.int

ctypedef np.float64_t DTYPE_t

ctypedef np.int_t DTYPE2_t

def csum(np.ndarray[DTYPE_t, ndim=2] m,np.ndarray[DTYPE2_t, ndim=1] x,np.ndarray[DTYPE2_t, ndim=1] y ):
	cdef int s = m.shape[0]
	cdef int s1 = x.shape[0]
	cdef int s2 = y.shape[0]
	cdef DTYPE_t sum = 0
	for i in range(s1):
		for j in range(s2):
			sum += m[x[i],y[j]]
	return sum
	

