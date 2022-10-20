import numpy
import math
import pandas as pd
from corefuns import hygetest as ht


# HygeCache class is for computing bulk hypergeometric tests.
# apply_hyge() function is called for computing multipe hypergeometric tests. Sample size, and case/control size is set when creating the object.
# INPUTS:
#	- g: genotype group size(sepcific class we're interested in draw)
#	- x: number of samples with specific genotype in case/control group
#	- case_flag: flag indicating which of the case_size or control_size is used.

class HygeCache:
	def __init__(self,sample_size,case_size,control_size):
		self.sample_size = sample_size
		self.case_size =  case_size
		self.control_size = control_size
		#self.cache = memoize(ht.hygetest)

	def hygetest_caller(self,input_row,case_flag):
		# input_row is : genotype_size,x
		# x is the number of phenotypes which have the genotype
		if case_flag:
			return ht.hygetest(self.sample_size,self.case_size,input_row[1],input_row[0])
		else:
			return ht.hygetest(self.sample_size,self.control_size,input_row[1],input_row[0])

	def call_hygetest(self,input_array,case_flag):
		# input_array is : sample_size * genotype_size,x
		# x is the number of phenotypes which have the genotype
		result = numpy.apply_along_axis(self.hygetest_caller,0,input_array,case_flag)
		return result

	def apply_hyge(self,g,x,case_flag):
		## conacatinate to apply call_hygetest
		y = numpy.stack((g,x))
		raw_out = self.call_hygetest(y,case_flag)
		return raw_out

