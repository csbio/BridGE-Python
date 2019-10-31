import numpy
import math
import pandas as pd
from hygetest import hygetest

class HygeCache:
	def __init__(self,sample_size,case_size,control_size):
		self.sample_size = sample_size
		self.case_size =  case_size
		self.control_size = control_size

	def hygetest_caller(self,input_row,case_flag):
		# input_row is : genotype_size,x
		# x is the number of phenotypes which have the genotype
		if case_flag:
			return hygetest(self.sample_size,self.case_size,input_row[1],input_row[0])
		else:
			return hygetest(self.sample_size,self.control_size,input_row[1],input_row[0])

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

