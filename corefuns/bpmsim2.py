# fast version

import pandas as pd
import numpy as np

import sys
import datetime

def time_log():
	print(datetime.datetime.now())
	sys.stdout.flush()

def intersection(l1,l2):
	tmp = set(l2)
	l3 = [value for value in l1 if value in tmp]
	return len(l3)


def bpmsim(BPMind1x, BPMind2x, BPMind1y, BPMind2y,proc_num=1):
	nx = BPMind1x.shape[0]
	ny = BPMind1y.shape[0]

	sim = np.zeros((nx,ny))
	test = BPMind1x.equals(BPMind1y)==True and BPMind2x.equals(BPMind2y)==True
	print('test check done...')
	time_log()

	for i in range(nx):
		for j in range(ny):
			print('i:'+str(i))
			print('j:'+str(j))
			time_log()
			if test and j < i:
				continue
			ind1_x = BPMind1x.iloc[0]
			ind2_x = BPMind2x.iloc[0]
			ind1_y = BPMind1y.iloc[0]
			ind2_y = BPMind2y.iloc[0]
			time_log()
			l = min(len(ind1_x)*len(ind2_x),len(ind1_y)*len(ind2_y))
			int1 = intersection(ind1_x,ind1_y)
			int2 = intersection(ind2_x,ind2_y)
			tmp1 = int1 * int2
			int3 = intersection(ind1_x,ind2_y)
			int4 = intersection(ind2_x,ind1_y)
			tmp2 = int3 * int4
			tmp = max(tmp1,tmp2)
			print('tmp computed')
			time_log()
			sim[i,j] = tmp / l

	if test:
		sim = sim + sim.T
	return sim
