import numpy as np
import pandas as pd
import sys



if __name__ == '__main__':
	genomefile = sys.argv[1]
	pi_hat = float(sys.argv[2])
#	s = '/'
#	tmp_dir = genomefile.split('/')
#	project_dir = s.join(tmp_dir[0:-1])
	data = pd.read_csv(genomefile,sep='\s+')
	data1 = data.iloc[:,0:2].to_numpy()
	data2 = data.iloc[:,2:4].to_numpy()
	varname = np.concatenate((data1,data2),axis=0)
	varname = tuple(map(tuple,varname))
	varname = list(set(varname))
	s = len(varname)
	n = data.shape[0]

	ibd = np.zeros((s,s))
	for i in range(n):
		v = data['PI_HAT'][i]
		x1 = tuple(data.iloc[i,0:2].to_numpy())
		x2 = tuple(data.iloc[i,2:4].to_numpy())
		id1 = next(j for j,(t1,t2) in enumerate(varname) if (t1,t2) == x1)
		id2 = next(j for j,(t1,t2) in enumerate(varname) if (t1,t2) == x2)
		ibd[id1,id2] = v
	ibd = np.maximum(ibd,ibd.T)
	sum_ibd = np.sum(ibd>=pi_hat,axis=0)
	sum_ibd_sorted = np.sum(ibd>=pi_hat,axis=0)
	sum_ibd_sorted[::-1].sort()
	ind2remove = []
	while sum_ibd_sorted[0] != 0:
		tmp = np.unique(np.argwhere(sum_ibd == sum_ibd_sorted[0]))
		idx = tmp.tolist()
		ibd[idx,:] = 0
		ibd[:,idx] = 0
		sum_ibd = np.sum(ibd>=pi_hat,axis=0)
		sum_ibd_sorted = np.sum(ibd>=pi_hat,axis=0)
		sum_ibd_sorted[::-1].sort()
		ind2remove.extend(idx)
	df = pd.DataFrame(varname, columns=['FID1','IID1'])
	df = df.iloc[ind2remove,:]
	df.to_csv('related_subject2remove.txt',index=False,sep='\t')



