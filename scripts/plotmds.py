import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys



# This script is used for plotting MDS plot on two dimmension for study population and multiple 1000 genome project populations.
#
# INPUTS:
#   mds_file: MDS data gemerated by PLINK
#	pop_id_file: population data(FID,IID) combined with one-hot encoded data of the sub-population membership. This is generated
#		by data_checkpopulation.sh script
#	prj1000Pop: Additional population from 1000 genome project to be included along with CEU, CHB, JPI and ASW.
#	outfile: output file name. pdf extenstion will be added
#
# OUTPUTS:
#   A outfile.pdf file showing the scatter plot on the top two MDS dimensions. Different populations are colored differently.
# 

# Example inputs
#mds_file ='gwas_subset_prj1000.mds'
#pop_id_file = 'allpopid.txt'
#prj1000Pop = 'GIH'
#outfile = 'MDS_gwas_subset'


if __name__ == '__main__':
	mds_file = sys.argv[1]
	pop_id_file = sys.argv[2]
	prj1000Pop = sys.argv[3]
	outfile = sys.argv[4]

	mds_data = pd.read_csv(mds_file,index_col=False,sep='\s+')
	pop_id_data = pd.read_csv(pop_id_file,header=None,skiprows=1,sep='\s+')
	tmp = pd.read_csv(pop_id_file,sep='\s+')
	tmp_header = tmp.columns.tolist()
	tmp_header.insert(0,'IID')
	tmp_header.insert(0,'FID')
	pop_id_data.columns = tmp_header

	mds_data = mds_data.merge(pop_id_data,how='inner')
	mds_data['population'] = (mds_data.iloc[:, 5:] == 1).idxmax(1)
	mds_data.sort_values('population')

	sns.scatterplot(data=mds_data,x='C1', y='C2', hue='population',marker="+",s=10)
	plt.title("MDS plot colored by population")
	plt.xlabel('C1')
	plt.ylabel('C2')
	plt.grid()
	plt.savefig(outfile+'.pdf')
	plt.close()

