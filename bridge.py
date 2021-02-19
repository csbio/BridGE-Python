import sys
from os import path
from datatools import plink2pkl as p2p
from datatools import bindataa as ba
from datatools import msigdb2pkl as msig2p
from datatools import mapsnp2gene as snp2gene
from datatools import snppathway as snpp
from datatools import bpmind as bpm
from corefuns import matrix_operations_par as ci
from corefuns import genstats as gs
from corefuns import fdrsampleperm as fdr
from corefuns import collectresults as clr

## main caller of bridge

## inputs:
# --job : desired job
# --plinkFile: plink file
# --genesets: geneset symbolsFile and entrezFile without extension

if __name__ == '__main__':
	## parse input arguments
	job = ''
	plinkfile = ''
	genesets='data/toy.genesets' ## default - #change
	gene_annotation = 'data/glist-hg19' ## default - #change
	mappingDistance = 50000
	minPath = 10
	maxPath = 300
	alpha1 = 0.05
	alpha2 = 0.05
	n_workers = 4
	sample_perms = 1
	binaryNetwork = False
	snpPerms = 100
	pval_cutoff = 0.05
	fdrcut = 0.4
	i = -1
	for arg in sys.argv:
		if '=' in arg and '--' in arg:
			o = arg.split('=')[0]
			a = arg.split('=')[1]
			if o == '--job':
				job = a
			elif o == '--plinkFile':
				plinkfile = a
			elif o == '--genesets':
				genesets = a
			elif o == '--geneAnnotation':
				gene_annotation = a
			elif o == '--mappingDistance':
				mappingDistance = int(a)
			elif o == '--maxPath':
				maxPath = int(a)
			elif o == '--minPath':
				minPath = int(a)
			elif o == '--model':
				model = a
			elif o == '--nWorker':
				n_workers = int(a)
			elif o == '--samplePerms':
				sample_perms = int(a)
			elif o == '--binaryNetwork':
				if int(a) == 1:
					binaryNetwork = True
			elif o == '--snpPerms':
				snpPerms = int(a)
			elif o == '--pvalueCutoff':
				pval_cutoff = float(a)
			elif o == '--fdrCutoff':
				fdrcut = flot(a)
			elif o == '--i':
				i = int(a)


	## run job
	if job == 'DataProcess':
		## data processing step

		## convert plinkfile to pickle
		if plinkfile == '':
			sys.exit('plinkFile not provided')
		rawfile = plinkfile + '.raw'
		bimfile = plinkfile + '.bim'
		famfile = plinkfile + '.fam'
		if not path.exists(rawfile) or not path.exists(bimfile) or not path.exists(famfile):
			sys.exit('plinkFiles do not exist')
		finalfile = plinkfile + '.pkl'
		p2p.plink2pkl(rawfile, bimfile, famfile, finalfile)
		## converting snp data based on the disease model
		ba.bindataa(finalfile,'r')
		ba.bindataa(finalfile,'d')
		symbolsFile = genesets + '.symbols.gmt'
		entrezFile = genesets + '.entrez.gmt'
		if not path.exists(symbolsFile) or not path.exists(entrezFile):
			sys.exit('genesets do not exist')
		## prepare gene set information
		msig2p.msigdb2pkl(symbolsFile, entrezFile)
		if not path.exists(gene_annotation):
			sys.exit('gene annotation file not found')
		sgmFile = 'snpgenemapping_' + str(int(mappingDistance/1000)) + 'kb.pkl'
		s = '/'
		dir_sgm = bimfile.split('/')
		dir_sgm[-1] = sgmFile
		sgmfile = s.join(dir_sgm)
		## build relationship between snps and genes
		snp2gene.mapsnp2gene(bimfile, gene_annotation, mappingDistance, 'matrix', sgmfile) # matrix mode - #change
		## extract snp-pathway information
		geneset_pkl = genesets + '.pkl'
		outfile = snpp.snppathway(finalfile, sgmfile, geneset_pkl, minPath, maxPath)
		bpm.bpmind(outfile)
	elif job == 'ComputeInteraction':
		## compute interaction network here
		## input parameters
		if not (model == 'RR' or model == 'RD' or model == 'DD' or model == 'combined'):
			sys.exit('wrong model')
		if not path.exists('data/SNPdataAD.pkl'):
			sys.exit('data/SNPdataAD.pkl not found')
		if not path.exists('data/SNPdataAR.pkl'):
			sys.exit('data/SNPdataAR.pkl not found')
		if model == 'combined':
			if i == -1:
				for R in range(sample_perms+1):
					ci.combine(alpha1,alpha2,n_workers,R)
			else:
				ci.combine(alpha1,alpha2,n_workers,i)
		else:
			if i == -1:
				for R in range(sample_perms+1):
					ci.run(model,alpha1,alpha2,n_workers,R)
			else:
				ci.run(model,alpha1,alpha2,n_workers,i)

	elif job == 'SamplePermutation':
		if not (model == 'RR' or model == 'RD' or model == 'DD' or model == 'combined'):
			sys.exit('wrong model')
		bpmfile = 'data/BPMind.pkl'
		if not path.exists(bpmfile):
			sys.exit('data/BPMind.pkl not found')
		if not path.exists('data/SNPdataAD.pkl'):
                	sys.exit('data/SNPdataAD.pkl not found')
		if not path.exists('data/SNPdataAR.pkl'):
                	sys.exit('data/SNPdataAR.pkl not found')
		for R in range(sample_perms+1):
			if model == 'combined':
				ci.combine(alpha1,alpha2,n_workers,R)
				# build ssM file name
				ssmfile = 'data/ssM_hygessi_combined_R'+ str(R) + '.pkl'
			else:
				ci.run(model,alpha1,alpha2,n_workers,R)	
				ssmfile = 'data/ssM_hygessi_' + model + '_R'+ str(R) + '.pkl'
			gs.genstats(ssmfile,bpmfile,binaryNetwork,snpPerms,minPath)

	elif job == 'Analysis':
		bpmfile = 'data/BPMind.pkl'
		if not path.exists(bpmfile):
			sys.exit('data/BPMind.pkl not found')
		if model == 'combined':
			ssmfile = 'data/ssM_hygessi_combined_R0.pkl'
		else:
			ssmfile = 'data/ssM_hygessi_' + model + '_R0.pkl'
		if not path.exists(ssmfile):
			sys.exit(ssmfile+' not found')
		fdr.fdrsampleperm(ssmfile, bpmfile, pval_cutoff, minPath, sample_perms)
		sgmFile = 'data/snpgenemapping_' + str(int(mappingDistance/1000)) + 'kb.pkl'
		snppathwayfile = 'data/snp_pathway_min'+str(minPath)+'_max'+str(maxPath)+'.pkl'
		if not path.exists(sgmFile):
			sys.exit('snp_gene mapping not found')
		if not path.exists(snppathwayfile):
			sys.exit('snp_pathway mapping not found')
		ssm_tmp = ssmfile.split('/')
		ssm_tmp[-1] = 'results_' + ssm_tmp[-1]
		resultsfile = '/'.join(ssm_tmp)
		clr.collectresults(resultsfile,fdrcut,ssmfile,bpmfile,snppathwayfile,sgmFile)
