## main function of BridGE

import sys
from os import path
from datatools import plink2pkl as p2p
from datatools import bindataa as ba
from datatools import msigdb2pkl as msig2p
from datatools import mapsnp2gene as snp2gene
from datatools import snppathway as snpp
from datatools import bpmind as bpm
from corefuns import matrix_operations_par as ci
from corefuns import genstats_perm as gs
from corefuns import fdrsampleperm as fdr
from corefuns import collectresults as cl
import datetime



if __name__ == '__main__':
	# Default parameters defined
	job = ''
	plinkfile = ''
	#project_dir = 'data'
	genesets='refdata/c2.cp.v7.1' 
	gene_annotation = 'refdata/glist-hg38' 
	mappingDistance = 50000
	minPath = 10
	maxPath = 300
	alpha1 = 0.05
	alpha2 = 0.05
	n_workers = 4
	sample_perms = 10
	binaryNetwork = False
	snpPerms = 10000
	i = -1
	r = -1
	pval_cutoff = 0.05
	fdrcut = 0.25
	densitycutoff = None
	ssmfile = None
	model = None
#	snppathwayfile = 'data/snp_pathway_min10_max300.pkl'
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
			elif o == '--i':
				i = int(a)
			elif o == '--fdrcut':
				fdrcut = float(a)
			elif o == '--snpPathFile':
				snppathwayfile = a
			elif o == '--projectDir':
				project_dir = a
			elif o == '--densityCutoff':
				densitycutoff = float(a)
			elif o == '--ssmfile':
				ssmfile = a
			elif o == '--R':
				r = int(a) 


	# First module: Data processing
	if job == 'DataProcess':
		print('data processing...')
		sys.stdout.flush()
		## data processing step

		## convert plinkfile to pickle
		if plinkfile == '':
			sys.exit('plinkFile not provided')
		rawfile = project_dir + '/intermediate/' + plinkfile + '.raw'
		bimfile = project_dir + '/intermediate/' + plinkfile + '.bim'
		famfile = project_dir + '/intermediate/' + plinkfile + '.fam'
		if not path.exists(rawfile) or not path.exists(bimfile) or not path.exists(famfile):
			sys.exit('plinkFiles do not exist')
		finalfile = project_dir + '/intermediate/' + plinkfile + '.pkl'
		p2p.plink2pkl(rawfile, bimfile, famfile, finalfile)
		## converting snp data assuming different disease models
		ba.bindataa(project_dir,finalfile,'r')
		ba.bindataa(project_dir,finalfile,'d')
		symbolsFile = project_dir + '/raw/' + genesets + '.symbols.gmt'
		entrezFile = project_dir + '/raw/'+ genesets + '.entrez.gmt'
		if not path.exists(symbolsFile) or not path.exists(entrezFile):
			sys.exit('genesets do not exist')
		## prepare gene set information
		msig2p.msigdb2pkl(symbolsFile, entrezFile)
		gene_annotation	= project_dir +	'/raw/' + gene_annotation
		sys.stdout.flush()
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
		geneset_pkl = project_dir + '/intermediate/' + genesets + '.pkl'
		outfile = snpp.snppathway(finalfile, sgmfile, geneset_pkl, minPath, maxPath)
		bpm.bpmind(outfile)
	elif job == 'ComputeInteraction':
		## Validating input parameters
		if not (model == 'RR' or model == 'RD' or model == 'DD' or model == 'combined'):
			sys.exit('wrong model')
		if not path.exists(project_dir+'/intermediate/SNPdataAD.pkl'):
			sys.exit(project_dir+'/intermediate/SNPdataAD.pkl not found')
		if not path.exists(project_dir+'/intermediate/SNPdataAR.pkl'):
			sys.exit(project_dir+'/intermediate/SNPdataAR.pkl not found')
		if r < 0 :
			if model == 'combined':
				ci.combine(project_dir,alpha1,alpha2,n_workers,i)
			else:
				ci.run(project_dir,model,alpha1,alpha2,n_workers,i)
		else:
			for i in range(r+1):
				if model == 'combined':
					ci.combine(project_dir,alpha1,alpha2,n_workers,i)
				else:
					ci.run(project_dir,model,alpha1,alpha2,n_workers,i)

	elif job == 'ComputeStats':
		if not (model == 'RR' or model == 'RD' or model == 'DD' or model == 'combined' or ssmfile != None):
			sys.exit('wrong model')
		bpmfile = project_dir+'/intermediate/BPMind.pkl'
		if not path.exists(bpmfile):
			sys.exit(project_dir+'/intermediate/BPMind.pkl not found')
		if not path.exists(project_dir+'/intermediate/SNPdataAD.pkl'):
			sys.exit(project_dir+'/intermediate/SNPdataAD.pkl not found')
		if not path.exists(project_dir+'/intermediate/SNPdataAR.pkl'):
			sys.exit(project_dir+'/intermediate/SNPdataAR.pkl not found')
		if ssmfile == None:
			if r < 0:
				if model == 'combined':
					#ci.combine(alpha1,alpha2,n_workers,R)
					# build ssM file name
					ssmfile = project_dir+'/intermediate/ssM_mhygessi_combined_R'+ str(i) + '.pkl'
				else:
					#ci.run(model,alpha1,alpha2,n_workers,R)	
					ssmfile = project_dir+'/intermediate/ssM_mhygessi_' + model + '_R'+ str(i) + '.pkl'
				gs.genstats(ssmfile,bpmfile,binaryNetwork,snpPerms,minPath,n_workers,densitycutoff)
			else:
				for i in range(r+1):
					if model == 'combined':
						#ci.combine(alpha1,alpha2,n_workers,R)
						# build ssM file name
						ssmfile = project_dir+'/intermediate/ssM_mhygessi_combined_R'+ str(i) + '.pkl'
					else:
						#ci.run(model,alpha1,alpha2,n_workers,R)	
						ssmfile = project_dir+'/intermediate/ssM_mhygessi_' + model + '_R'+ str(i) + '.pkl'
					gs.genstats(ssmfile,bpmfile,binaryNetwork,snpPerms,minPath,n_workers,densitycutoff)
		else:
			ssmfile = project_dir+'/intermediate/'+ssmfile
			gs.genstats(ssmfile,bpmfile,binaryNetwork,snpPerms,minPath,n_workers,densitycutoff)

	elif job == 'ComputeFDR':
		bpmfile = project_dir+'/intermediate/BPMind.pkl'
		if not path.exists(bpmfile):
			sys.exit(project_dir+'/intermediate/BPMind.pkl not found')
		if ssmfile == None:
			if model == 'combined':
				ssmfile = project_dir+'/intermediate/ssM_mhygessi_combined_R0.pkl'
			else:
				ssmfile = project_dir+'/intermediate/ssM_mhygessi_' + model + '_R0.pkl'
		else:
			ssmfile = project_dir+'/intermediate/'+ssmfile
		if not path.exists(ssmfile):
			sys.exit(ssmfile+' not found')
		fdr.fdrsampleperm(ssmfile, bpmfile, pval_cutoff, minPath, sample_perms)
	elif job == 'Summarize':
		bpmfile = project_dir+'/intermediate/BPMind.pkl'
		snppathwayfile = project_dir+ '/intermediate/' + snppathwayfile
		#snpgenemappingfile  = 'data/snpgenemapping_50kb.pkl'
		snpgenemappingfile = project_dir+'/intermediate/snpgenemapping_' + str(int(mappingDistance/1000)) + 'kb.pkl'
		validationfile = None
		if not path.exists(bpmfile):
			sys.exit('bpm file not found at:'+bpmfile)
		if not path.exists(snppathwayfile):
			sys.exit('snp-pathway mapping file not found')
		if not path.exists(snpgenemappingfile):
			sys.exit('snpgenemappingfile not found, check mapping Distance arg')
		if ssmfile == None:
			if model == 'combined':
				ssmfile = project_dir+'/intermediate/ssM_mhygessi_combined_R0.pkl'
				resultsfile = project_dir+'/intermediate/results_ssM_mhygessi_combined_R0.pkl'
			else:
				ssmfile = project_dir+'/intermediate/ssM_mhygessi_' + model + '_R0.pkl'
				resultsfile = project_dir+'/intermediate/results_ssM_mhygessi_'+model +'_R0.pkl'
		else:
			resultsfile = project_dir+'/intermediate/results_'+ ssmfile
			ssmfile = project_dir+'/intermediate/'+ssmfile
		if not path.exists(ssmfile):
			sys.exit('interaction file not found, check model arg')
		if not path.exists(resultsfile):
			sys.exit('results file from analysis module not found')
		cl.collectresults(resultsfile,fdrcut,ssmfile,bpmfile,snppathwayfile,snpgenemappingfile, densitycutoff)
