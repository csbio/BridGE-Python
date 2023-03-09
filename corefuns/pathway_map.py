import numpy as np
import pandas as pd
import pickle
from classes import bpmindclass
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from classes import fdrresultsclass
from corefuns import check_BPM_WPM_redundancy as cbwr
from networkx.drawing.nx_agraph import to_agraph
from matplotlib.backends.backend_pdf import PdfPages


def draw_map(project_dir,fdrcut,resultsfile,BPM_group_tmp,WPM_group_tmp,PATH_group_tmp,bpm_limit=20):
	
	# load results file(FDRs)
	#fdr_file = project_dir + '/' + resultsfile
	pklin = open(resultsfile,'rb')
	rfd = pickle.load(pklin)
	pklin.close()
	fdrBPM, fdrWPM, fdrPATH = rfd.fdrbpm2, rfd.fdrwpm2, rfd.fdrpath2
	bpm_pv, wpm_pv, path_pv = rfd.bpm_pv, rfd.wpm_pv, rfd.path_pv
	bpm_ranksum, wpm_ranksum, path_ranksum = rfd.bpm_ranksum, rfd.wpm_ranksum, rfd.path_ranksum
	
	# filter by fdrcut
	ind_bpm = fdrBPM[fdrBPM<=fdrcut].dropna()
	ind_wpm = fdrWPM[fdrWPM<=fdrcut].dropna()
	ind_path = fdrPATH[fdrPATH<=fdrcut].dropna()
	fdrBPM = fdrBPM.iloc[ind_bpm.index]
	fdrWPM = fdrWPM.iloc[ind_wpm.index]
	fdrPATH = fdrPATH.iloc[ind_path.index]
	
	# load BPMind.pkl file
	bpm_file = project_dir + '/BPMind.pkl'
	pklin = open(bpm_file,'rb')
	bpm_data = pickle.load(pklin)
	pklin.close()
	bpm = bpm_data.bpm
	wpm = bpm_data.wpm
	wpm_size = wpm.shape[0]
	bpm_size = bpm.shape[0]

	# find relevant fdrcut id
	if fdrcut > 0.4:
		fdrcut = 0.4
	fdr_group = int(fdrcut/0.05)
	if fdr_group * 0.05 == fdrcut:
		fdr_group = fdr_group - 1
	fdr_th = (fdr_group+1) * 0.05

	# Find and add nodes
	to_draw = []
	used_pathways = []
	## add WPMs
	wpms = fdrWPM[fdrWPM['wpm2'] <= fdr_th ]
	wpms = wpms.sort_values(kind='stable',by='wpm2')
	wpm_groups = np.array(WPM_group_tmp[fdr_group])
	unique_groups = np.unique(wpm_groups).shape[0]
	for g in range(unique_groups):
		tmp = np.where(wpm_groups==g)
		if len(tmp[0]) > 0:
			risk_type = False
			wpm_id = tmp[0][0]
			#wpm_row = wpms.iloc[bpm_id,:]
			xid = wpms.index[wpm_id]
			if xid > wpm_size:
				risk_type = True
				xid = xid - wpm_size
			p1 = wpm['pathway'][xid]
			used_pathways.append(p1)
			if risk_type:
				to_draw.append((p1,p1,'risk'))
			else:
				to_draw.append((p1,p1,'protective'))

	## add PATHs
	path_results = fdrPATH[fdrPATH['path2'] <= fdr_th ]
	path_results = path_results.sort_values(kind='stable',by='path2')
	path_groups = np.array(PATH_group_tmp[fdr_group])
	unique_groups = np.unique(path_groups).shape[0]
	for g in range(unique_groups):
		tmp = np.where(path_groups==g)
		if len(tmp[0]) > 0:
			risk_type = False
			path_id = tmp[0][0]
			#wpm_row = wpms.iloc[bpm_id,:]
			xid = path_results.index[path_id]
			if xid > wpm_size:
				risk_type = True
				xid = xid - wpm_size
			p1 = wpm['pathway'][xid]
			used_pathways.append(p1)
			if risk_type:
				to_draw.append((p1,None,'risk'))
			else:
				to_draw.append((p1,None,'protective'))

	## add BPMs based on the bpm limit
	used_pathways = np.unique(used_pathways)
	bpms = fdrBPM[fdrBPM['bpm2'] <= fdr_th ]
	bpms = bpms.sort_values(kind='stable',by='bpm2')
	bpm_groups = np.array(BPM_group_tmp[fdr_group])
	unique_groups = np.unique(bpm_groups).shape[0]
	### priority is with WPM-PATH associated pathways
	for p in used_pathways:
		x = np.where(bpm['path1names'] == p)
		x = x[0]
		y = np.where(bpm['path2names'] == p)
		y = y[0]
		tmp = np.union1d(x,y)
		for bpm_id in tmp:
			if bpm_id in bpms.index:
				risk_type = False
				if bpm_id > bpm_size:
					risk_type = True
					bpm_id = bpm_id - bpm_size
				p1 = bpm['path1names'][bpm_id]
				p2 = bpm['path2names'][bpm_id]
				if risk_type:
					to_draw.append((p1,p2,'risk'))
				else:
					to_draw.append((p1,p2,'protective'))

	if unique_groups > bpm_limit:
		unique_groups = bpm_limit
	### draw exactly one bpm from each group
	for g in range(unique_groups):
		tmp = np.where(bpm_groups==g)
		if len(tmp[0]) > 0:
			risk_type = False
			bpm_id = tmp[0][0]
			#bpm_row = bpms.iloc[bpm_id,:]
			xid = bpms.index[bpm_id]
			if xid > bpm_size:
				risk_type = True
				xid = xid - bpm_size
			p1 = bpm['path1names'][xid]
			p2 = bpm['path2names'][xid]
			if risk_type:
				to_draw.append((p1,p2,'risk'))
			else:
				to_draw.append((p1,p2,'protective'))

	# Create adjacency matrix
	## first find all pathways
	used_pathways = []
	for t in to_draw:
		p1 = t[0]
		p2 = t[1]
		if p2 != None and p1 != p2:
			used_pathways.append(p2)
		used_pathways.append(p1)
	used_pathways = list(set(used_pathways))
	inds = {}
	nodes_labels = {}
	simple_labels = {}
	for i in range(len(used_pathways)):
		p = used_pathways[i]
		inds[p] = i
		nodes_labels[i] = p
		simple_labels[i] = i
	## now create the matrix
	adj_matrix = np.zeros((len(used_pathways),len(used_pathways)))
	path_array = []
	for t in to_draw:
		p1 = t[0]
		p2 = t[1]
		int_type = t[2]
		if p2 == None:
			path_array.append((p1,int_type))
			continue
		if p1 == p2:
			xid = inds[p1]
			if int_type == 'protective':
				adj_matrix[xid,xid] = 1
			else:
				adj_matrix[xid,xid] = -1
		else:
			xid = inds[p1]
			yid = inds[p2]
			if int_type == 'protective':
				adj_matrix[xid,yid] = 1
				adj_matrix[yid,xid] = 1
			else:
				adj_matrix[xid,yid] = -1
				adj_matrix[yid,xid] = -1

	# Draw
	## Graph(map)
	G = nx.DiGraph()
	added_node_flag = np.zeros((len(used_pathways),))
	### add PATHs with color
	for t in path_array:
		p = t[0]
		int_type = t[1]
		x_id = inds[p]
		if added_node_flag[x_id] == 1: # added before!
			continue
		if int_type == 'protective':
			G.add_nodes_from([(x_id, {"color": "green"})])
		else:
			G.add_nodes_from([(x_id, {"color": "green"})])
		added_node_flag[x_id] = 1

	### add all other the nodes
	for i in range(len(used_pathways)):
		if added_node_flag[i] == 0:
			p = used_pathways[i]
			G.add_nodes_from([(inds[p], {"color": "green"})])
	### add BPMs by adding edges
	for i in range(len(used_pathways)):
		for j in range(i,len(used_pathways)):
			val = adj_matrix[i,j]
			if val == 1:
				p1 = used_pathways[i]
				p2 = used_pathways[j]
				G.add_edge(inds[p1],inds[p2],color='yellow')
				G.add_edge(inds[p2],inds[p1],color='yellow')
			elif val == -1:
				p1 = used_pathways[i]
				p2 = used_pathways[j]
				G.add_edge(inds[p1],inds[p2],color='blue')
				G.add_edge(inds[p2],inds[p1],color='blue')
	### drawing 
	fig_title = 'Non-redundant network map with FDR threshold=' + str(int(fdr_th*100))
	nodes = G.nodes
	node_colors = [nodes[u]['color'] for u in nodes]
	edges = G.edges()
	colors = [G[u][v]['color'] for u,v in edges]
	plt.rcParams["figure.autolayout"] = True
	plt.rcParams["figure.figsize"] = [12, 12]
	fig1 = plt.figure()
	pos = nx.spring_layout(G,k=1.5, iterations=200)
	nx.draw(G, pos, edgelist=edges, edge_color=colors,style='-',node_color=node_colors,nodelist=nodes)
	nx.draw_networkx_labels(G,pos,labels=simple_labels,font_size=10)
	plt.axis('off')
	axis = plt.gca()
	axis.set_xlim([1.5*x for x in axis.get_xlim()])
	axis.set_ylim([1.5*y for y in axis.get_ylim()])
	## add legend
	legend_elements = [Line2D([0], [0], marker='o', color='w', label='pathway',markerfacecolor='g', markersize=15)]
	#legend_elements.append(Line2D([0], [0], marker='o', color='w', label='protective PATH pwathway',markerfacecolor='yellow', markersize=15))
	#legend_elements.append(Line2D([0], [0], marker='o', color='w', label='risk PATH pwathway',markerfacecolor='blue', markersize=15))
	legend_elements.append(Line2D([0], [0], color='yellow', lw=4, label='protective interaction'))
	legend_elements.append(Line2D([0], [0], color='blue', lw=4, label='risk interaction'))
	legend_elements.append(Line2D([0], [0], marker='o', color='w', label='self loops indicate WPM interactions',markerfacecolor='w', markersize=1))
	plt.legend(handles=legend_elements, loc='lower left')



	## add Pathway names to the next page
	fig2 = plt.figure()
	txt = ''
	for i in range(len(used_pathways)):
		tmp = str(i)+': '+nodes_labels[i]
		txt = txt + tmp + '\n'
	plt.axis('off')
	plt.text(0.05,0.05,txt, transform=fig2.transFigure, size=12)
	# find output file name based on the resultsfile
	tmp = resultsfile.split('results_')
	ssmfile = tmp[1]
	tmp = ssmfile.split('.pkl')
	ssmfile = tmp[0]
	pp = PdfPages(project_dir+'/network-map-'+ssmfile+'.pdf')
	fig_nums = plt.get_fignums()
	figs = [plt.figure(n) for n in fig_nums]
	for fig in figs:
		fig.savefig(pp, format='pdf')
	pp.close()






