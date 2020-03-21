from scipy import io
from classes import InteractionNetwork
from classes import snpsetclass as snps
from classes import SNPdataclass as snpc
import pickle
import numpy as np
import pandas as pd
from classes import bpmindclass as bpmc

input_bpm = 'bpm.mat'


bb = io.loadmat(input_bpm)
ind1 = bb['ind1']
path1idx = bb['path1idx']
ind2 = bb['ind2']
path2idx = bb['path2idx']
ind1size = bb['ind1size']
ind2size = bb['ind2size']
size = bb['size']

d = {'ind1': [ind1], 'ind2': [ind2], 'path1idx': [path1idx], 'path2idx': [path2idx], 'ind1size': [ind1size], 'ind2size': [ind2size], 'size': [size] }
bpm = pd.DataFrame(d)

input_wpm = 'wpm.mat'
wb = io.loadmat(input_wpm)

ind = wb['ind']
indsize = wb['indsize']
pathway = wb['pathway']
size = wb['size']

d={'ind': [ind], 'indsize': [indsize] , 'pathway': [pathway], 'size': [size]}
wpm = pd.DataFrame(d)


b = bpmc.bpmindclass(bpm,wpm)
final = open('bpmind.pkl', 'wb')
pickle.dump(b, final)
final.close()



