import sys, os, math
import numpy as np

from collections import defaultdict
from glob import glob
from matplotlib import pyplot as plt

from os.path import abspath, dirname, split
sys.path.append(dirname(dirname(abspath(__file__))))

from NucApi import Nucleus

q = 5

track_name = 'RNAseq_Nuc'

nuc_path  = '../data/referenceNuc/Mouse_Hap_ESC_mm10.nuc'

nuc = Nucleus(nuc_path)

group, source = nuc._getDataTrackSource(track_name)

track_group = group[track_name]

all_values = []
for chromo in track_group:
  values  = nuc.getDataTrackValues(track_name, chromo)[:,0]
  values = values[values.nonzero()]
  all_values.append(values)

all_values = np.concatenate(all_values)
sort_vals = all_values[all_values.argsort()]

n = len(sort_vals)
thresh = [int(n*(i+1)/5)-1 for i in range(q)]
thresh = [sort_vals[i] for i in thresh]

print thresh, all_values.mean()

quant_region_dicts = [{} for i in range(q)]
quant_value_dicts =  [{} for i in range(q)]

for chromo in track_group:
  
  v = np.array(nuc.getDataTrackValues(track_name, chromo))
  r = np.array(nuc.getDataTrackRegions(track_name, chromo))
  values  = v[:,0]
   
  for i in range(q):
    
    idx = (values <= thresh[i]).nonzero()
    idx2 = (values > thresh[i]).nonzero()
    
    #print chromo, i, len(idx[0]), len(idx2[0])
    
    quant_value_dicts[i][chromo]  = v[idx]
    quant_region_dicts[i][chromo] = r[idx]
    
    values = values[idx2]
    v = v[idx2]
    r = r[idx2]

for i in range(q):
  code = '%s_Q%d' % (track_name, i+1)
  region_dict = quant_region_dicts[i]
  value_dict = quant_value_dicts[i]
  nuc.setDataTrack(code, source, region_dict, value_dict)

nuc.save()
