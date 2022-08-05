import sys, math
import numpy as np
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from NucApi import Nucleus, EXTERNAL
from matplotlib import pyplot as plt


def make_contact_enrichment_tracks(nuc_paths, pop_nuc_path, chromosomes, track_name, bin_size):

  pop_nuc = Nucleus(pop_nuc_path)
  
  group_name = 'all' # pop_nuc.getDefaultContactGroup()
  
  sz = float(bin_size)
  chromo_limits = {}
  chromo_pop_data = {}
  for chromo in chromosomes:
    start, end = pop_nuc.getChromosomeLimits(chromo)
    num_bins = int(end/sz) - int(start/sz)+ 1
    chromo_pop_data[chromo] = np.zeros(num_bins, float)
    chromo_limits[chromo] = start, end 
  
  for i, chr_a in enumerate(chromosomes):
    region_a = chromo_limits[chr_a]
    
    for chr_b in chromosomes[i:]:
      region_b = chromo_limits[chr_b]
     
      matrix = pop_nuc.getContactMatrix(chr_a, chr_b, bin_size, group_name,
                                        region_a, region_b)
      
      if chr_a == chr_b:
        for j in range(len(matrix)):
          matrix[j,j] = 0
      
      chromo_pop_data[chr_a] += matrix.sum(axis=1)
      chromo_pop_data[chr_b] += matrix.T.sum(axis=1)
      
  for chromo in chromosomes:
    chromo_pop_data[chromo] /= chromo_pop_data[chromo].sum()
        
  for nuc_path in nuc_paths:
    nuc = Nucleus(nuc_path)
    group_name = nuc.getDefaultContactGroup()
    
    chromo_data = {}
    for chromo in chromosomes:
      chromo_data[chromo] = np.zeros(chromo_pop_data[chromo].shape, float)
    
    for i, chr_a in enumerate(chromosomes):
      region_a = chromo_limits[chr_a]
 
      for chr_b in chromosomes[i:]:
        region_b = chromo_limits[chr_b]
       
        matrix = nuc.getContactMatrix(chr_a, chr_b, bin_size, group_name,
                                      region_a, region_b)
        
        chromo_data[chr_a] += matrix.sum(axis=1)
        chromo_data[chr_b] += matrix.T.sum(axis=1)
    
    b_min = 1e99
    b_max = -1e99
    chromo_bias = {}
    for chromo in chromosomes:
      obs = chromo_data[chromo]
      
      obs /= obs.sum() or 1.0
      
      exp = chromo_pop_data[chromo]
      
      idx = (obs*exp).nonzero()
      
      bias = np.zeros(obs.shape, float)
  
      bias[idx] = obs[idx] + np.log2(obs[idx]/exp[idx])
      
      chromo_bias[chromo] = bias
      
      if b_min > bias.min():
        b_min = bias.min()

      if b_max < bias.max():
        b_max = bias.max()
      
    
    region_dict = {}
    value_dict = {}
    
    for chromo in chromosomes:
      bias = chromo_bias[chromo] 
      start, end = chromo_limits[chromo]
      
      start = sz * int(start/sz)
      starts = np.arange(start, start+bin_size*len(bias), bin_size)
      ends = starts + bin_size - 1
      
      norm_bias = (bias-b_min)/(b_max-b_min)
      
      value_dict[chromo] = np.array([bias, norm_bias]).T
      region_dict[chromo] = np.array([starts, ends]).T
      
    nuc.setDataTrack(track_name, EXTERNAL, region_dict, value_dict, 
                     color=(0.0, 0.5, 1.0), scale=1.0, threshold=0.0,
                     showText=False, shape=0)
   
    nuc.save()
    
   
if __name__ == '__main__':
  
  from glob import glob
  
  #pop_nuc_path = '/data/hi-c/pop_nuc_files/SLX-7671_hapsort_pop_mm10.nuc'
  pop_nuc_path = 'combined2.nuc'
  
  nuc_paths = glob('paper_structs/*_ambig_10x_100kb.nuc')
  
  chromosomes = [str(i) for i in range(1,20)] + ['X']
    
  make_contact_enrichment_tracks(nuc_paths, pop_nuc_path, chromosomes, 'cont_enrich', bin_size=int(4e5))
