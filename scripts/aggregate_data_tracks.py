import sys, time
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from glob import glob
from NucApi import Nucleus, EXTERNAL

import numpy as np
from matplotlib import pyplot as plt
from cUtil import dataLayer, apiUtil

# Look at structure contact distance distribution (eff. violation)
# Get regions of abnormally high restraint violation.

# For several cells
# get structures, get original contacts
# calculate contact distance distribution


def combine_data_tracks(nuc_paths, comb_nuc_path, chromosomes, viol_thesh=3.0, particle_sep=100000):

  comb_nuc = Nucleus(comb_nuc_path)

  
  for nuc_path in nuc_paths:

    # Calc region violation track
    nuc = Nucleus(nuc_path)  
    
    name = split(nuc_path)[1].split('_')[1]
    
    viol_track = 'a_cell_%s_region_viol' % name
    contact_track = 'a_cell_%s_contacts' % name
    
    nuc.calcRegionViolDataTrack(group_name=None, track_name=viol_track, structure=None,
                                win_size=1, threshold=viol_thesh, dist_power_law=-0.33)
    
    nuc.calcContactDataTrack(groupName=None, cis=True, trans=True, binSize=particle_sep, trackName=contact_track, countInBin=False)         
    
    nuc.save()
    
    comb_nuc._copyObject(comb_nuc.externalData, nuc.externalData[viol_track], viol_track)
    comb_nuc._copyObject(comb_nuc.externalData, nuc.derivedData[contact_track], contact_track)
     
  comb_nuc.save()
      

if __name__ == '__main__':

  nuc_paths = glob('paper_structs/Cell_*_ambig_10x_100kb.nuc') 
  comb_nuc_path = 'combined2.nuc'
    
  chromosomes = [str(i) for i in range(1,20)] + ['X']
  combine_data_tracks(nuc_paths, comb_nuc_path, chromosomes)
  
