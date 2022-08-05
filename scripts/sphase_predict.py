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


def combine_rt_with_violation(nuc_paths, rt_nuc_path, rt_track, chromosomes, rt_quantile=0.9,
                              viol_thesh=3.0, viol_track='region_viols', out_track='rt_viol'):

  rt_nuc = Nucleus(rt_nuc_path)

  
  # Get absolute RT threshold
  rt_vals = [rt_nuc.getDataTrackValues(rt_track, chromo)[:,0] for chromo in chromosomes]
  rt_vals = np.concatenate(rt_vals, axis=0)
  rt_vals = rt_vals[rt_vals.argsort()]
  ref_idx = int(rt_quantile*len(rt_vals))
  ref_val = rt_vals[ref_idx]
  
  print "ref RT^ val", ref_val
  
  for nuc_path in nuc_paths:
    
    print nuc_path
    
    region_dict = {}
    values_dict = {}
    
    # Calc region violation track
    nuc = Nucleus(nuc_path)  
    nuc.calcRegionViolDataTrack(group_name=None, track_name=viol_track, structure=None,
                                win_size=1, threshold=viol_thesh, dist_power_law=-0.33)
    
    bin_size = nuc.getBackboneSpacing()
    
    for chromo in chromosomes:
      
      print chromo
      
      # Get chromo RT values above threshold
      rt_vals = rt_nuc.getDataTrackValues(rt_track, chromo)[:,0]
      rt_regions = rt_nuc.getDataTrackRegions(rt_track, chromo)
      idx = (rt_vals > ref_val).nonzero()[0]
      
      if not len(idx):
        continue
      
      rt_vals = rt_vals[idx]
      rt_regions = rt_regions[idx]
      
      # Bin hight RT values to beads, take non-zeros
      s, e = nuc.getChromosomeLimits(chromo)
      s = bin_size * int(s/bin_size)
      e = bin_size * int(e/bin_size) + bin_size
      binned_rt = dataLayer.regionBinValues(rt_regions, rt_vals, np.int32(bin_size), s, e)
      bin_starts = np.arange(s, e, bin_size)
      
      idx = binned_rt.nonzero()[0]
      bin_starts  = bin_starts[idx]
      
      if (idx is None) or not len(idx):
        continue
      
      # Find intersection non-zero high RT beads with violation regions     
      viol_regions = nuc.getDataTrackRegions(viol_track, chromo)
      
      if viol_regions is None:
        continue
      
      bin_regions = np.array([bin_starts, bin_starts+bin_size], np.int32).T
      idx = apiUtil. pairRegionsIntersection(bin_regions, viol_regions,
                                             exclude=False, allow_partial=True,
                                             region_indices=False)
      
      if (idx is None) or not len(idx):
        continue
       
      region_dict[chromo] = bin_regions[idx]
      values_dict[chromo] = np.ones(region_dict[chromo].shape)

    # Store data track
    if values_dict:
      nuc.setDataTrack(out_track, EXTERNAL, region_dict, values_dict, color='#FF00E0', threshold=0.0)          
      nuc.save()
      

def plot_restraint_distance_distrib(nuc_paths, structure=None, particle_size=int(1e5), win_size=5,
                                    dist_power_law=-0.33, dist_lower=0.8, dist_upper=1.2):
  
  # See API calcRegionContactDistDataTrack
  chromosomes = [str(x) for x in range(1,20)] + ['X']
  
  for nuc_path in nuc_paths:
    
    print nuc_path
    
    dists = []
    nuc = Nucleus(nuc_path)
    
    particle_size = nuc.getBackboneSpacing(structure)
    num_models = nuc.getNumModels()
    models = range(num_models)
    partic_group = nuc._getParticleGroup(structure)
    group_name = nuc.getDefaultContactGroup()
    contact_dict = nuc.getCachedContacts(group_name)
    
    restrDict, posDict, ambigDict, bboneDict, binLstDict = apiUtil.calcRestraints(chromosomes, [], contact_dict, True,
                                                                          [particle_size, particle_size], {}, 1.0, # scale not used
                                                                           dist_power_law, dist_lower, dist_upper, 1)
    for chr_a in contact_dict:
      chr_pos_a = partic_group[chr_a]['positions']
      coords_a = nuc.getChromoCoords(chr_a)
      
      if coords_a is None:
        continue

      n_a = coords_a.shape[1]
      
      for chr_b in contact_dict[chr_a]:
        chr_pos_b = partic_group[chr_b]['positions']
        coords_b = nuc.getChromoCoords(chr_b)
        
        if coords_b is None:
          continue
          
        n_b = coords_b.shape[1]             
        restraints = restrDict[chr_a][chr_b] # [i, j, weight, target, lower, upper] 
        
        idx_a = np.array(restraints[0], int)
        idx_b = np.array(restraints[1], int)
        target = restraints[3]
        
        if not len(idx_a):
          continue
        
        bead_coords_a = coords_a[:,idx_a]
        bead_coords_b = coords_b[:,idx_b]

        for m in models:
        
          deltas = bead_coords_a[m] - bead_coords_b[m]
          model_dists = np.sqrt((deltas*deltas).sum(axis=1))
          
          restraint_deltas = model_dists - target
          
          region_means = []
          
          for i in range(win_size):
            j = win_size - i
            region_means.append(restraint_deltas[i:-j])
          
          region_means = np.max(region_means, axis=0)
          
          dists.append(region_means)

    dists = np.concatenate(dists, axis=0)

    #anchor = np.argsort(dists)[int(0.05 * len(dists))]

    h, e = np.histogram(dists, bins=200, normed=True)

    plt.plot(e[1:], h, label=split(nuc_path)[1])
  
  plt.legend()
  plt.show()

if __name__ == '__main__':

  nuc_paths = glob('paper_structs/*_ambig_10x_100kb.nuc') 
  
  #plot_restraint_distance_distrib(nuc_paths)
  
  rt_nuc_path = 'data/referenceNuc/rep_timing_mm10.nuc'
  rt_track = 'rep_timing'
    
  chromosomes = [str(i) for i in range(1,20)] + ['X']
  
  combine_rt_with_violation(nuc_paths, rt_nuc_path, rt_track, chromosomes, 0.9)
  
