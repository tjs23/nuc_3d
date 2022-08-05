import sys
from os.path import abspath, dirname
sys.path.append(dirname(dirname(abspath(__file__))))

from glob import glob
from matplotlib import pyplot as plt
import numpy as np

from cUtil import apiUtil
from NucApi import Nucleus

def getChromoRegionGenes(chromo, regions, gene_nuc_file='../data/referenceNuc/gene_text.nuc', gene_data_track='gene'):
  
  regions = np.array(regions, np.int32)
    
  nuc = Nucleus(gene_nuc_file)
  
  gene_regions = nuc.getDataTrackRegions(gene_data_track, chromo)
  
  gene_annos = nuc.getDataTrackAnnotations(gene_data_track, chromo)
  
  idx = apiUtil.pairRegionsIntersection(gene_regions, regions, exclude=False, allow_partial=True)
  
  annos = []
  
  if len(idx):
    gene_regions = gene_regions[idx]
    gene_annos = gene_annos[idx]
 
    for region in regions:
      idx2 = apiUtil.pairRegionsIntersection(gene_regions, np.array([region,]), exclude=False, allow_partial=True)
 
      gene_annos2 = gene_annos[idx2]
 
      annos.append(gene_annos2)

  return annos
  

if __name__ == '__main__':

  getChromoRegionGenes('12', [[65000000,72000000],])
