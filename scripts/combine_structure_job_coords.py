import numpy as np
import sys
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from glob import glob
from NucApi import Nucleus


def combine_job_coords(nuc_file_path, in_file_paths):
  
  chromosomes = sorted([str(x) for x in range(1,20)] + ['X'])
    
  nuc = Nucleus(nuc_file_path)
    
  coords = np.array([np.load(f) for f in sorted(in_file_paths)])
    
  nuc.setAllCoords(coords, chromosomes, structure='0')      
  
  nuc.modelAlign(chromosomes=chromosomes, structure='0')
  
  nuc.save()


if __name__ == '__main__':
  
  #nuc_file_paths = sorted(glob('/home/tjs23/nucleus/P2-E8-PA_test/P2-E8-PA_test_*_20x_400kb.nuc'))
  nuc_file_paths = sorted(glob('/home/tjs23/nucleus/P2-E8-PA_test/P2-E8-PA_test_noise_*_10x_100kb.nuc'))
  #nuc_file_paths = sorted(glob('/home/tjs23/nucleus/P2-E8-PA_test/P2-E8-PA_test_*_10x_100kb.nuc'))
  #nuc_file_paths = sorted(glob('/home/tjs23/nucleus/paper_structs/*_20x_400kb.nuc'))
  #nuc_file_paths = sorted(glob('/home/tjs23/nucleus/paper_structs/*_ambig_10x_100kb.nuc'))
  
  for i, nuc_file_path in enumerate(nuc_file_paths):
    
    dir_name, file_name = split(nuc_file_path)
    file_root, file_ext = splitext(file_name)
    
    coord_path = join(dir_name, 'coords', file_root)
    
    in_file_paths = sorted(glob(coord_path + '_job_coords_*.npy'))
  
    print i+1, nuc_file_path
    for p in in_file_paths:
      print '  ', p
  
    combine_job_coords(nuc_file_path, in_file_paths)
  
