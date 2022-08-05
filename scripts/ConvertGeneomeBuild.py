import sys, os
import numpy as np
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from NucApi import Nucleus, VL_str
from cUtil.apiUtil import pairRegionsIntersection, pointRegionsIntersection

# Get chain file

# chain score chrA sizeA strandA startA endA chrB sizeB strandB startB endB chainId
#   alignLen gapLenA gapLenB
#   alignLen gapLenA gapLenB
#   alignLen

# Go through each chain, indexed by chromo
#   Get main range
#   Step through each equivalent range
#   Store pairs of ranges for each build

# Go through lists of chromo, pos
#  find chain where chromo matches and pos is in range


def convertDataTracksGeneomeBuild(in_nuc_path, chain_file_path, out_nuc_path):

  chain_file_obj = open(chain_file_path)
  
  region_map = {}
  
  for line in chain_file_obj:
    data = line.split()
    
    if not data:
      pos_a = None
      pos_b = None
    
    elif data[0] == 'chain':
      chromo = data[2][3:] # no "chr"
      pos_a = int(data[5])
      pos_b = int(data[10])
      
      if chromo not in region_map:
        region_map[chromo] = []
      
    elif len(data) == 3:
      delta = int(data[0])
      gap_a = int(data[1])
      gap_b = int(data[2])
      region_map[chromo].append((pos_a, pos_a+delta, pos_b, pos_b+delta))
       
      pos_a += gap_a + delta
      pos_b += gap_b + delta
    
    elif len(data) == 1:
      delta = int(data[0])
      region_map[chromo].append((pos_a, pos_a+delta, pos_b, pos_b+delta))

      pos_a = None
      pos_b = None
  
  for chromo in region_map:
    region_map[chromo] = np.array(sorted(region_map[chromo]), np.int32)
  
  nuc = Nucleus(in_nuc_path)
  
  if os.path.exists(out_nuc_path):
    os.unlink(out_nuc_path)
  
  nuc.saveAs(out_nuc_path)
  
  for source in nuc.dataTracks:
    source_group = nuc.dataTracks[source]
  
    for code in source_group:
      code_group = source_group[code]
      print "Processing  ", code
    
      for chromo in code_group:
        print ' ', chromo,
        
        if chromo not in region_map:
          print '(?)',
          continue
        
        chromo_group = code_group[chromo]
        
        orig_regions = np.array(chromo_group['regions'], np.int32)
        n = len(orig_regions)
        
        if not n:
          continue

        map_regions = region_map[chromo]
        map_a = map_regions[:,:2]
        
        idx = pairRegionsIntersection(orig_regions, map_a, exclude=False, allow_partial=False)
        orig_regions = orig_regions[idx] # can be mapped over
        
        if len(idx) != n:
          nuc._setData('values', chromo_group, np.int32, np.array(chromo_group['values'])[idx])
          
          if 'annotations' in chromo_group:
            annos = np.array(chromo_group['annotations'])
            if annos.shape:
              nuc._setData('annotations', chromo_group, VL_str, annos[idx])
        
        map_a_ends  = map_a[:,1]
        idx_starts  = np.searchsorted(map_a_ends, orig_regions[:,0])
        idx_ends    = np.searchsorted(map_a_ends, orig_regions[:,1])
        new_regions = np.empty(orig_regions.shape, np.int32)
        
        for i, j in enumerate(idx_starts):
          k = idx_ends[i]
          p1, p2 = orig_regions[i]
          o1, o2 = p1, p2
          d1 = p2-p1
          
          a1, a2, b1, b2 = map_regions[j]
          p1 = b1 + (p1-a1)

          a1, a2, b1, b2 = map_regions[k]
          p2 = b1 + (p2-a1)
          
          p1, p2 = sorted([p1,p2])
                   
          d2 = p2-p1
          
          if (j != k) and (d2 > 1.5 * d1):
            p2 = min(b1 + d1, map_regions[j][3])
          
          new_regions[i,0] = p1
          new_regions[i,1] = p2
          
        
        nuc._setData('regions', chromo_group, np.int32, new_regions)
        
      print "  "
  
  for code in nuc.interactions:
    print "Processing  interactions", code
    
    for chrA in nuc.interactions[code]:
      for chrB in nuc.interactions[code][chrA]:
        print '%s-%s' % (chrA, chrB),
        
        chromo_group = nuc.interactions[code][chrA][chrB]
        orig_regions = np.array(chromo_group['regions'], np.int32)
        
        map_regions_a = region_map[chrA]
        map_a = map_regions_a[:,:2]
        map_regions_b = region_map[chrB]
        map_b = map_regions_b[:,:2]
        n = len(orig_regions)
        
        mask = np.zeros(n, int)
        
        idx = pairRegionsIntersection(orig_regions[:,:2], map_a, exclude=False, allow_partial=False)
        mask[idx] += 1
        idx = pairRegionsIntersection(orig_regions[:,2:], map_b, exclude=False, allow_partial=False)
        mask[idx] += 1
        
        idx = (mask == 2).nonzero()
        orig_regions = orig_regions[idx] # ChrA and ChrB can be mapped over
        
        if len(idx) != n:
          nuc._setData('values', chromo_group, np.int32, np.array(chromo_group['values'])[idx])
          
          if 'annotations' in chromo_group:
            annos = np.array(chromo_group['annotations'])
            if annos.shape:
              nuc._setData('annotations', chromo_group, VL_str, annos[idx])
        
        new_regions = np.empty(orig_regions.shape, np.int32)
        
        # ChrA
        
        map_ends  = map_a[:,1]
        idx_starts = np.searchsorted(map_ends, orig_regions[:,0])
        idx_ends   = np.searchsorted(map_ends, orig_regions[:,1])
        
        for i, j in enumerate(idx_starts):
          k = idx_ends[i]
          p1, p2 = orig_regions[i,:2]
          o1, o2 = p1, p2
          d1 = p2-p1
          
          a1, a2, b1, b2 = map_regions[j]
          p1 = b1 + (p1-a1)

          a1, a2, b1, b2 = map_regions[k]
          p2 = b1 + (p2-a1)
          
          p1, p2 = sorted([p1,p2])
                   
          d2 = p2-p1
          
          if (j != k) and (d2 > 1.5 * d1):
            p2 = min(b1 + d1, map_regions[j][3])
          
          new_regions[i,0] = p1
          new_regions[i,1] = p2
        
        # ChrB  
        
        map_ends  = map_b[:,1]
        idx_starts = np.searchsorted(map_ends, orig_regions[:,2])
        idx_ends   = np.searchsorted(map_ends, orig_regions[:,3])
        
        for i, j in enumerate(idx_starts):
          k = idx_ends[i]
          p1, p2 = orig_regions[i,2:]
          o1, o2 = p1, p2
          d1 = p2-p1
          
          a1, a2, b1, b2 = map_regions[j]
          p1 = b1 + (p1-a1)

          a1, a2, b1, b2 = map_regions[k]
          p2 = b1 + (p2-a1)
          
          p1, p2 = sorted([p1,p2])
                   
          d2 = p2-p1
          
          if (j != k) and (d2 > 1.5 * d1):
            p2 = min(b1 + d1, map_regions[j][3])
          
          new_regions[i,2] = p1
          new_regions[i,3] = p2
          
        nuc._setData('regions', chromo_group, np.int32, new_regions)
           
    print ''
  
  
  for cat in nuc.contacts:
    for name in nuc.contacts[cat]:
      print "Processing  contacts", name
      group = nuc.contacts[cat][name]
      
      for chrA in group:
        if chrA not in region_map:
          continue
        map_regions_a = region_map[chrA]
        map_a = map_regions_a[:,:2]
          
        for chrB in group[chrA]:
          if chrB not in region_map:
            continue
         
          print '%s-%s' % (chrA, chrB),
          contacts = np.array(group[chrA][chrB])
          new_contacts = nuc._convertContactChainMapping(contacts, region_map[chrA], region_map[chrB])
  
          nuc._setData(chrB, group[chrA], np.int32, new_contacts)
      
      print ''


        
if __name__ == '__main__':

  import sys, glob
  
  chain_file_path = 'data/convert_chains/mm9ToMm10.over.chain'
  
  convertDataTracksGeneomeBuild('/home/tjs23/nucleus/data/referenceNuc/rep_timing_mm9.nuc', chain_file_path, '/home/tjs23/nucleus/data/referenceNuc/rep_timing_mm10.nuc')
  
  #for nuc_file_path in glob.glob('Nextera_01/mm10_100kb/*10x_100kb.nuc'):
  
  #  convertDataTracksGeneomeBuild(nuc_file_path, chain_file_path, nuc_file_path[:-4]+'_mm10.nuc')

 
  #bam_file = '/data/hi-c/SLX-7671_72/SLX-7671_uniq.bam'
  
  #nuc = Nucleus('SLX-7671_pop_mm9.nuc')  
  
  #chain_file = 'nucleus/data/convert_chains/mm10ToMm9.over.chain'
  
  #nuc.importContacts(bam_file, 'SAM', 'SLX-7671', 50000, isSingleCell=False, updateChromos=True, chainMapFile=chain_file)
  
  #nuc.save()
