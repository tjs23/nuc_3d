import sys
from numpy import outer, diag, zeros, cov, ones, dstack, arange, uint32, vstack, int32, array, log, clip
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from util.Image import pixmapToImage
from util.Cluster import kMeansSpread
from cUtil import dataLayer, apiUtil

from NucApi import Nucleus, DERIVED, EXTERNAL

def makeAbCompartmentTracks(nuc, active='Active', binSize=int(5e5), min_obs=0.2):

  groupName = nuc.getDefaultContactGroup()    
  cacheDict = nuc.getCachedContacts(groupName)
  chromosomes = cacheDict.keys()
  
  regionDictA = {}
  valueDictA = {}
  
  regionDictB = {}
  valueDictB = {}
 
  markerTrack = nuc.getRefDataTrackGroup(EXTERNAL, active)
  
  nuc.calcContactVoidRegions(trackName='void', source=EXTERNAL)

  voidTrack = nuc.getRefDataTrackGroup(EXTERNAL, 'void')

  #nuc.normaliseContacts(groupName)
  
  for chromo in chromosomes:

    if chromo not in cacheDict[chromo]:
      continue

    if chromo not in markerTrack:
      continue
    
    print 'Chr', chromo
    
    startPoint, endPoint = nuc.getChromosomeLimits(chromo)
    obs = nuc.getContactMatrix(chromo, chromo, binSize, groupName).astype(float)
     
    obs -= diag(diag(obs))
    obs /= obs.sum()
    n = len(obs)
    
    counts = zeros(n)
    sig = zeros(n)
    
    for i in range(n):
      for j in range(i,n):
        d = j-i
        sig[d] += obs[i,j]
        counts[d] += 1.0
    
    for c, j in enumerate(counts):
      if c:
        sig[j] /= c
    
    sig /= sig.sum()  
    
    exp = zeros((n,n), float)
    for i in range(n):
      exp[i,:i+1] = sig[:i+1][::-1]
      exp[i,i:] = sig[:n-i]
    
    vals = obs.sum(axis=0)
    vals /= vals.sum()
    
    exp *= outer(vals, vals)
    exp /= exp.sum()  
   
    pos = startPoint + arange(0, n+1) * binSize
    pos = vstack([pos[:-1],  pos[:-1]+1]).T.astype(int32)
    idx_nv = ones(len(pos))   
    
    void_regions = array(voidTrack[chromo]['regions'], int32)
    idx_v = apiUtil.pairRegionsIntersection(pos, void_regions, exclude=False, allow_partial=False)
    idx_nv[idx_v] = 0
    idx_v = apiUtil.pairRegionsIntersection(pos+int32(binSize/3), void_regions, exclude=False, allow_partial=False)
    idx_nv[idx_v] = 0
    idx_v = apiUtil.pairRegionsIntersection(pos+int32(2*binSize/3), void_regions, exclude=False, allow_partial=False)
    idx_nv[idx_v] = 0
    idx_v = apiUtil.pairRegionsIntersection(pos+int32(binSize)-1, void_regions, exclude=False, allow_partial=False)
    idx_nv[idx_v] = 0
    
    
    idx_nv = idx_nv.nonzero()[0]
    
    # Null the void
    #obs[idx_v] = 0 
    #obs[:,idx_v] = 0 
    
    obs = obs[idx_nv]
    obs = obs[:,idx_nv]
    
    exp = exp[idx_nv]
    exp = exp[:,idx_nv]
    
    idx = (exp * obs).nonzero()
    z = ((exp * obs) == 0.0).nonzero()
    
    h = obs
    h[idx] /= exp[idx]
    h[idx] = log(h[idx])
    
    h = clip(h, -4.0, 4.0)
    h = cov(h)
    h[z] = 0.0 # Covarience of xero points is irrelevant
    
    #print h.min(), h.mean(), h.max(), vals.min(), vals.mean(), vals.std(), vals.max()
          
    centers, clusters, labels = kMeansSpread(h, 2, verbose=True)
    labels += 1
    
    #v_lim = min_obs * vals.max()
    #labels[(vals < v_lim).nonzero()] = 0
    
    regions = array(markerTrack[chromo]['regions'], int32) # start, end
    values = array(markerTrack[chromo]['values'], float)[:,1] # origValue, normValue
    binned_markers = dataLayer.regionBinValues(regions, values, int32(binSize), startPoint, endPoint)
    
    pos = startPoint + arange(0, n+1) * binSize
    pos = vstack([pos[:-1],  pos[1:]-1]).T.astype(uint32)
    #pos += startPoint - binSize/2
    
    pos = pos[idx_nv]
    binned_markers = binned_markers[idx_nv]
    
    
    idx_a = (labels == 1).nonzero()[0]
    idx_b = (labels == 2).nonzero()[0]
    
    f_a = binned_markers[idx_a].sum()/len(idx_a)
    f_b = binned_markers[idx_b].sum()/len(idx_b)
    
    if f_b > f_a:
      f_a, f_b = f_b, f_a
      idx_a, idx_b = idx_b, idx_a
    
    pos_a = pos[idx_a]
    pos_b = pos[idx_b]
    
    pos_a_merge = []
    pos_b_merge = []
    
    for x, y in pos_a:
      if pos_a_merge:
        if pos_a_merge[-1][1] == x-1:
          pos_a_merge[-1][1] = y
        
        else:
          pos_a_merge.append([x,y])
        
      else:
        pos_a_merge.append([x,y])

    for x, y in pos_b:
      if pos_b_merge:
        if pos_b_merge[-1][1] == x-1:
          pos_b_merge[-1][1] = y
        
        else:
          pos_b_merge.append([x,y])
        
      else:
        pos_b_merge.append([x,y])
    
    
    regionDictA[chromo] = array(pos_a_merge, uint32)
    valueDictA[chromo] = ones((len(pos_a_merge), 2))
 
    regionDictB[chromo] = array(pos_b_merge, uint32)
    valueDictB[chromo] = ones((len(pos_b_merge), 2))
      
    """   
    
    r = zeros(obs.shape, float)
    g = zeros(obs.shape, float)
    b = zeros(obs.shape, float)

    #i = (h > 0.0).nonzero()
    #j = (h < 0.0).nonzero()
    #k = (h == 0.0).nonzero()
    
    i = (labels == 0).nonzero()
    j = (labels == 1).nonzero()   
    k = (labels == 2).nonzero()   
    
    r[i] = 1.0
    g[j] = 1.0
    b[k] = 1.0
    
    #r = g = b = h
                 
    pixmap  = dstack([r,g,b])
    pixmap -= pixmap.min()
    pixmap /= pixmap.max()
    pixmap **= 0.25
    
    img = pixmapToImage(255*pixmap, 'RGB')
    img.show()
    """
    
  nuc.setDataTrack('Comp_A', EXTERNAL, regionDictA, valueDictA, 
                   color=(0.0, 0.5, 1.0), scale=1.0, threshold=0.0,
                   showText=False, shape=0)
                   
  nuc.setDataTrack('Comp_B', EXTERNAL, regionDictB, valueDictB, 
                   color=(1.0, 0.5,0.0), scale=1.0, threshold=0.0,
                   showText=False, shape=0)

  nuc.save()

  
if __name__ == '__main__':

  filePath = 'SLX-7671_hapsort_pop_mm10.nuc'
  
  nuc = Nucleus(filePath)
    
  makeAbCompartmentTracks(nuc, binSize=2.5e5)
