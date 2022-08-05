import sys

"""
Script to make chromosome colour data tracks

Input 4/5 histone marker data tracks in a nuc file

Quantile normalise the data track with max() to tie-break

Bin at e.g. 100 kb

Use k-means clustering to split the data into 4/5 best-separated classes


"""
from os.path import abspath, dirname, join, splitext, split
sys.path.append(abspath(dirname(dirname(__file__))))

from math import ceil
from numpy import array, arange, zeros, int32, uint32, concatenate, empty, log
from scipy.cluster import vq

from NucApi import Nucleus, EXTERNAL
from cUtil import dataLayer
from util.Cluster import kMeansSpread, kMeans

from random import randint, sample
from numpy import array, cov, diag, dot, linalg, ones
from numpy import outer, random, sqrt, vstack, zeros

COLORS = ((0.0, 0.6, 1.0), # Active
          (0.85, 0.85, 0.0), # Silenced
          (1.0, 0.0, 0.0), # Inactive
          (0.5, 0.5, 0.5)) # Null

LABELS = ('Active', 'Polycomb', 'Inactive', 'Null')

def makeChromoColorDataTracks(nuc, active, silenced, inactive, k=4, binSize=1e5, normChromo=False):
  
  from scipy import stats
  dataTrackCodes = active + silenced + inactive
  
  particGroup = nuc._getParticleGroup()
  
  binSize = int(binSize)
  
  # Collect chromosomes common to all data tracks
  chromos = set(nuc.chromosomes.keys())
  for code in dataTrackCodes:
    refLayer = nuc.getRefDataTrackGroup(EXTERNAL, code)   
    chromos = chromos & set(refLayer.keys())
  
  # Collect binned track data as a dingle array
  data = None
  posDict = {}
  chromoRanges = {}
  
  t = len(dataTrackCodes)
  
  i = 0    
  for chromo in chromos:
    startPoint, endPoint = nuc.getChromosomeLimits(chromo)
     
    if not endPoint:
      continue
    
    vectors = None
    
    for j, code in enumerate(dataTrackCodes):
      refLayer = nuc.getRefDataTrackGroup(EXTERNAL, code)
      regions = array(refLayer[chromo]['regions'], int32) # start, end
      values = array(refLayer[chromo]['values'], float)[:,1] # origValue, normValue

      binned = dataLayer.regionBinValues(regions, values, int32(binSize), startPoint, endPoint)
      
      if normChromo:
        binned = stats.rankdata(binned, method='min')
      
      binned /= binned.max()
        
      if vectors is None:
        vectors = zeros((t, len(binned)), float)
      
      vectors[j] = binned
    
    if data is None:
      data = vectors
    else:
      data = concatenate([data, vectors], axis=1)
    
    n = len(vectors[0])
    pos = binSize * arange(0, n+1, 1, uint32)
    pos += startPoint - binSize/2
    pos = vstack([pos[:-1],  pos[1:]-1]).T
    posDict[chromo] = pos
    chromoRanges[chromo] = (i, i+n)
    i += n  
  
  m = len(data[0])
  quantData = empty((m, t), float)
  if normChromo:
    for j in range(t):
      quantData[:,j] = data[j]
  
  else:
    for j in range(t):
      ranks = stats.rankdata(data[j], method='min')
      ranks /= ranks.max()
      quantData[:,j] = ranks
   
  centers, clusters, labels = kMeansSpread(quantData, k, verbose=True)
  cat_map = {}
  
  for i, cat in enumerate((active, silenced, inactive)):
    data = quantData[:,dataTrackCodes.index(cat[-1])]
    
    dps = []
    for c in range(k):
      cat_vec = array((labels == c), float)
      dps.append( dot(data, cat_vec) )

    j = dps.index(max(dps) )
    print 'Colour', i, 'Cluster', j, dps
    cat_map[j] = i
    
  
  for j in range(k):
    cat = cat_map.get(j, 3)
    regionDict = {}
    valueDict = {}
    
    for chromo in posDict:
      a, b = chromoRanges[chromo]
      idx = (labels[a:b] == j).nonzero()
      pos = posDict[chromo][idx]
    
      regionDict[chromo] = pos
      valueDict[chromo] = ones((len(pos), 2))

    code = LABELS[cat]
    nuc.setDataTrack(code, EXTERNAL, regionDict, valueDict, 
                     color=COLORS[cat], scale=1.0, threshold=0.0,
                     showText=True, shape=0)
  
  nuc.save()
  

if __name__ == '__main__':

  filePath = 'data/referenceNuc/Mouse_Hap_ESC_mm10.nuc'
  
  nuc = Nucleus(filePath)
    
  makeChromoColorDataTracks(nuc, active=['Laue_H3K4me3_hap'], silenced=['Laue_H3K27me3_hap'], inactive=['Laue_H3K9me3_hap'], binSize=100000)
