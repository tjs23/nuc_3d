from numpy import array, sqrt, hstack, histogram, log
from matplotlib import pyplot

"""
Null hypothesis from
* permuting data along the sequence
* keeps 1D relationship, clumpiness etc.

Metrics including
* [Done] Radial distribution, wrt nucleus
* Anisotropy (radial)
* Nearest k-neighbour distance
* Density distribution
* Paired distances, e.g. promoter-enhancer
* Corellation coefficient for paired voxel properties
"""

def getRadialDistanceDistrib(nuc, dataTracks, chromosomes=None, models=None,
                             normed=False, nullSamples=100, nBins=25):
  
  if not models:
    models = list(range(nuc.getNumModels()))
  
  deltas = []
  centers = []
  for model in models:
    coords = nuc.getModelCoords(model, chromosomes)
    deltas.append( (coords.max(axis=0) - coords.min(axis=0)).max() )
    centers.append(coords.mean(axis=0))
  
  centers = array(centers)

  rMax = 1.2 * max(deltas)/2.0
  barWidth = 1.0/(nBins*len(dataTracks))
  barOff = 0.0
  
  for typ, code in dataTracks:
    chromoData = nuc.getDataTrack(typ, code, chromosomes)
    nullDists = array([], float)
    allDists = array([], float)
    allVals = array([], float)

    for chromo in chromoData:
      regions, valPairs, annos = chromoData[chromo]
      start, end = nuc.getChromosomeLimits(chromo)
      posWidth = end-start
      nullStep = posWidth/float(nullSamples)
 
      if normed:
        vals = valPairs[:,1]
      else:
        vals = valPairs[:,0]
 
      pos = regions.mean(axis=1)
      
      indices = (pos>start).nonzero() 
      pos = pos[indices]
      vals = vals[indices]
      
      indices = (pos<end).nonzero()
      pos = pos[indices]
      vals = vals[indices]
      
      modelDists = []
 
      for i, model in enumerate(models):
        coords = nuc.getPositionCoords(model, pos, chromo)
        deltas = coords - centers[i]
        dists2 = (deltas*deltas).sum(axis=1)
        modelDists.append(sqrt(dists2)/rMax)

      modelDists = array(modelDists).mean(axis=0)
      allDists = hstack([allDists, modelDists])
      allVals = hstack([allVals, vals])
      
      for j in range(nullSamples):
        pos += nullStep
        pos[(pos>=end).nonzero()] -= posWidth
        
        modelDists = []
        for i, model in enumerate(models):
          coords = nuc.getPositionCoords(model, pos, chromo)
          deltas = coords - centers[i]
          dists2 = (deltas*deltas).sum(axis=1)
          modelDists.append(sqrt(dists2)/rMax)

        modelDists = array(modelDists).mean(axis=0)     
        nullDists = hstack([nullDists, modelDists])
      
 
    print '%s using %d data points' % (code, len(allDists))
    
    color = nuc.getDataTrackColor(typ, code, dType=str, alpha=False)
    
    histData, edgesData = histogram(allDists, bins=nBins, range=(0.0, 1.0),
                                    weights=allVals, normed=True)
                                    
    histNull, edgesNull = histogram(nullDists, bins=nBins, range=(0.0, 1.0),
                                    weights=None, normed=True)
    
    idx = histData.nonzero()
    histData = histData[idx]
    histNull = histNull[idx]
    edgesData = edgesData[idx]

    idx = histNull.nonzero()
    histData = histData[idx]
    histNull = histNull[idx]
    edgesData = edgesData[idx]
    
    #pyplot.hist(allDists, bins=100, range=(0.0, rMax), label=code, 
    #            color=color, normed=True, weights=allVals,
    #            histtype='step', rwidth=dOff)
    
    divergence = histData * log(histData/histNull)
    #divergence = log(histData/histNull)

    pyplot.bar(edgesData + barOff, divergence, width=barWidth,
               color=color, label=code)
               
    barOff += barWidth           
  
  pyplot.legend()
  pyplot.show()
  

def getSphericalAnisotropyDistrib(nuc, dataTrack=None, chromosomes=None,
                                  models=None, normed=True, nBins=45):

  from numpy import zeros, arctan2, arccos, pi, log, arange, meshgrid, vstack, dot
  from mpl_toolkits.basemap import Basemap
      
  if not models:
    models = list(range(nuc.getNumModels()))
  
  deltas = []
  centers = []
  for model in models:
    coords = nuc.getModelCoords(model, chromosomes)
    deltas.append( (coords.max(axis=0) - coords.min(axis=0)).max() )
    centers.append(coords.mean(axis=0))
  
  centers = array(centers)

  if dataTrack is None:
    colorMap = 'RdBu'
    code = 'Sequence'
    if not chromosomes:
      chromosomes = nuc.getChromosomes()
      
    chromoData = {}
    for chromo in chromosomes:
      start, end = nuc.getChromosomeLimits(chromo)
      
      if start == end:
        continue
        
      size = end - start
      step = size / 1000.0
      
      pos = arange(start, end, step) + step/2.0
      values = 2.0 * pos / pos.max()
      values -= 1.0
      
      chromoData[chromo] = pos, values
      
    
  else:
    typ, code = dataTrack
    chromoData = nuc.getDataTrack(typ, code, chromosomes)
    colorMap = 'coolwarm'
  
  
  matrix = zeros((nBins,nBins), float)
  moment = zeros(3, float)

  for chromo in chromoData:
    start, end = nuc.getChromosomeLimits(chromo)
    
    if dataTrack:
      regions, valPairs = chromoData[chromo][:2]
      if normed:
        vals = valPairs[:,1]
      else:
        vals = valPairs[:,0]
      
      pos = regions.mean(axis=1)
    
    else:
      pos, vals = chromoData[chromo]
    
    indices = (pos>start).nonzero() 
    pos = pos[indices]
    vals = vals[indices]
    
    indices = (pos<end).nonzero()
    pos = pos[indices]
    vals = vals[indices]
     
    for i, model in enumerate(models):
      coords = nuc.getPositionCoords(model, pos, chromo)
      deltas = coords - centers[i]
      x, y, z = deltas.T
      
      psi = arccos(z/sqrt(x*x + y*y + z*z))
      phi = arctan2(y, x)
      
      b = array(nBins * psi/pi, int)
      b = b % nBins
 
      a = array(nBins*(1.0 + phi/pi)/2.0, int)
      a = a % nBins
      
      #moment += (coords * vals.reshape(len(coords), 1)).sum(axis=0)
      
      matrix[a,b] += vals
  
  #moment /= sqrt(dot(moment, moment))
  
  matrix -= matrix.min()
  matrix /= matrix.max() 

  fig = pyplot.figure()
  step = 180.0/nBins     
  bm = Basemap(projection='hammer',lon_0=0)
  
  x = arange(-180.0,180.0+step,2*step)
  y = arange(-90.0, 90.0+step,step)

  X, Y = meshgrid(x, y)

  im1 = bm.pcolormesh(X, Y, matrix.T, cmap=colorMap, latlon=True)
  
  bm.drawparallels(arange(-90.0,90.0,30.0), color='k')
  bm.drawmeridians(arange(-180.0,180.0,60.0), color='k')
  
  cb = bm.colorbar(im1, "bottom", size="5%", pad="2%")
  
  pyplot.title('%s spherical anisotropy' % code)
  pyplot.show()
  
  
def getNearestNeighbourDistrib(nuc, typ, code, chromosomes=None, models=None,
                               normed=False, nNeighbours=1, nullSamples=100):

  chromoData = nuc.getDataTrack(typ, code, chromosomes)
  
  
  
def getSpatialDensityDistrib(nuc, typ, code, chromosomes=None,
                             models=None, normed=False, nullSamples=100):

  chromoData = nuc.getDataTrack(typ, code, chromosomes)
  
  
  
def getDensityCorrelation(nuc, typA, codeA, typB, codeB,
                          chromosomes=None, models=None, normed=False):

  chromoDataA = nuc.getDataTrack(typA, codeA, chromosomes)
  chromoDataB = nuc.getDataTrack(typB, codeB, chromosomes)


if __name__ == '__main__':

  import sys
  from os.path import dirname, abspath, join

  thisDir = dirname(abspath(__file__))
  sys.path.remove(thisDir)

  nucDir = dirname(thisDir)
  sys.path.append(nucDir)

  from NucApi import Nucleus
  
  nuc = Nucleus(join(nucDir, 'data', 'examples', 'NXT-33_50k.nuc'))
  #nuc = Nucleus('/home/tjs23/Desktop/old_nuc/NXT-43_calc_1M.nuc')
  
  dataTracks = [#('external', 'RNA-Seq'),
                #('external', 'LaminB'),
                ('external', 'ExpressQ1'),
                ('external', 'ExpressQ2'),
                ('external', 'ExpressQ3'),
                ('external', 'ExpressQ4'),
                ]
  getRadialDistanceDistrib(nuc, dataTracks)
  
  #getSphericalAnisotropyDistrib(nuc, ('external', 'RNA-Seq'))
  #getSphericalAnisotropyDistrib(nuc, None)
  
  
  
