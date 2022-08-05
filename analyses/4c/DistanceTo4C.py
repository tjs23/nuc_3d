from math import log, sqrt
import os, sys

sys.path.append('/home/tjs23/chromoStruct')
  
from memops.qtgui.Application import Application
from memops.qtgui.BasePopup import BasePopup
from memops.qtgui.Graph import Graph, GraphDataSet, GraphAxis, DEFAULT_COLORS
from memops.qtgui.PixmapPlot import PixmapPlot

from chromoStruct.analyses.Util import getMapabilityDict, readVariableWigFile
from chromoStruct.core.FileIO import loadDataSet

CSM_DIR = '/home/tjs23/chromoVista/data/paperCells/'
MAPABILITY = '/home/tjs23/chromoStruct/data/mapability/mm_ref_MGSCv37_chrX_readFrags.dat'

BIN_SIZE = 1e5
#LOCUS = 100626856 # Tsix
#LOCUS = 100655714 # Xist - Atrx (DNA repair),  NUP62CL (Nucleoporin)
#LOCUS = 103460373
#geneName = 'Xist'

LOCUS = 71271931 # MECP2 74026821
geneName = 'MECP2'
dataFile4c = '/home/tjs23/chromoStruct/data/4c/4C_Splinter_GSM730551_MeCP2_Xa.wig'

LOCUS = 148700000 # Jarid1C mm37
geneName = 'Jarid1C'
dataFile4c = '/home/tjs23/chromoStruct/data/4c/4C_Splinter_GSM730554_Jarid1C_Xa.wig'

CHROMOSOME = 'X'
markers = [(148700000/170e6, 'Jarid1C'), (71271931/170e6, 'MeCP2')]


# Load the structure data

fileNames = [os.path.join(CSM_DIR, f) for f in os.listdir(CSM_DIR)]
fileNames = [f for f in fileNames if os.path.isfile(f)]
#fileNames = [f for f in fileNames if f.endswith('_ADC_X_200.csm')]
fileNames = [f for f in fileNames if f.endswith('X_200.csm')]
fileNames.sort()

for fileName in fileNames:
  print(fileName)

print("Reading mapability")

mapability = getMapabilityDict(MAPABILITY, BIN_SIZE, CHROMOSOME)

print("Reading 4c data")

dataList4c = readVariableWigFile(dataFile4c)[CHROMOSOME]
dataDict4c = {}

for pos, value in dataList4c:
  bin = BIN_SIZE * int(pos//BIN_SIZE)
  
  if bin in dataDict4c:
    dataDict4c[bin] += value
  else:
    dataDict4c[bin] = value

maxVal = log(float(max(dataDict4c.values())))
for bin in dataDict4c:
  value = log(dataDict4c[bin])
  dataDict4c[bin] = value/maxVal

#

groupDists = {}
sampleDict = {}

for fileName in fileNames:
  print("Loading %s..." % fileName)
  sample = loadDataSet(fileName)
  sampleDict[fileName] = sample

#

d = set()
distDict = {}
cisDict = {}
distDictCell = {}
cisDictCell = {}

for fileName in fileNames:
  distDictCell[fileName] = {}
  cisDictCell[fileName] = {}

  sample = sampleDict[fileName]
  models = range(sample.getNumModels())
  m = float(len(models))
  nodeDict = sample.chromoNodeDict

  cisLoci = sample.getCisInteractionLoci(LOCUS, 5e6, CHROMOSOME)
  for locus in cisLoci:
    key = BIN_SIZE * int(locus//BIN_SIZE)
 
    if key in cisDict:
      cisDict[key] += 1
    else:
      cisDict[key] = 1
 
    if key in cisDictCell[fileName]:
      cisDictCell[fileName][key] += 1
    else:
      cisDictCell[fileName][key] = 1
  
  nodes = [n for n in sample.chromoNodes if n.chromosome == CHROMOSOME]
  nodes = [n for n in nodes if (n.locus) % BIN_SIZE == 0] # Backbone only
  nodes = [n for n in nodes if mapability[n.locus // BIN_SIZE] > 0.5] # Only mapped

  refNode = nodeDict[CHROMOSOME][int(BIN_SIZE * (LOCUS//BIN_SIZE))]
  dists = []

  for node in nodes:
    dist = 0.0

    for model in models:
      x1, y1, z1 = refNode.coords[model]
      x2, y2, z2 = node.coords[model]

      dx = x1-x2
      dy = y1-y2
      dz = z1-z2

      dist += sqrt(dx*dx + dy*dy + dz*dz)

    dist /= m
    
    locus = int(node.locus)
    
    if locus not in distDict:
      distDict[locus] = []
    
    if locus not in distDictCell[fileName]:
      distDictCell[fileName][locus] = []
    
    distDict[locus].append(dist)
    distDictCell[fileName][locus].append(dist)
    d.add(dist)

dMeanDict = {}
dists = []
dAveDict = {}
loci = distDict.keys()
loci.sort()
for locus in loci:
  dist = distDict[locus]
  dAve = sum(dist)/float(len(dist))
  dists.append((locus, dAve))
  dAveDict[locus] = dAve

for fileName in distDictCell:
  for locus in distDictCell[fileName]:
    dists = distDictCell[fileName][locus]
    dAve = sum(dists)/float(len(dists))
    distDictCell[fileName][locus] = dAve
  values = distDictCell[fileName].values()
  dMeanDict[fileName] = sum(values)/float(len(values))

# Make graphs

app = Application()
popup = BasePopup(title='Pseudo-4C Graph')
popup.setSize(1600, 600)

from numpy import array, empty, ones, uint8

nCells = len(fileNames)
n = int(170e6/BIN_SIZE)
m = 62 * (nCells + 1) - 1


# 4C band

data4c = empty((10 , n, 4), float)
data4c[:,:] = (0.8, 0.8, 0.8, 1.0)

for i in range(n):
  v = dataDict4c.get(i*BIN_SIZE)
  
  if v is None:
    data4c[:,i] = (0.0, 0.0, 0.1, 1.0)

  else:
    data4c[:,i] = (0.0, v, v, 1.0)

#

data = empty((m , n, 4), float)
data[:,:] = (0.5, 0.5, 0.5, 1.0)

for f, fileName in enumerate(fileNames):
  p0 = f * 62
  p3 = p0 + 50
  data[p3+1:p3+11,:] = data4c
  
  d = distDictCell[fileName].values()
  dMin = min(d)
  dMax = max(d)
  

  for i in range(n):
    v = distDictCell[fileName].get(i*BIN_SIZE)
    
    if v is None:
      data[p0:p3,i] = (0.0, 0.0, 0.15, 1.0)
    
    else:
      v = 1.0 - (v - dMin)/dMax
      v = v * v
      data[p0:p3,i] = (v, v*0.67, 0.0, 1.0)

  p1 = p0+15
  p2 = p0+36
  for i in range(n):
    v = cisDictCell[fileName].get(i*BIN_SIZE, 0)
    if v:
      data[p1:p2,i] = (0.0, 0.1, 1.0, 1.0)
  
p0 = nCells * 62
p3 = p0 + 50
data[p3+1:p3+11,:] = data4c

d = dAveDict.values()
dMin = min(d)
dMax = max(d)

for i in range(n):
  v = dAveDict.get(i*BIN_SIZE)
  
  if v is None:
    data[p0:p3,i] = (0.0, 0.0, 0.1, 1.0)

  else:
    v = 1.0 - (v - dMin)/dMax
    #v = v * v
    data[p0:p3,i] = (v, v*0.67, 0.0, 1.0)

p1 = p0+15
p2 = p0+36
for i in range(n):
  v = cisDict.get(i*BIN_SIZE, 0)
  if v:
    data[p1:p2,i] = (0.0, 0.1, 1.0, 1.0)

pixmapArray = array(255 * data, uint8, order='C')

scale = [x*10 for x in range(18)]
plot = PixmapPlot(popup, pixmapArray, scale, markers)

app.start()
