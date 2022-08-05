import sys, os, numpy

from math import log
from random import shuffle

from os.path import abspath, dirname
sys.path.append(abspath(dirname(dirname(__file__))))

from gui.qtgui.Graph import GraphAxis, GraphDataSet, Graph
from gui.qtgui.BasePopup import BasePopup
from gui.qtgui.Application import Application

DEFAULT_COLORS = ['#F00000','#F0F000','#00F000',
                  '#00F0F0','#0000F0','#F000F0',
                  '#800000','#008000','#0000A0']

from NucApi import Nucleus

# 2D trans distance

# Cis vs trans p-Value
#   Get total trans for each - calc p(trans|chromo) - fraction of total
#   Random p(trans|chrA, chrB) = p(trans|chrA) * p(trans|chrB)
#   Compare with observed



dirName = '/data/hi-c/sc2' # sys.argv[1]

fileNames = os.listdir(dirName)

fileNames = [os.path.join(dirName, f) for f in fileNames if f.endswith('.nuc')]

fileNames.sort()

dataSets = []

header = ['Cell','Contacts','Cis<10 Kb','Cis','Trans','%Cis<10 Kb', '%Trans','%Isolated','%Trans Iso.','G-score']
header = ['%16.16s' % header[0]] + ['%10.10s' % x for x in header[1:]]

h = ' '.join(header)
print(h)
print(('-' * len(h)))

rVals = numpy.array([])

for f, filePath in enumerate(fileNames):
    
  dirName, fileName = os.path.split(filePath)
  libName = os.path.splitext(fileName)[0]

  nuc = Nucleus(filePath, 'r')
  
  groupName = 'singleCell'
  
  n = nuc.getNumContacts(groupName)
  cc1 = nuc.getNumContacts(groupName, trans=False, maxSep=1e4)
  c = nuc.getNumContacts(groupName, trans=False)
  t = nuc.getNumContacts(groupName, cis=False)
  pcc = 100.0 * nuc.getNumContacts(groupName, trans=False, maxSep=1e4, asFraction=True)
  pt = 100.0 * nuc.getNumContacts(groupName, cis=False, asFraction=True)
  pi = 100.0 * nuc.getNumIsolatedContacts(groupName, threshold=int(1e6), asFraction=True)
  pit = 100.0 * nuc.getNumIsolatedContacts(groupName, threshold=int(1e6), cis=False, asFraction=True)
  
  contactDictTrans = nuc.getContacts(groupName, cis=False)
  contactDictAll = nuc.getContacts(groupName)
  
  transDict = {}
  cisTransDict = {}
  
  for pairKey in contactDictAll:
    chrA, chrB = pairKey
    contacts = contactDictAll[pairKey]
    
    posA = contacts[0]
    posB = contacts[1]
  
    if chrA in cisTransDict:
      cisTransDict[chrA].update(posA)
    else:
      cisTransDict[chrA] = set(posA)

    if chrB in cisTransDict:
      cisTransDict[chrB].update(posB)
    else:
      cisTransDict[chrB] = set(posB)
  
  for pairKey in contactDictTrans:
    chrA, chrB = pairKey
    contacts = contactDictTrans[pairKey]
  
    posA = contacts[0]
    posB = contacts[1]
    
    if chrA in transDict:
      transDict[chrA].update(posA)
    else:
      transDict[chrA] = set(posA)

    if chrB in transDict:
      transDict[chrB].update(posB)
    else:
      transDict[chrB] = set(posB)
    
  g = nuc.calcTransContactGscore()
    
  print('%16.16s %10d %10d %10d %10d %10.2f %10.2f %10.2f %10.2f %10.3f' % (libName, n, cc1, c, t, pcc, pt, pi, pit, g))
  
  vals = numpy.array([])
  valsAll = numpy.array([])
  
  for chromo in cisTransDict:
    nSel = len(transDict.get(chromo, []))
    if not nSel:
      continue
    
    pos = list(cisTransDict[chromo])
    
    for i in range(10):
      shuffle(pos)
      posSel = pos[:nSel]
      posSel = numpy.array(sorted(posSel), int)
      deltas = posSel[1:] - posSel[:-1]
      valsAll = numpy.append(valsAll, deltas)
    
  
  for chromo in transDict:
    pos = numpy.array(sorted(transDict[chromo]), int)
    deltas = pos[1:] - pos[:-1]
    vals = numpy.append(vals, deltas)
    
    start, end = nuc.getChromosomeLimits(chromo)
    rPos = numpy.random.uniform(start, end, len(pos)).tolist()
    rPos = numpy.array(sorted(rPos))
    
    rDeltas = rPos[1:] - rPos[:-1]
    rVals = numpy.append(rVals, rDeltas)
        
  vals = numpy.array(vals, float) / 1000.0
  vals = numpy.log(vals)
  hist, edges = numpy.histogram(vals, 34, (-7.0, 10.0))
  hist = numpy.array(hist, float) / hist.sum()
  
  """
  xyValues = [((edges[i]+edges[i+1])/2.0, hist[i]) for i in range(len(hist))]
  ds = GraphDataSet(xyValues, libName + ' Trans', DEFAULT_COLORS[f % len(DEFAULT_COLORS)],
                    plotType='line', symbol='circle', symbolSize=2)
  dataSets.append(ds)
  """
        
  valsAll = numpy.array(valsAll, float) / 1000.0
  valsAll = numpy.log(valsAll)
  hist2, edges = numpy.histogram(valsAll, 34, (-7.0, 10.0))
  hist2 = numpy.array(hist2, float) / hist2.sum()
  
  """
  xyValues = [((edges[i]+edges[i+1])/2.0, hist2[i]) for i in range(len(hist2))]
  ds = GraphDataSet(xyValues, libName+' All', DEFAULT_COLORS[(f+3) % len(DEFAULT_COLORS)],
                    plotType='line', symbol='circle', symbolSize=2)
  dataSets.append(ds)
  """
  
  idx = hist.nonzero()
  hist = hist[idx]
  hist2 = hist2[idx]
  
  idx = hist2.nonzero()
  hist = hist[idx]
  hist2 = hist2[idx]
  
  histG = hist * numpy.log2(hist/hist2)

  
  xyValues = [((edges[i]+edges[i+1])/2.0, histG[i]) for i in range(len(histG))]
  ds = GraphDataSet(xyValues, libName+' G', DEFAULT_COLORS[(f) % len(DEFAULT_COLORS)],
                    plotType='line', symbol='circle', symbolSize=2)
  dataSets.append(ds)
  
  
"""
rVals = numpy.array(rVals, float) / 1000.0 
rVals = numpy.log(rVals)
hist, edges = numpy.histogram(rVals, 40)
hist = numpy.array(hist, float) / hist.sum()
xyValues = [((edges[i]+edges[i+1])/2.0, hist[i]) for i in range(len(hist))]

ds = GraphDataSet(xyValues, 'Random', '#000000', plotType='line', symbol='circle', symbolSize=2)
dataSets.append(ds)
"""

app = Application()
  
popup = BasePopup(title='Popup graph')
xAxis = GraphAxis('Log Separation (kb)', labels=None, ticks=True)
yAxis = GraphAxis('Count', labels=None, ticks=True,)
graph = Graph(popup, (xAxis, yAxis), dataSets, title='Contact sequence spacing') 

app.start()
