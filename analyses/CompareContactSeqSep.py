import sys, numpy

sys.path.append('/home/tjs23/nuc/')

from NucApi import Nucleus

from gui.qtgui.Application import Application
from gui.qtgui.BasePopup import BasePopup
from gui.qtgui.Graph import GraphAxis, GraphDataSet, Graph

fileName1 = '../GrowingAll_Chr18_40K.nuc'
fileName2 = '../SenescentAll_Chr18_40K.nuc'

fLadFile = '/home/tjs23/chromoVista/data/senescence/layers/fLADs.dat'
cLadFile = '/home/tjs23/chromoVista/data/senescence/layers/cLADs.dat'

nuc1 = Nucleus(fileName1)
nuc2 = Nucleus(fileName2)

nuc1.saveAs(fileName1 + '.temp')
nuc2.saveAs(fileName2 + '.temp')

nuc1.importGenomeData(fLadFile, 'fLAD')
nuc2.importGenomeData(fLadFile, 'fLAD')
nuc1.importGenomeData(cLadFile, 'cLAD')
nuc2.importGenomeData(cLadFile, 'cLAD')

app = Application()

popup = BasePopup(title='Graph Popup')
popup.setSize(600, 400)


dataSets = []
dataSetsC = []
dataSetsF = []
dataSetsCF = []

def getHistDatset(values, numBins=50):
  
  hist, edges = numpy.histogram(values, 50)
  n = len(hist)
  xyValues = [((edges[i]+edges[i+1])/2.0, hist[i]) for i in range(n)]
  
  return xyValues

data = [(nuc1, 'Growing'  ,('#800000','#800000','#008000','#000080'), 'circle'),
        (nuc2, 'Senescent',('#000080','#FF8080','#80FF80','#8080FF'), 'star')]

for nuc, prefix, colors, symbol in data:
  
  color1, color2, color3, color4, = colors
  contactDict = nuc.getContacts(cis=True, trans=False, singleCell=True)
  
  for key in contactDict:
    chromo = key[0]
    name = prefix + ' Chr' + chromo[:-1]
    
    fGenomeData = nuc.getGenomeData('fLAD', chromosomes=[chromo])[chromo]
    cGenomeData = nuc.getGenomeData('cLAD', chromosomes=[chromo])[chromo]
 
    fRegions = fGenomeData[0].tolist()
    cRegions = cGenomeData[0].tolist()
 
    contacts = contactDict[key]
 
    points = [(x[0], x[1]) for x in contacts]
    
    sepsF  = []
    sepsC  = []
    sepsCF = []
 
    for p1, p2 in points:
      
      
      isF1 = False
      isF2 = False
      isC1 = False
      isC2 = False
      
      for f1, f2 in fRegions:
        if f1 < p1 < f2:
          isF1 = True
          break
      
      for f1, f2 in fRegions:
        if f1 < p2 < f2:
          isF2 = True
          break
      
      for c1, c2 in cRegions:
        if c1 < p1 < c2:
          isC1 = True
          break
      
      for c1, c2 in cRegions:
        if c1 < p2 < c2:
          isC2 = True
          break
      
      sep = abs(p2-p1) / 1e6
      
      if isF1 and isF2:
        sepsF.append(sep)
      
      elif isC1 and isC2:
        sepsC.append(sep)
        
      elif isC1 and isF2:
        sepsCF.append(sep)
      
      elif isF1 and isC2:
        sepsCF.append(sep)
      
    pointsA = contacts[:,0]
    pointsB = contacts[:,1]
    
    seps   = abs(pointsA - pointsB) / 1e6
    
    xyValues = getHistDatset(seps)
    xyValuesF = getHistDatset(sepsF)
    xyValuesC = getHistDatset(sepsC)
    xyValuesCF = getHistDatset(sepsCF)
 
    ds   = GraphDataSet(xyValues,   name, color1, plotType='line', symbol='circle', symbolSize=3)
    dsF  = GraphDataSet(xyValuesF,  prefix + ' Fac', color2, plotType='line', symbol=symbol, symbolSize=3)
    dsC  = GraphDataSet(xyValuesC,  prefix + ' Con', color3, plotType='line', symbol=symbol, symbolSize=3)
    dsCF = GraphDataSet(xyValuesCF, prefix + ' Con/Fac', color4, plotType='line', symbol=symbol, symbolSize=3)

    dataSets.append(ds)
    dataSetsF.append(dsF)
    dataSetsC.append(dsC)
    dataSetsCF.append(dsCF)


dataSetsB = dataSetsC + dataSetsF + dataSetsCF

xAxis = GraphAxis('Seq. sepratation (Mb)', labels=None, ticks=True)
yAxis = GraphAxis('Count', labels=None, ticks=True,)
graph = Graph(popup, (xAxis, yAxis), dataSets, title='Contact sequence separation', grid=(0,0))
  
xAxis = GraphAxis('Seq. sepratation (Mb)', labels=None, ticks=True)
yAxis = GraphAxis('Count', labels=None, ticks=True,)
graph = Graph(popup, (xAxis, yAxis), dataSetsB, title='Contact sequence separation', grid=(0,1))  


app.start()
