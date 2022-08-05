import sys
from os.path import dirname, abspath
from math import log

thisDir = dirname(abspath(__file__))
sys.path.remove(thisDir)

nucDir = dirname(thisDir)
sys.path.append(nucDir)

from gui.qtgui.BasePopup import BasePopup
from gui.qtgui.Graph import GraphAxis, GraphDataSet, Graph
from cUtil.samread import pairedSamSeqSepHist

def graphMeanContactSeqSep(dataSets, bamFileName, color):

  step = 0.1
  
  print "Reading BAM file"
  sepDict = pairedSamSeqSepHist(bamFileName, 'rb')
  
  chromosomes = sorted(sepDict.keys())
  yVals = {}

  for chromo in chromosomes:
    if chromo in ('Y','MT'):
      continue
    
    hist = sepDict[chromo]
    n = float(hist.sum())
    
    if not n:
      continue

    for bin, val in enumerate(hist):
      if not val:
        continue
    
      x = (bin*step) + step/2.0
      #x = log(bin*1e4+1.0 , 10)
      y = log(val/n, 10.0)
      
      if x in yVals:
        yVals[x].append(y)
      else:
        yVals[x] = [y]
  
    
  xVals = sorted(yVals.keys())
  yVals = [array(yVals[x]).mean() for x in xVals]
  
  
  xyValues = zip(xVals, yVals)
  
  ds = GraphDataSet(xyValues, chromo, color, plotType='line', symbol='circle', symbolSize=1)
  dataSets.append(ds)
 


def graphContactSeqSep(dataSets, bamFileName, color):

  step = 0.1
  
  print "Reading BAM file"
  sepDict = pairedSamSeqSepHist(bamFileName, 'rb')
  
  chromosomes = sorted(sepDict.keys())

  for chromo in chromosomes:
    if chromo in ('Y','MT'):
      continue
    
    hist = sepDict[chromo]
    n = float(hist.sum())
    
    if not n:
      continue
    
    xyValues = []
    for bin, val in enumerate(hist):
      if not val:
        continue
    
      x = (bin*step) + step/2.0
      #x = log(bin*1e4+1.0 , 10)
      y = log(val/n, 10.0)
      xyValues.append((x,y))
        
    if not xyValues:
      continue
    
    ds = GraphDataSet(xyValues, chromo, color, plotType='line', symbol='circle', symbolSize=1)
    dataSets.append(ds)
 

if __name__ == '__main__':
  
  from PySide import QtCore, QtGui
 
  qtApp = QtGui.QApplication(('Graph',))
  
  popup = BasePopup(title='Graph popup')  
  xAxis = GraphAxis('Seq. sepratation Log10(bp) binned',
                    labels=None, ticks=True)
  yAxis = GraphAxis('Log10 proportion',
                    labels=None, ticks=True)

  fileNameP = '/data/sample_892_L007_ESC6_Stefan_HiC_sort.bam'
  fileNameS = '/home/tjs23/nuc/data/SiCUP/sample_1527/merged_processing/sample_1527_L001_L002_merged_CTCTCTAC_CTAAGCCT_SiCUPPED.bam'
  
  dataSets = []
  
  graphContactSeqSep(dataSets, fileNameP, '#008000')
  graphContactSeqSep(dataSets, fileNameS, '#FF0000')

  graph = Graph(popup, (xAxis, yAxis), dataSets, title='Contact sequence separation')  
  
  sys.exit(qtApp.exec_())
