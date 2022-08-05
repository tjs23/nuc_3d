from numpy import array, histogram, int32
from math import ceil
from PySide import QtCore, QtGui
Qt = QtCore.Qt

from gui.qtgui.ButtonList import ButtonList
from gui.qtgui.CheckButton import CheckButton
from gui.qtgui.Colors import ColorDialog
from gui.qtgui.SpinBox import FloatSpinBox
from gui.qtgui.LabelFrame import LabelFrame
from gui.qtgui.Label import Label
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.MessageDialog import showOkCancel
from gui.qtgui.Table import ObjectTable, Column
from gui.qtgui.Graph import Graph, GraphAxis, GraphDataSet

from cUtil import dataLayer
from NucApi import DATA_TRACK_SYMBOLS, DERIVED, INNATE, EXTERNAL

SOURCE_NAME_DICT = {DERIVED:'Derived', INNATE:'Genome', EXTERNAL:'Experimental'}

class DataTrackPanel(QtGui.QWidget):

  def __init__(self, mainApp, parent=None):
  
    QtGui.QWidget.__init__(self, parent)
    
    self.mainApp = mainApp
    self.chromosome = None
    
    row = 0
    frame = LabelFrame(self, 'Available datasets', grid=(row,0), stretch=(2,1))
    
    columns = [
               Column('Name', self.getDataName),
               Column('Class', self.getDataSource),
               Column('Colour', None,
                      getColor=self.getDataColor,
                      getEditValue=self.getDataColor,
                      setEditValue=self.setDataColor,
                      editClass=ColorDialog),
               Column('Show?', self.getDataShown,
                      setEditValue=self.setDataShown),
               Column('Show\ntext?', self.getDataShowText,
                      setEditValue=self.setDataShowText),
               Column('Min. value', self.getDataValMin,
                      setEditValue=self.setDataValMin,
                      editStep=0.1),
               Column('Max. value', self.getDataValMax,
                      setEditValue=self.setDataValMax,
                      editStep=0.1),
               Column('Scale', self.getDataScale,
                      setEditValue=self.setDataScale,
                      editStep=0.1),
               Column('Symbol', self.getDataShape,
                      getEditValue=self.getEditDataShape,
                      setEditValue=self.setDataShape),
               Column('Reference file', self.getDataRef),
               Column('Num.\nChrs.', self.getDataNumChrs),
               Column('Num.\nValues', self.getDataNumValues),
               ]
    
    self.data = None
    self.dataTable = ObjectTable(frame, columns, [], callback=self.selectData,
                                 multiSelect=True, grid=(0,0))
    
    texts = ['Remove selected', 'Import dataset', 'Export dataset']
    callbacks = [self.removeData,
                 self.mainApp.importDataTrack,
                 self.exportData]
    buttons = ButtonList(frame, texts, callbacks, grid=(1,0))
    
    texts = ['Move to experiment\nreference file', 'Move to genome\nreference file', 'Move to local file']
    callbacks = [self.moveToExperimentRef,
                 self.moveToGenomeRef,
                 self.moveToLocal]
    buttons = ButtonList(frame, texts, callbacks, grid=(2,0))
    
    frame = LabelFrame(self, 'Value distribution', grid=(row,1), stretch=(2,1))
    
    xAxis = GraphAxis('Value', labels=None, ticks=True)
    yAxis = GraphAxis('Count', labels=None, ticks=True, valRange=(0.0, 1.0))
  
    axes = (xAxis, yAxis)
  
    self.valGraph = Graph(frame, axes, [], title='Value distribution',
                          size=(300, 200), callback=None,
                          fgColor=QtGui.QColor(255, 255, 255, 255),
                          bgColor=QtGui.QColor(0, 0, 0, 255),
                          mgColor=QtGui.QColor(64, 64, 64, 255),
                          grid=(0,0))  
   
    
    row += 1
    frame = LabelFrame(self, 'Normalisation', grid=(row,0), stretch=(0,1), gridSpan=(1,2))
    
    
    self.clipValsCheck = CheckButton(frame, text='Clip min/max values', selected=False, grid=(0,0))
    
    texts = ['Quantile', 'Max', 'Sqrt Max', 'Log Max', 
             'Unit', 'Inverse', 'Reset']
    callbacks = [self.normaliseQuant,
                 self.normaliseMax,
                 self.normaliseSqrt,
                 self.normaliseLog,
                 self.normaliseUnit,
                 self.normaliseInv,
                 self.normaliseOrig]
    buttons = ButtonList(frame, texts, callbacks, grid=(1,0), gridSpan=(1,2))
    
    row += 1
    frame = LabelFrame(self, 'Chromosomal distribution', grid=(row,0), stretch=(0,1), gridSpan=(1,2))
    self.layout().setRowMinimumHeight(row, 200)
    
    label = Label(frame, ' Chromosome: ', grid=(0,0))
    
    self.chromoPulldown = PulldownList(frame, callback=self.selectChromosome, grid=(0,1))
    
    label = Label(frame, ' Bin size: ', grid=(0,2))
    
    self.binSizeEntry = FloatSpinBox(frame, value=1.0, callback=self._redrawGraphs,
                                     minValue=1e-5, maxValue=10.0, suffix=' Mb',
                                     multiplier=2.0, grid=(0,3))
    self.binSizeEntry.setDecimals(4)                                         
    
    xAxis = GraphAxis('Position (Mb)', labels=None, ticks=True)
    yAxis = GraphAxis('% Chromosome Max.', labels=None, ticks=True, valRange=(0.0, 100.0))
  
    axes = (xAxis, yAxis)
    self.yAxis = yAxis
  
    self.posGraph = Graph(frame, axes, [], title='Data track value distribution',
                          size=(700, 150), callback=None,
                          fgColor=QtGui.QColor(255, 255, 255, 255),
                          bgColor=QtGui.QColor(0, 0, 0, 255),
                          mgColor=QtGui.QColor(64, 64, 64, 255),
                          grid=(1,0), gridSpan=(1,5))  
  
    frame.layout().setColumnStretch(4, 2)
  
  
  def selectChromosome(self, chromo):
    
    if chromo != self.chromosome:
      start, end = self.mainApp.nuc.getChromosomeLimits(chromo)
      
      binSize = self.binSizeEntry.get()
      idealSize = (end-start)/1e8
      
      if binSize > idealSize:
        self.binSizeEntry.set(idealSize)
      
      self.chromosome = chromo
      self._redrawGraphs()
  
  
  def updateChromosomes(self):
    
    chromo = self.chromosome
    nuc = self.mainApp.nuc
    
    if nuc and nuc.getChromosomes():
      chromos = nuc.getChromosomes()
      
      if chromo in chromos:
        index = chromos.index(chromo)
      
      else:
        index = 0
        chromo = chromos[0]
      
    else:
      chromos = []
      chromo = None
      index = 0 
    
    self.chromoPulldown.setData(chromos, chromos, index)
    self.selectChromosome(chromo)
 
  
  def resizeEvent(self, event):
  
    w = event.size().width() - 150
    
    prev, h= self.posGraph.size
    
    self.posGraph.setSize(w, h)
    
    QtGui.QWidget.resizeEvent(self, event)
  
  
  def _normaliseData(self, method):
    
    nuc = self.mainApp.nuc
    selected = self.dataTable.getSelectedObjects()
    clipValues = self.clipValsCheck.get()
    
    if selected:
      for typ, code in selected:
        if typ in nuc.dataTracks:
          if code in nuc.dataTracks[typ]:
            
            if clipValues:
              minMaxVals = nuc.getDataTrackThresholds(typ, code)
            else:
              minMaxVals = None
             
            nuc.calcNormDataTrack(code, typ, None, method, minMaxVals)
                
      self.mainApp.updateContents()
  
  
  def normaliseLog(self):
   
    self._normaliseData('log')


  def normaliseInv(self):
   
    self._normaliseData('inverse')
    
    
  def normaliseOrig(self):
   
    self._normaliseData('orig')
    
 
  def normaliseMax(self):
   
    self._normaliseData('max')
    
  def normaliseQuant(self):
   
    self._normaliseData('quantile')
    

  def normaliseSqrt(self):
   
    self._normaliseData('sqrt')


  def normaliseUnit(self):
  
    self._normaliseData('unity')
    
  
  def removeData(self):
    
    nuc = self.mainApp.nuc
    selected = self.dataTable.getSelectedObjects()
    
    if selected:
      nSel = len(selected)
      
      if nSel > 1:
        names = ', '.join([s[1] for s in selected])
        msg = 'Remove %d datasets?\n(%s)' % (nSel,names)
      else:
        msg = 'Remove "%s" dataset?' % (selected[0][1],) 
      
      if not showOkCancel('Query', msg, parent=self):
        return
        
      for typ, code in selected:
        nuc.removeDataTrack(typ, code)
            
      self.mainApp.updateContents()
  
  
  def exportData(self):
    
    selected = self.dataTable.getCurrentObject()
    
    if selected:
      self.mainApp.exportDataTrack(selected)
      
   
  def _isDataLocal(self, typ, code):
     
    nuc = self.mainApp.nuc.getRefDataTrackNuc(typ, code)
    
    return nuc is self.mainApp.nuc
  
  
  def moveToExperimentRef(self):
  
    nuc = self.mainApp.nuc
    selected = self.dataTable.getSelectedObjects()
    
    selected = [obj for obj in selected if obj[0] == EXTERNAL]
    selected = [obj for obj in selected if self._isDataLocal(*obj)]
    
    if selected:
      nSel = len(selected)
      
      if nSel > 1:
        names = ', '.join([s[1] for s in selected])
        msg = 'Move %d local datasets to experiment reference file?\n(%s)' % (nSel,names)
      else:
        msg = 'Move local dataset "%s" to experiment reference file?' % (selected[0][1],) 
      
      if not showOkCancel('Confirm', msg, parent=self):
        return
      
      rNuc = nuc._experimentRefNuc
      
      if not rNuc:
        msg = 'Experiment reference file not set. Set now?'
        if showOkCancel('Confirm', msg, parent=self):
          self.mainApp.setExperimentRef()
        else:
          return
          
      rNuc = nuc.getExperimentRefNuc(mode='a')
      if not rNuc:
        return
        
      for typ, code in selected:
        rNuc.importNucDataTrack(nuc, typ, code)
      
      for typ, code in selected:
        nuc.removeDataTrack(typ, code)
      
      nuc.updateProxyDataTracks(rNuc, EXTERNAL)
      rNuc.save()
            
      self.mainApp.updateContents()
   
      
  def moveToGenomeRef(self):
  
    nuc = self.mainApp.nuc
    selected = self.dataTable.getSelectedObjects()

    selected = [obj for obj in selected if obj[0] == INNATE]
    selected = [obj for obj in selected if self._isDataLocal(*obj)]
    
    if selected:
      nSel = len(selected)
      
      if nSel > 1:
        names = ', '.join([s[1] for s in selected])
        msg = 'Move %d local datasets to genome reference file?\n(%s)' % (nSel,names)
      else:
        msg = 'Move local dataset "%s" to genome reference file?' % (selected[0][1],) 
      
      if not showOkCancel('Confirm', msg, parent=self):
        return
             
      rNuc = nuc._genomeRefNuc
      
      if not rNuc:
        msg = 'Genome reference file not set. Set now?'
        if showOkCancel('Confirm', msg, parent=self):
          self.mainApp.setGenomeRef()
        else:
          return
          
      rNuc = nuc.getGenomeRefNuc(mode='a')
      if not rNuc:
        return
        
      for typ, code in selected:
        rNuc.importNucDataTrack(nuc, typ, code)
      
      for typ, code in selected:
        nuc.removeDataTrack(typ, code)
      
      nuc.updateProxyDataTracks(rNuc, INNATE)
      rNuc.save()
      
      self.mainApp.updateContents()

 
  def moveToLocal(self):
 
    nuc = self.mainApp.nuc
    selected = self.dataTable.getSelectedObjects()
    
    selected = [obj for obj in selected if not self._isDataLocal(*obj)]
    
    if selected:
      nSel = len(selected)
      
      if nSel > 1:
        names = ', '.join([s[1] for s in selected])
        msg = 'Copy %d  datasets to local .nuc file?\n(%s)' % (nSel,names)
      else:
        msg = 'Copy dataset "%s" to local .nuc file?' % (selected[0][1],) 
      
      if not showOkCancel('Confirm', msg, parent=self):
        return
              
      for typ, code in selected:
        rNuc = nuc.getRefDataTrackNuc(typ, code)
        nuc.importNucDataTrack(rNuc, typ, code)
      
      nuc.save()
      
      self.mainApp.updateContents()

  
  def showEvent(self, event):
  
    self.updateContents()
            
            
  def updateContents(self):
    
    if self.isVisible():
      nuc = self.mainApp.nuc
 
      codes = []
 
      if nuc:
        for typ in nuc.dataTracks:
          for code in nuc.dataTracks[typ]:
            if 'options' in nuc.dataTracks[typ][code].attrs:
              codes.append((typ, code))
 
        codes.sort()

 
      self.dataTable.setObjects(codes)
      self.updateChromosomes()
      self._redrawGraphs()
  
  
  def _redrawGraphs(self, *args):
  
    nuc = self.mainApp.nuc
    chromo = self.chromosome
    
    binSize = int(self.binSizeEntry.get() * 1e6)
    
    if nuc and chromo:
      dataSets = []
      dataSetsV = []
      xMin = 1e10
      xMax = 0
      vMax = 0
      
      for obj in self.dataTable.getSelectedObjects():
 
        typ, code = obj
        group = nuc.getRefDataTrackGroup(typ, code)

        if chromo in group:
          regions = array(group[chromo]['regions'], int32)
          if not len(regions):
            continue
          
          values = array(group[chromo]['values'], float)[:,1]
          
          hist = dataLayer.regionBinValues(regions, values, int32(binSize))
          hMax = hist.max()
          
          if not hMax:
            continue
          
          hist *= 100.0/hMax

          xyValues = [(i*1e-6, v) for i, v in enumerate(hist)]
          color = self.getDataColor(obj)
          
          ds = GraphDataSet(xyValues, code, color, plotType='histogram', symbolSize=2)
          dataSets.append(ds)
          
          counts, edges = histogram(values, 50)
          vMax = max(vMax, counts.max())
          xyValues = [((edges[i]+edges[i+1])/2.0, counts[i]) for i in range(len(counts))]
          
          ds = GraphDataSet(xyValues, code, color, plotType='histogram', symbolSize=2)
          dataSetsV.append(ds)
      
      xMin = binSize * int(xMin/binSize) *1e-6
      xMax = binSize * int(ceil(xMax/binSize)) *1e-6
      
      if dataSets:
        xAxis = GraphAxis('Position (Mb)', labels=None, ticks=True, valRange=(xMin, xMax))
        axes = (xAxis, self.yAxis) 
        self.posGraph.updateData(dataSets, 'Chromosome %s' % chromo, axes)
        
        xAxis = GraphAxis('Value', labels=None, ticks=True)
        yAxis = GraphAxis('Count', labels=None, ticks=True, valRange=(0, vMax))       
        self.valGraph.updateData(dataSetsV, 'Chromosome %s' % chromo, axes=(xAxis, yAxis))
        
        
  def selectData(self, obj, row, col):
  
    self.data = obj
    
    if obj:
      self._redrawGraphs()
    
  
  def _getDataTrackAttr(self, key, typ, code, index):
    
    nuc = self.mainApp.nuc    
    val = nuc.dataTracks[typ][code].attrs[key][index]
    
    if key == 'options':
      return int(val)
    else:
      return float(val)


  def getDataSource(self, obj):
  
    typ, code = obj
    return SOURCE_NAME_DICT[typ]


  def getDataName(self, obj):
  
    typ, code = obj
    return code
  
  
  def setDataShown(self, obj, isShown):
  
    typ, code = obj
    
    self.mainApp.nuc.setDataTrackDisplayed(typ, code, isShown)
  
    
  def getDataShown(self, obj):
  
    typ, code = obj
    
    if self._getDataTrackAttr('options', typ, code, 0):
      return True
      
    else:
      return False  


  def getDataNumChrs(self, obj):
  
    typ, code = obj

    nuc = self.mainApp.nuc
    group = nuc.getRefDataTrackGroup(typ, code)
    
    if group:
      chromos = [c for c in group]
    
      return len(chromos)
 
 
  def getDataNumValues(self, obj):
  
    typ, code = obj
    
    nuc = self.mainApp.nuc
    group = nuc.getRefDataTrackGroup(typ, code)
    n = 0
    
    if group:
      for chromo in group:
        n += len(group[chromo]['values'])
    
    return n
    
  def setDataColor(self, obj, colorObj):
    
    if colorObj: # not cancelled
      typ, code = obj
      
      if self.mainApp.nuc.setDataTrackColor(typ, code, colorObj.getRgbF()):
        self.mainApp.updateContents()


  def getDataColor(self, obj):
  
    typ, code = obj
    nuc = self.mainApp.nuc
    
    color = nuc.getDataTrackColor(typ, code).rgb24 or (255, 255, 0)
    
    return QtGui.QColor(*color)
  
  
  def setDataShowText(self, obj, value):
    
    typ, code = obj
    
    return self.mainApp.nuc.setDataTrackLabelled(typ, code, value)
    
    
  def getDataShowText(self, obj):
  
    typ, code = obj
    
    if self._getDataTrackAttr('options', typ, code, 1): 
      return True
    
    else:
      return False  
      
      
  def getEditDataShape(self, obj):
  
    typ, code = obj
    index = self._getDataTrackAttr('options', typ, code, 2)
    
    texts = DATA_TRACK_SYMBOLS
    objects = list(range(len(texts)))
    
    if index >= len(texts):
      index = 0
  
    return texts, objects, index
    
    
  def setDataShape(self, obj, value):
  
    typ, code = obj
    
    return self.mainApp.nuc.setDataTrackSymbol(typ, code, value)
    
    
  def getDataRef(self, obj):
  
    typ, code = obj
    nuc = self.mainApp.nuc.getRefDataTrackNuc(typ, code)

    if nuc is self.mainApp.nuc:
      return None
      
    elif nuc:  
      return nuc.root.filename
    
    
  def getDataShape(self, obj):
  
    typ, code = obj
    
    i = self._getDataTrackAttr('options', typ, code, 2)
  
    if i >= len(DATA_TRACK_SYMBOLS):
      i = 0
    
    return DATA_TRACK_SYMBOLS[i]
    
    
  def setDataScale(self, obj, value):
  
    typ, code = obj

    return self.mainApp.nuc.setDataTrackScale(typ, code, value)
    
    
  def getDataScale(self, obj):
  
    typ, code = obj

    return self._getDataTrackAttr('display', typ, code, 4)


  def setDataValMin(self, obj, value):
  
    typ, code = obj
        
    return self.mainApp.nuc.setDataTrackThresholds(typ, code, minVal=value)


  def getDataValMin(self, obj):
  
    typ, code = obj
        
    return self._getDataTrackAttr('display', typ, code, 5)


  def setDataValMax(self, obj, value):
  
    typ, code = obj
        
    return self.mainApp.nuc.setDataTrackThresholds(typ, code, maxVal=value)


  def getDataValMax(self, obj):
  
    typ, code = obj
        
    return self._getDataTrackAttr('display', typ, code, 6)
   
