from PySide import QtCore, QtGui

from gui.qtgui.ButtonList import ButtonList
from gui.qtgui.Frame import Frame
from gui.qtgui.LabelFrame import LabelFrame
from gui.qtgui.Label import Label
from gui.qtgui.MessageDialog import showOkCancel
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.Table import ObjectTable, Column
from gui.qtgui.Colors import ColorDialog

from numpy import array

Qt = QtCore.Qt

class LabelsPanel(Frame):
  
  def __init__(self, parent, mainApp):
  
    Frame.__init__(self, parent)
    
    self.mainApp = mainApp
    self.labelData = []
    self.isModified = False
    
    frame = LabelFrame(self, 'Positional Labels', grid=(0,0))
    
    row = 0
    
    texts = ['Add new', 'Remove selected']
    callbacks = [self.addLabel, self.removeLabels] 
    buttons = ButtonList(frame, texts, callbacks, grid=(row,0))
               
    row += 1

    texts = ['Label chr. starts', 'Label chr. ends']
    callbacks = [self.addStarts, self.addEnds] 
    buttons = ButtonList(frame, texts, callbacks, grid=(row,0))

    row += 1

    texts = ['Set same colour', 'Set same label']
    callbacks = [self.sameColor, self.sameLabel] 
    buttons = ButtonList(frame, texts, callbacks, grid=(row,0))
               
    row += 1
    
    columns = [Column('Chr.', self.getChromosome,
                      getEditValue=self.getEditChromosome,
                      setEditValue=self.setChromosome),
               Column('Position\n(bp)', self.getPosition,
                      setEditValue=self.setPosition, format='{:,}'),
               Column('Label', self.getLabelText,
                      setEditValue=self.setLabelText),
               Column('Color', None, 
                      getEditValue=self.getColor,
                      getColor=self.getColor,
                      setEditValue=self.setColor,
                      editClass=ColorDialog)]
    self.table = ObjectTable(frame, columns, None,
                             callback=None, sortCol=None,
                             multiSelect=True, grid=(row,0))

 
  def getEditChromosome(self, obj):
  
    i, chromo, datum = obj
    
    chromos = self.mainApp.nuc.getChromosomes()
    
    
    if chromo in chromos:
      idx = chromos.index(chromo)
    else:
      idx = 0
    
    return chromos, chromos, idx  
  
  
  def getChromosome(self, obj):
    
    nuc = self.mainApp.nuc
    i, chromo, datum = obj
    
    return chromo
    
  
  def setChromosome(self, obj, val):
  
    nuc = self.mainApp.nuc
    i, chromo, datum = obj
    
    self.isModified = True
    self.labelData[i] = (i, val, datum)
    self.updateTable()
    
  
  def getPosition(self, obj):
    
    nuc = self.mainApp.nuc
    i, chromo, datum = obj
    
    return int(datum[1])
    
  
  def setPosition(self, obj, val):
    
    nuc = self.mainApp.nuc
    i, chromo, datum = obj
    
    datum = list(datum)
    datum[1] = val
     
    
    self.isModified = True
    self.labelData[i] = (i, chromo, tuple(datum))
    self.updateTable()
  
  
  def getLabelText(self, obj):
  
    nuc = self.mainApp.nuc
    i, chromo, datum = obj
    
    return datum[0]
  
  
  def setLabelText(self, obj, val):
  
    nuc = self.mainApp.nuc
    i, chromo, datum = obj
    
    datum = list(datum)
    datum[0] = val

    self.isModified = True
    self.labelData[i] = (i, chromo, tuple(datum))
    self.updateTable()
    

  def getColor(self, obj):
    
    nuc = self.mainApp.nuc
    i, chromo, datum = obj
    
    rgb = [int(x*255) for x in datum[2:5]]
    
    return QtGui.QColor(*rgb)
    
    
  def setColor(self, obj, colorObj):
    
    nuc = self.mainApp.nuc
     
    if colorObj: # not cancelled
      i, chromo, datum = obj
      datum = list(datum)
      datum[2:5] = colorObj.getRgbF()[:3]
      self.isModified = True
      self.labelData[i] = (i, chromo, tuple(datum))
      self.updateTable()
  
  
  def leaveEvent(self, event):
    
    if self.isModified:
      self.saveLabels()
  
  
  def showEvent(self, event):
  
    self.loadLabels()
    self.updateTable()
  
    QtGui.QWidget.showEvent(self, event)
 
 
  def saveLabels(self):
    
    nuc = self.mainApp.nuc
    dataDict = {}
    for i, chromo, datum in self.labelData:
      if chromo in dataDict:
        dataDict[chromo].append(datum)
  
      else:
        dataDict[chromo] = [datum,]
    
    for chromo in nuc.getChromosomes():
      nuc.setChromoLabels(chromo, dataDict.get(chromo, []))

    self.isModified = False
   
  
  def loadLabels(self):
    
    self.isModified = False
    self.labelData = []
    i = 0
    nuc = self.mainApp.nuc
    for chromo in nuc.chromosomes:
      data = nuc.getChromoLabels(chromo)
      
      if data is not None:
        for datum in data:
          self.labelData.append((i, chromo, tuple(datum)))
          i += 1
        
  
  def addLabel(self):
  
    nuc = self.mainApp.nuc
    chromos = nuc.getChromosomes()
    
    if chromos:
      chromo = chromos[0]
      i = len(self.labelData)
      datum = (i, chromo, ('Label', 0, 0.0, 0.5, 1.0, 0, 0 ,0)) # Pos, r,g,b, opts
      self.isModified = True
      self.labelData.append(datum)
      self.updateTable()


  def addStarts(self):
  
    nuc = self.mainApp.nuc
    chromos = nuc.getChromosomes()
    particGroup = nuc._getParticleGroup()
    
    for chromo in chromos:
      start = particGroup[chromo]['positions'][0]
      r, g, b = nuc.getChromoColor(chromo).rgb
      
      i = len(self.labelData)
      datum = (i, chromo, ('Chr %s start' % chromo, start, r, g, b, 0, 0 ,0)) # Pos, r,g,b, opts
      self.isModified = True
      self.labelData.append(datum)
      self.updateTable()


  def addEnds(self):
  
    nuc = self.mainApp.nuc
    chromos = nuc.getChromosomes()
    particGroup = nuc._getParticleGroup()
    
    for chromo in chromos:
      end = particGroup[chromo]['positions'][-1]
      r, g, b = nuc.getChromoColor(chromo).rgb
      
      i = len(self.labelData)
      datum = (i, chromo, ('Chr %s end' % chromo, end-1, r, g, b, 0, 0 ,0)) # Pos, r,g,b, opts
      self.isModified = True
      self.labelData.append(datum)
      self.updateTable()
  
  
  def sameLabel(self):
  
    nuc = self.mainApp.nuc
    objs = self.table.getSelectedObjects()
    
    if not objs:
      return
      
    label = objs[0][2][0]
    
    for i, chromo, datum in objs[1:]:
      datum = list(datum)
      datum[0] = label
      self.labelData[i] = (i, chromo, tuple(datum))
    
    self.isModified = True
    self.updateTable()
  
  
  def sameColor(self):
  
    nuc = self.mainApp.nuc
    objs =self.table.getSelectedObjects()
    
    if not objs:
      return
      
    r, g, b = objs[0][2][2:5]
    
    for i, chromo, datum in objs[1:]:
      datum = list(datum)
      datum[2:5] = [r,g,b]
      self.labelData[i] = (i, chromo, tuple(datum))
    
    self.isModified = True
    self.updateTable()   
    
  
  def removeLabels(self):
  
    nuc = self.mainApp.nuc
    objs = set(self.table.getSelectedObjects())
    
    if not objs:
      return
    
    msg = 'Remove selected labels?'
    if not showOkCancel('Confirm', msg, parent=self):
      return
    
    self.isModified = True
    self.labelData = [x for x in self.labelData if x not in objs]
    
    self.updateTable()
      

  def updateTable(self):
  
    
    if self.isVisible():
      self.table.setObjects(self.labelData)
      self.table.update()
    
    self.saveLabels()
    self.mainApp.updateContents()
    
    
  def updateContents(self):
  
    if self.isVisible():
      self.loadLabels()
      self.table.setObjects(self.labelData)
      self.table.update()
    
