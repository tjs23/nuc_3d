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

COLOR_OPTS = ['Unaltered', 'Chromosomal', 'User specified', 'Gradient']

class RoiPanel(Frame):
  
  def __init__(self, parent, mainApp):
  
    Frame.__init__(self, parent)
    
    self.mainApp = mainApp
    self.region_data = []
    self.isModified = False
    
    frame = LabelFrame(self, 'Regions of Interest', grid=(0,0))
    
    row = 0
    
    texts = ['Add new', 'Remove selected']
    callbacks = [self.addRegion, self.removeRegions] 
    buttons = ButtonList(frame, texts, callbacks, grid=(row,0))
               
    row += 1

    texts = ['Set same colour', 'Set same width']
    callbacks = [self.sameColor, self.sameSize] 
    buttons = ButtonList(frame, texts, callbacks, grid=(row,0))
               
    row += 1
    
    columns = [Column('Chr.', self.getChromosome,
                      getEditValue=self.getEditChromosome,
                      setEditValue=self.setChromosome),
               Column('Display?', self.getRegionShown,
                      setEditValue=self.setRegionShown),
               Column('Start\n(bp)', self.getStartPosition,
                      setEditValue=self.setStartPosition, format='{:,}'),
               Column('End\n(bp)', self.getEndPosition,
                      setEditValue=self.setEndPosition, format='{:,}'),
               Column('Color mode', self.getColorOpt, 
                      getEditValue=self.getColorOpts,
                      setEditValue=self.setColorOpt),
               Column('User color', None, 
                      getEditValue=self.getColor,
                      getColor=self.getColor,
                      setEditValue=self.setColor,
                      editClass=ColorDialog)]
                      
    self.table = ObjectTable(frame, columns, None,
                             callback=None, sortCol=None,
                             multiSelect=True, grid=(row,0))
  
  def _mod_col(self, i, j, val):
  
     row = list(self.region_data[i])
     row[j] = val
     self.region_data[i] = tuple(row)
    
    
  def getRegionShown(self, obj):
  
    i, chromo, start, end, display_mode, selected, color = obj
    
    return bool(selected)


  def setRegionShown(self, obj, isShown):
   
    isShown = int(isShown)
    i, chromo, start, end, display_mode, selected, color = obj
    
    if isShown != selected:
      self.isModified = True
      self._mod_col(i, 5, isShown)
      self.updateTable()
 
 
  def getEditChromosome(self, obj):
  
    i, chromo, start, end, display_mode, selected, color = obj
    chromos = self.mainApp.nuc.getChromosomes()
    
    if chromo in chromos:
      idx = chromos.index(chromo)
    else:
      idx = 0
    
    return chromos, chromos, idx  
  
  
  def getChromosome(self, obj):
    
    nuc = self.mainApp.nuc
    i, chromo, start, end, display_mode, selected, color = obj
    
    return chromo
    
  
  def setChromosome(self, obj, val):
  
    nuc = self.mainApp.nuc
    i, chromo, start, end, display_mode, selected, color = obj
    
    if val != chromo:
      self.isModified = True
      self._mod_col(i, 1, val)
      self.updateTable()
    
  
  def getStartPosition(self, obj):
    
    nuc = self.mainApp.nuc
    i, chromo, start, end, display_mode, selected, color = obj
    
    return int(start)
    
  
  def setStartPosition(self, obj, val):
    
    nuc = self.mainApp.nuc
    i, chromo, start, end, display_mode, selected, color = obj
    
    if start != val: 
      self.isModified = True
      self._mod_col(i, 2, val)
      self.updateTable()


  def getEndPosition(self, obj):
    
    nuc = self.mainApp.nuc
    i, chromo, start, end, display_mode, selected, color = obj
    
    return int(end)
    
  
  def setEndPosition(self, obj, val):
    
    nuc = self.mainApp.nuc
    i, chromo, start, end, display_mode, selected, color = obj
    
    if end != val:
      self.isModified = True
      self._mod_col(i, 3, val)
      self.updateTable()
  

  def getColorOpts(self, obj):

    i, chromo, start, end, display_mode, selected, color = obj
    
    return COLOR_OPTS, [0, 1, 2, 3], display_mode  


  def getColorOpt(self, obj):
   
    i, chromo, start, end, display_mode, selected, color = obj
    
    
    return COLOR_OPTS[display_mode]


  def setColorOpt(self, obj, val):
    
    i, chromo, start, end, display_mode, selected, color = obj
    
    if  display_mode != val:
      self.isModified = True
      self._mod_col(i, 4, val)
      self.updateTable()


  def getColor(self, obj):
    
    i, chromo, start, end, display_mode, selected, color = obj
    
    color = [int(x) for x in color]
    
    return QtGui.QColor(*color)
    
    
  def setColor(self, obj, colorObj):
    
    nuc = self.mainApp.nuc
     
    if colorObj: # not cancelled
      i, chromo, start, end, display_mode, selected, color = obj

      self.isModified = True
      self._mod_col(i, 6, tuple([255*x for x in colorObj.getRgbF()[:3]]))
      self.updateTable()
  
  
  def leaveEvent(self, event):
    
    if self.isModified:
      self.save_regions()
  
  
  def showEvent(self, event):
  
    self.load_regions()
    self.updateTable()
  
    QtGui.QWidget.showEvent(self, event)
 
 
  def save_regions(self):
    
    nuc = self.mainApp.nuc
    
    dataDict = {}
    for i, chromo, start, end, display_mode, selected, color in self.region_data:
      datum = [start, end, display_mode, 1, selected] + list(color)
      
      if chromo in dataDict:
        dataDict[chromo].append(datum)
  
      else:
        dataDict[chromo] = [datum,]
    
    for chromo in nuc.getChromosomes():
      chromoGroup = nuc.chromosomes[chromo]
      
      rois = dataDict.get(chromo, [])
      nuc._setData('regions', chromoGroup, int, array(rois, int))
      
    self.isModified = False
   
  
  def load_regions(self):
    
    self.isModified = False
    self.region_data = region_data = []
    
    i = 0
    nuc = self.mainApp.nuc
    for chromo in nuc.chromosomes:
      data = nuc.getChromoRegionsOfInterest(chromo)
      
      if data is not None:
        for row in data:
          start, end, display_mode, thick, selected = row[:5]
          color = tuple(row[5:8])
          region_data.append((i, chromo, start, end, display_mode, selected, color))
          i += 1
        
  
  def addRegion(self):
  
    nuc = self.mainApp.nuc
    chromos = nuc.getChromosomes()
    
    if chromos:
      obj = self.table.getCurrentObject()
      
      if obj:
        i, chromo, start, end, display_mode, selected, color = obj
      
      else:
        chromo = chromos[0]
        display_mode = 0
        color = [0, 128, 255]
        
      i = len(self.region_data)
      datum = (i, chromo, 0, 0, display_mode, 0, color)
       
      self.isModified = True
      self.region_data.append(datum)
      self.updateTable()
        
  
  def sameSize(self):
  
    nuc = self.mainApp.nuc
    objs = self.table.getSelectedObjects()
    
    if not objs:
      return
    
    if len(objs) < 2:
      return
    
    i, chromo, start, end, display_mode, color = self.table.getCurrentObject() or objs[0]
    half_width = abs(end-start) / 2
    
    for j, chromo, start, end, display_mode, selected, color in objs:
      if i == j:
        continue
    
      mid = (start + end)/2
      self._mod_col(j, 2, mid - half_width)
      self._mod_col(j, 3, mid + half_width)

    
    self.isModified = True
    self.updateTable()
  
  
  def sameColor(self):
  
    nuc = self.mainApp.nuc
    objs = self.table.getSelectedObjects()
    
    if not objs:
      return
    
    if len(objs) < 2:
      return
    
    i, chromo, start, end, display_mode, color = self.table.getCurrentObject() or objs[0]
    
    for j, chromo, start, end, display_mode, selected, color_old in objs:
      if i == j:
        continue
    
      self._mod_col(j, 6, color)
      
    self.isModified = True
    self.updateTable()

 
  def removeRegions(self):
  
    nuc = self.mainApp.nuc
    objs = set(self.table.getSelectedObjects())
    
    if not objs:
      return
    
    msg = 'Remove selected regions?'
    if not showOkCancel('Confirm', msg, parent=self):
      return
    
    idx = set([x[0] for x in objs])
    
    self.isModified = True
    self.region_data = [x for x in self.region_data if x[0] not in idx]
    
    self.updateTable()
      

  def updateTable(self):
  
    
    if self.isVisible():
      self.table.setObjects(self.region_data)
      self.table.update()
    
    self.save_regions()
    self.mainApp.updateContents()
    
    
  def updateContents(self):
  
    if self.isVisible():
      self.load_regions()
      self.table.setObjects(self.region_data)
      self.table.update()
    
