from PySide import QtCore, QtGui

from gui.qtgui.ButtonArray import ButtonArray
from gui.qtgui.Button import Button
from gui.qtgui.ButtonList import ButtonList
from gui.qtgui.CheckButton import CheckButton
from gui.qtgui.Colors import ColorDialog
from gui.qtgui.Entry import FloatEntry
from gui.qtgui.LabelFrame import LabelFrame
from gui.qtgui.Label import Label
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.SpinBox import FloatSpinBox, IntSpinBox
from gui.qtgui.Table import ObjectTable, Column

from colorsys import hsv_to_rgb

from numpy import array

COLOR_MODES = ['Seq. Position','Chromsome ID','Density','Model number','Data track', 'Faint']
DISPLAY_MODES = ['Ball and Stick', 'Line', 'Tube', 'Surface']

class DisplayOptionsPanel(QtGui.QWidget):

  def __init__(self, mainApp, parent=None):
  
    QtGui.QWidget.__init__(self, parent)
    
    self.mainApp = mainApp
    
    row = 0
    frame = LabelFrame(self, 'Colouring:', grid=(row, 0))
    
    buttons = ButtonArray(frame, COLOR_MODES, range(6), icons=None,
                          callback=self.setColorMode,
                          toggled=True, radio=True, selected=0, tipTexts=None, 
                          maxCols=2, grid=(0,0))
    self.colorButtons = buttons
    
    row += 1
    frame = LabelFrame(self, 'Display mode:', grid=(row, 0))
     
    buttons = ButtonArray(frame, DISPLAY_MODES, range(4), icons=None,
                          callback=self.setDisplayMode,
                          toggled=True, radio=True, selected=0, tipTexts=None, 
                          maxCols=5, grid=(0,0), gridSpan=(1, 7))
    self.displayButtons = buttons
    
    label = Label(frame, 'Ball radius:', grid=(1,0))
    self.sphSizeEntry = FloatSpinBox(frame, 0.0, minValue=0.1, maxValue=10.0, step=0.5,
                                     callback=self.setSphereSize, multiplier=1.1,
                                     grid=(1,1), sticky='ew')

    label = Label(frame, 'Stick width:', grid=(1,2))
    self.stickSizeEntry = FloatSpinBox(frame, 0.5, minValue=0.0, maxValue=1.0, step=0.25,
                                       callback=self.setStickSize, grid=(1,3), sticky='ew')

    label = Label(frame, 'Sphere detail:', grid=(1,4))
    self.sphDetailEntry = IntSpinBox(frame, 0.5, minValue=0, maxValue=4, step=1,
                                    callback=self.setSphereDetail, grid=(1,5), sticky='ew')
    
    
    label = Label(frame, 'Line width:', grid=(2,0))
    self.lineWidthEntry = FloatSpinBox(frame, 1.0, minValue=0.5, maxValue=5.0, step=0.5,
                                       callback=self.setLineWidth, grid=(2,1), sticky='ew')

    label = Label(frame, 'Line smooth:', grid=(2,2))
    self.lineSmoothEntry = IntSpinBox(frame, 0, minValue=0, maxValue=4, step=1,
                                      callback=self.setLineSmooth, grid=(2,3), sticky='ew')
    
    
    label = Label(frame, 'Tube width:', grid=(3,0))
    self.tubeSizeEntry = FloatSpinBox(frame, 0.0, minValue=0.1, maxValue=10.0, step=0.5,
                                     callback=self.setTubeSize, multiplier=1.1,
                                     grid=(3,1), sticky='ew')
                                     
    label = Label(frame, 'Tube smooth:', grid=(3,2))
    self.tubeSmoothEntry = IntSpinBox(frame, 0, minValue=0, maxValue=4, step=1,
                                      callback=self.setTubeSmooth, grid=(3,3), sticky='ew')

    label = Label(frame, 'Tube detail:', grid=(3,4))
    self.tubeDetailEntry = IntSpinBox(frame, 0, minValue=0, maxValue=4, step=1,
                                      callback=self.setTubeDetail, grid=(3,5), sticky='ew')
    
    
    frame.layout().setColumnStretch(6, 2)
    
    
    row += 1
    frame = LabelFrame(self, 'Extras:', grid=(row, 0), stretch=(0,1))
    
    self.showCisCheck = CheckButton(frame, 'Cis restraints', grid=(0,0),
                                    callback=self.setRestraintLinesCis)
                                          
    self.showTransCheck = CheckButton(frame, 'Trans restraints', grid=(0,1),
                                      callback=self.setRestraintLinesTrans)
                                            
    self.showLabelCheck = CheckButton(frame, 'Text labels', callback=self.setTextLabels, grid=(0,2))
    self.showScaleCheck = CheckButton(frame, 'Scalebar', callback=self.setScalebar, grid=(0,3))
                 
    row += 1
    self.layout().setRowStretch(row, 2)
    frame = LabelFrame(self, 'Chromosome specific options:', grid=(row, 0), stretch=(2,1))

    columns = [Column('Chromosome', self.getChromoName,
                      tipText='Chromosome identifier/name'),
               Column('Show?', self.getChromoShown,
                      getEditValue=None,
                      setEditValue=self.setChromoShown,
                      tipText='Toggle whether chromosome is shown'),
               Column('ID Colour', None,
                      getColor=self.getChromoColor,
                      getEditValue=self.getChromoColor,
                      setEditValue=self.setChromoColor,
                      editClass=ColorDialog,
                      tipText='The identifying single colour for the chromsosome'),
               Column('Colour mode', self.getChromoColorMode,
                      getEditValue=self.chooseChromoColorMode,
                      setEditValue=self.setChromoColorMode,
                      tipText='How the chromosome is coloured in a chromosome specific structure display'),
               Column('Display mode', self.getChromoDisplayMode,
                      getEditValue=self.chooseChromoDisplayMode,
                      setEditValue=self.setChromoDisplayMode,
                      tipText='How chromosome is rendered in a chromosome specific structure display'),
               ]
    
    self.chromosome = None
    self.chromoTable = ObjectTable(frame, columns, [], callback=self.selectChromo,
                                   multiSelect=True, grid=(0,0), stretch=(2,1))
  
    texts = ['Toggle selected', 'Show all']
    callbacks = [self.toggleSelected, self.showAll]    
  
    button = ButtonList(frame, texts=texts, callbacks=callbacks)

    texts = ['Highlight selected', 'Restore defaults']
    callbacks = [self.highlightSelected, self.restoreDefaults]    
  
    button = ButtonList(frame, texts=texts, callbacks=callbacks)

    row = 0
    frame = LabelFrame(self, 'Models:', grid=(row, 1), gridSpan=(4, 1))
    columns = [Column('Model', self.getModelName,
               tipText='Help text'),
               Column('Show?', self.getModelShown,
               setEditValue=self.setModelShown, 
               tipText='Help text'),
               Column('Colour', None,
               getColor=self.getModelColor,
               tipText='Help text'),
               ]
    
    self.model = None
    self.modelTable = ObjectTable(frame, columns, [], multiSelect=True,
                                  callback=self.selectModel, grid=(0,0)) 
    texts = ['Hide selected', 'Show selected', 'Show all']
    callbacks = [self.hideModels, self.showModels, self.showAllModels]
    butonList = ButtonList(frame, texts=texts, callbacks=callbacks, grid=(1,0))

  
  # Detail levels
  
  def setSphereDetail(self, value):
    
    if self.mainApp.nuc.setDisplayDetailLevels(ballDetail=value):
      self.mainApp.updateContents()


  def setLineSmooth(self, value):
  
    if self.mainApp.nuc.setDisplayDetailLevels(lineSmooth=value):
      self.mainApp.updateContents()


  def setTubeSmooth(self, value):
  
    if self.mainApp.nuc.setDisplayDetailLevels(tubeSmooth=value):
      self.mainApp.updateContents()
  
  
  def setTubeDetail(self, value):
  
    if self.mainApp.nuc.setDisplayDetailLevels(tubeDetail=value):
      self.mainApp.updateContents()
  
  
  # Sizes
            
  def setSphereSize(self, radius):
  
    if self.mainApp.nuc.setDisplaySizes(ball=radius):
      self.mainApp.updateContents()


  def setStickSize(self, frac):
       
    if self.mainApp.nuc.setDisplaySizes(stick=frac):
      self.mainApp.updateContents()


  def setLineWidth(self, value):

    if self.mainApp.nuc.setDisplaySizes(line=value):
      self.mainApp.updateContents()
    

  def setTubeSize(self, radius):
  
    if self.mainApp.nuc.setDisplaySizes(tube=radius):
      self.mainApp.updateContents()
 
  
  # Other options
  
  def setColorMode(self, value):
    
    refresh = False
    
    nuc = self.mainApp.nuc
    for chromo in nuc.getChromosomes():
      if nuc.setChromoDisplayParams(chromo, colorMode=value):
        refresh = True
    
    if refresh:
      self.mainApp.updateContents()
    
    
  def setDisplayMode(self, value):
    
    refresh = False
    
    nuc = self.mainApp.nuc
    for chromo in nuc.getChromosomes():
      if nuc.setChromoDisplayParams(chromo, displayMode=value):
        refresh = True
    
    if refresh:
      self.mainApp.updateContents()
  
  
  def setRestraintLinesCis(self, isSelected):
        
    if self.mainApp.nuc.setRestraintsDisplayed(cis=isSelected):
      self.mainApp.updateContents()
  
  
  def setRestraintLinesTrans(self, isSelected):
    
    if self.mainApp.nuc.setRestraintsDisplayed(trans=isSelected):
      self.mainApp.updateContents()
  
  
  def setTextLabels(self, isSelected):
    
    if self.mainApp.nuc.setTextDisplayed(isSelected):
      self.mainApp.updateContents()
  
    
  def setScalebar(self, isSelected):
  
    if self.mainApp.nuc.setScaleDisplayed(isSelected):
      self.mainApp.updateContents()
  
  
  def showEvent(self, event):
  
    self.updateContents()
    
    
  def _updateRadioButtons(self):
  
    nuc = self.mainApp.nuc
 
    if nuc:
      chromoColors = set()
      chromoDisplay = set()
      
      for chromo in nuc.getChromosomes():
         isShown, useLabels, colMode, dispMode = nuc.getChromoDisplayParams(chromo)
         chromoColors.add(colMode)
         chromoDisplay.add(dispMode)
      
      if len(chromoColors) == 1:
        self.colorButtons.setSelected(chromoColors)
        
      else:
        self.colorButtons.setSelected([])
      
      if len(chromoDisplay) == 1:
        self.displayButtons.setSelected(chromoDisplay)
        
      else:
        self.displayButtons.setSelected([]) 


  def updateContents(self):
    
    if self.isVisible():
      nuc = self.mainApp.nuc
 
      if nuc:
        attrs = nuc.display.attrs
        models = list(range(nuc.getNumModels()))

        sphDetail  = attrs['detailLevel'][0]
        lineSmooth = attrs['detailLevel'][1]
        tubeSmooth = attrs['detailLevel'][2]
        tubeDetail = attrs['detailLevel'][3]
 
        sphereSize = attrs['sizes'][0]
        stickSize  = attrs['sizes'][1]
        lineWidth  = attrs['sizes'][2]
        tubeSize   = attrs['sizes'][3]
        
        showCis    = attrs['options'][2]
        showTrans  = attrs['options'][3]
        showLabels = attrs['options'][4]
        showScale  = attrs['options'][5]
         
        self._updateRadioButtons()
        
      else:
        models = []
 
        lineSmooth = 2
        tubeSmooth = 2
        sphDetail  = 2
 
        sphereSize = 1.0
        stickSize  = 0.5
        lineWidth  = 1.0
        tubeSize   = 1.0
 
      self.modelTable.setObjects(models)
 
      self.sphSizeEntry.set(sphereSize, doCallback=False)
      self.stickSizeEntry.set(stickSize, doCallback=False)
      self.sphDetailEntry.set(sphDetail, doCallback=False)
 
      self.lineWidthEntry.set(lineWidth, doCallback=False)
      self.lineSmoothEntry.set(lineSmooth, doCallback=False)
 
      self.tubeSizeEntry.set(tubeSize, doCallback=False)
      self.tubeSmoothEntry.set(tubeSmooth, doCallback=False)
      self.tubeDetailEntry.set(tubeDetail, doCallback=False)
      
      self.showCisCheck.set(showCis, doCallback=False)
      self.showTransCheck.set(showTrans, doCallback=False)
      self.showLabelCheck.set(showLabels, doCallback=False)
      self.showScaleCheck.set(showScale, doCallback=False)
      
      if nuc:
        chromosomes = nuc.getChromosomes()
      else:
        chromosomes = []
 
      self.chromoTable.setObjects(chromosomes)
  
  def hideModels(self):
  
    hideModels = set(self.modelTable.selectedObjects)
    
    if hideModels:
      showModels1 = set(self.mainApp.nuc.display.attrs['models'])
      showModels2 = showModels1 - hideModels
      
      if self.mainApp.nuc.setModelsDisplayed(showModels2):
        self.mainApp.updateContents()
  
  
  def showAllModels(self):
    
    nuc = self.mainApp.nuc
    showModels = set(range(nuc.getNumModels()))
    
    if nuc.setModelsDisplayed(showModels):
      self.mainApp.updateContents()
    

  def showModels(self):
  
    showModels = set(self.modelTable.selectedObjects)
    
    if showModels:
      showModels1 = set(self.mainApp.nuc.display.attrs['models'])
      showModels2 = showModels1 | showModels
     
      if self.mainApp.nuc.setModelsDisplayed(showModels2):
        self.mainApp.updateContents()
 
  
  def selectModel(self, obj, row, col):
  
    self.model = obj
    
    
  def getModelName(self, model):
  
    return str(model)  


  def getModelShown(self, model):
  
    showModels =  set(self.mainApp.nuc.display.attrs['models'])
    
    return model in showModels
      
    
  def setModelShown(self, model, value):
    
    showModels =  set(self.mainApp.nuc.display.attrs['models'])
    
    if value and (model not in showModels):
      showModels.add(model)
      
      self.mainApp.nuc.setModelsDisplayed(showModels)
      self.mainApp.updateContents()
     
    elif not value and (model in showModels):
      showModels.remove(model)

      self.mainApp.nuc.setModelsDisplayed(showModels)
      self.mainApp.updateContents()
  
  
  def getModelColor(self, model):
    
    nuc = self.mainApp.nuc
    prop = model/float(nuc.getNumModels() or 1.0)
    r, g, b = [int(x*255.0) for x in hsv_to_rgb(prop*0.85, 1.0, 1.0)]
    
    return QtGui.QColor(r, g, b, 255)
      
      
  def selectChromo(self, obj, row, col):
    
    self.chromosome = obj
   
   
  def getChromoName(self, chromo):
  
    return chromo
  
  
  def setChromoColor(self, chromo, colorObj):
    
    if colorObj: # not cancelled
      rgba =  colorObj.getRgbF()
      
      if self.mainApp.nuc.setChromoColor(chromo, rgba):
        self.mainApp.updateContents()
      
      
  def getChromoColor(self, chromo):
    
    nuc = self.mainApp.nuc
    
    if nuc:
      color = nuc.getChromoColor(chromo)
      color = array(color*255, int) # uint8 , or set bounds
      return QtGui.QColor(*color)
      
    else:
      return QtGui.QColor(255, 255, 255, 255)
  
  
  def setChromoShown(self, chromo, boolean):
  
    self.mainApp.nuc.setChromoDisplayed(chromo, boolean)
    self.mainApp.updateContents()
    
  
  def getChromoShown(self, chromo):

    return self.mainApp.nuc.getChromoDisplayed(chromo)
     
     
  def getChromoLabels(self, chromo):
    
    nuc = self.mainApp.nuc
    isShown, useLabels, colMode, dispMode = nuc.getChromoDisplayParams(chromo)
    
    return useLabels
    
    
  def getChromoColorMode(self, chromo):
    
    nuc = self.mainApp.nuc
    isShown, useLabels, colMode, dispMode = nuc.getChromoDisplayParams(chromo)
    
    return COLOR_MODES[colMode]
  
  
  def chooseChromoColorMode(self, chromo):
  
    nuc = self.mainApp.nuc
    
    objects = range(len(COLOR_MODES))
    index = nuc.getChromoDisplayParams(chromo)[2]
    
    return COLOR_MODES, objects, index
  
  
  def setChromoColorMode(self, chromo, value):
  
    self.mainApp.nuc.setChromoDisplayParams(chromo, colorMode=value)
    self._updateRadioButtons()
    
    
  def getChromoDisplayMode(self, chromo):
    
    nuc = self.mainApp.nuc
    isShown, useLabels, colMode, dispMode = nuc.getChromoDisplayParams(chromo)
    
    return DISPLAY_MODES[dispMode]


  def chooseChromoDisplayMode(self, chromo):
  
    nuc = self.mainApp.nuc
    
    objects = range(len(DISPLAY_MODES))
    index = nuc.getChromoDisplayParams(chromo)[3]
    
    return DISPLAY_MODES, objects, index
  
  
  def setChromoDisplayMode(self, chromo, value):
  
    self.mainApp.nuc.setChromoDisplayParams(chromo, displayMode=value)
    self._updateRadioButtons()
 
 
  def highlightSelected(self):
  
    nuc = self.mainApp.nuc
    
    if nuc:
      chromosomes = nuc.getChromosomes()
      selection = set(self.chromoTable.getSelectedObjects())
      
      for chromo in chromosomes:
        if chromo in selection:
          nuc.setChromoDisplayParams(chromo, 1, 0, 1, 2)
        else:
          nuc.setChromoDisplayParams(chromo, 1, 0, 5, 1)
      
      self.mainApp.updateContents()
     
     
  def restoreDefaults(self):
  
    nuc = self.mainApp.nuc
    
    if nuc:
      chromosomes = nuc.getChromosomes()
      
      for chromo in chromosomes:
        nuc.setChromoDisplayParams(chromo, 1, 0, 1, 1)
      
      nuc.resetChromoColors()
      
      self.mainApp.updateContents()
  
  
  def toggleSelected(self):
  
    nuc = self.mainApp.nuc
    
    if nuc:
      selection = set(self.chromoTable.getSelectedObjects())
      
      for chromo in selection:
        nuc.setChromoDisplayed(chromo, not nuc.getChromoDisplayed())
      
      self.mainApp.updateContents()
  
  
  def showAll(self):
  
    nuc = self.mainApp.nuc
    
    if nuc:
      chromosomes = nuc.getChromosomes()
      
      for chromo in chromosomes:
        nuc.setChromoDisplayed(chromo, True)
        
      self.mainApp.updateContents()
        
  

