from PySide import QtCore, QtGui
from time import time

from math import log, exp
from numpy import array

import multiprocessing

from gui.qtgui.ButtonArray import ButtonArray
from gui.qtgui.Button import Button
from gui.qtgui.ButtonList import ButtonList
from gui.qtgui.CheckButton import CheckButton
from gui.qtgui.Entry import FloatEntry
from gui.qtgui.Frame import Frame
from gui.qtgui.LabelFrame import LabelFrame
from gui.qtgui.Label import Label
from gui.qtgui.MessageDialog import showInfo
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.RadioButtons import RadioButtons
from gui.qtgui.SpinBox import FloatSpinBox, IntSpinBox
from gui.qtgui.Table import ObjectTable, Column


class StructureCalcPanel(Frame):

  def __init__(self, parent, mainApp, updateFunc=None):
  
    Frame.__init__(self, parent)
    
    self.mainApp = mainApp
    self.updateFunc = updateFunc
    self.structure = None
    self.groupName = None
    self.timer = QtCore.QTimer()
    self.timer.setInterval(1000)
    self.timer.timeout.connect(self.checkStatus)    
    self.calcJob = None
    self.calcChromos = None
    
    row = 0
    frame0 = LabelFrame(self, 'Structure calculation settings', grid=(row,0))
    
    Label(frame0, 'Target structure: ', grid=(0,0), hAlign=self.right, stretch=(0,1))
    self.structPulldown = PulldownList(frame0, grid=(0, 1), stretch=(0,0),
                                       callback=self.changeStructure)
    
    
    Label(frame0, 'Contact group: ', grid=(0,2), hAlign=self.right, stretch=(0,1))
    self.contGroupPulldown = PulldownList(frame0, grid=(0, 3), stretch=(0,0),
                                          callback=self.changeContactGroup)
    
    
    
    row += 1
    frame1 = LabelFrame(self, 'Calculation presets', grid=(row,0))
    
    texts = ['High-res anneal', 'Quick low-res anneal', 'Minimisation']
    callbacks = [self.fullAnneal, self.quickAnneal, self.minimise]
    buttons1 = ButtonList(frame1, texts, callbacks, grid=(0,0), hAlign=self.left)
     
    
    row += 1
    frame2 = LabelFrame(self, 'Model parameters', grid=(row,0))
    
    sRow = 0
    Label(frame2, 'Number of models: ', grid=(sRow,0))
    self.numModelsEntry = IntSpinBox(frame2, value=1, minValue=1, maxValue=1000, step=1,
                                     multiplier=2, grid=(sRow, 1))
    
    Label(frame2, 'Power law: ', hAlign=self.right, grid=(sRow,2))
    self.powerLawEntry = FloatSpinBox(frame2, value=-1.0, minValue=-4.0, maxValue=0.0, step=1.0,
                                      callback=None, grid=(sRow,3))
                                    
    sRow += 1
    Label(frame2, 'Regular backbone: ', hAlign=self.right, grid=(sRow,0))
    self.bboneRegCheck = CheckButton(frame2, '', callback=self.checkSpacing,
                                       selected=False, grid=(sRow,1))
  
    Label(frame2, 'Restraint distance: ', hAlign=self.right, grid=(sRow,2))
    self.restrDistEntry = FloatSpinBox(frame2, value=0.5, minValue=1e-8, maxValue=128.0, step=1,
                                       multiplier=2.0, grid=(sRow,3))
    
    sRow += 1
    label2 = Label(frame2, 'Binned restraints? ', hAlign=self.right, grid=(sRow,0))
    self.restrBinnedCheck = CheckButton(frame2, '', selected=False, grid=(sRow, 1),
                                        callback=self.checkBinned)

    Label(frame2, 'Distance error: ', hAlign=self.right, grid=(sRow,2))
    self.restrErrEntry = FloatSpinBox(frame2, value=20.0, minValue=0.0, maxValue=100.0, step=5.0,
                                      suffix='%', grid=(sRow,3))
    
    
    sRow += 1
    Label(frame2, 'Backbone spacing: ', hAlign=self.right, grid=(sRow,0))
    self.bboneSpaceEntry = IntSpinBox(frame2, value=100, minValue=1, maxValue=int(2e5), step=1,
                                      suffix=' kb', multiplier=2, grid=(sRow,1))

    Label(frame2, 'Unbinned scale: ', hAlign=self.right, grid=(sRow,2))
    self.seqScaleEntry = IntSpinBox(frame2, value=10, minValue=1, maxValue=10000, step=1,
                                    callback=None, suffix=' kb', prefix=None, multiplier=2,
                                    grid=(sRow,3))
                                   
    #frame2.layout().setColumnStretch(4, 2)
       

    row += 1
    frame3 = LabelFrame(self, 'Annealing parameters', grid=(row,0))
    
    sRow += 0
    Label(frame3, 'Max. temperature: ', hAlign=self.right, grid=(sRow,0))
    self.tempMaxEntry =  IntSpinBox(frame3, value=5000, minValue=1, maxValue=10000,
                                    callback=self.checkTemps, step=1, multiplier=2, grid=(sRow,1))
                                    
    Label(frame3, 'Hierarchical\nscaling? ', hAlign=self.right, grid=(sRow,2))
    self.hierProtocolCheck = CheckButton(frame3, '', selected=True, grid=(sRow, 3))
    
    sRow += 1
    Label(frame3, 'Min. temperature: ', hAlign=self.right, grid=(sRow,0))
    self.tempMinEntry = IntSpinBox(frame3, value=10, minValue=0, maxValue=10000,
                                   callback=self.checkTemps, step=1, multiplier=2, grid=(sRow,1))
    
    Label(frame3, ' Starting scale: ', hAlign=self.right, grid=(sRow,2))
    self.hierStartEntry = IntSpinBox(frame3, value=4, minValue=1, maxValue=100,
                                     step=1, suffix=' Mb', multiplier=2, grid=(sRow,3))

    sRow += 1
    Label(frame3, 'Cooling steps: ', hAlign=self.right, grid=(sRow,0))
    self.tempStepsEntry = IntSpinBox(frame3, value=1024, minValue=1, maxValue=10000,
                                     step=1, multiplier=2, grid=(sRow,1))
    
    Label(frame3, 'Scaling steps: ', hAlign=self.right, grid=(sRow,2))
    self.hierStepsEntry = IntSpinBox(frame3, value=4, minValue=1, maxValue=10,
                                     step=1, grid=(sRow,3))

    sRow += 1
    Label(frame3, 'Dynamics steps: ', hAlign=self.right, grid=(sRow,0))
    self.dynStepsEntry = IntSpinBox(frame3, value=64, minValue=1, maxValue=10000,
                                    step=1, multiplier=2, grid=(sRow,1))
    
    #frame3.layout().setColumnStretch(4, 2)
    
    row += 1
    frame4 = LabelFrame(self, 'Starting structure', grid=(row,0))

    texts = ['Random sphere volume', 'Random walk', 'Previous coordinates']
    self.startStrucRadio = RadioButtons(frame4, texts=texts, callback=None, direction='v',
                                        selectedInd=0, grid=(0,0), gridSpan=(3,1))
                                        
    Label(frame4, 'Sphere radius: ', hAlign=self.right, grid=(0,1))
    self.randRadEntry = FloatSpinBox(frame4, value=100.0, minValue=0.1, maxValue=1e4,
                                     step=1, grid=(0,2))

    Label(frame4, 'Random number seed: ', hAlign=self.right, grid=(1,1))
    self.randSeedEntry = IntSpinBox(frame4, value=0, minValue=0,
                                    step=1, grid=(1,2))
                                    
    
    row += 1
    frame5 = LabelFrame(self, 'Launch calculation', grid=(row,0))

    Label(frame5, 'Number of CPU cores: ', grid=(0,0), hAlign=self.right, stretch=(0,1))
    cpuOpts = list(range(1,1+multiprocessing.cpu_count()))
    self.numCpusPulldown = PulldownList(frame5, [str(x) for x in cpuOpts], cpuOpts,
                                        grid=(0, 1), stretch=(0,0))
    
    Label(frame5, 'Visualised calculation: ', grid=(1,0), hAlign=self.right)
    self.startButton = Button(frame5, 'Go!',
                              callback=self.startCalculation,
                              icon=None, # QtGui.QIcon('gui/icons/dialog-ok.png'),
                              grid=(1,1), bgColor='#F0D0B0')

    Label(frame5, 'Background calculation: ', grid=(2,0), hAlign=self.right)
    self.startButtonPara = Button(frame5, 'Go!',
                              callback=self.startCalculationPara,
                              icon=None, # QtGui.QIcon('gui/icons/dialog-ok.png'),
                              grid=(2,1), bgColor='#B0E0B0')

    Label(frame5, 'Individiual chromosome calculation: ', grid=(3,0), hAlign=self.right)
    self.startButtonChromo = Button(frame5, 'Go!',
                              callback=self.startCalculationChromo,
                              icon=None, # QtGui.QIcon('gui/icons/dialog-ok.png'),
                              grid=(3,1), bgColor='#B0B0F0')

    row += 1
    self.layout().setRowStretch(row, 2)
    self.updateContents()
    
    
  def updateContents(self):

    nuc = self.mainApp.nuc

    if nuc:
      self.loadParams()
      
      names = nuc.getContactGroupNames()
      
      if names:
        if self.groupName in names:
          index = names.index(self.groupName)
        
        elif 'singleCell' in names:
          self.groupName = 'singleCell'
          index = names.index(self.groupName)
          
        else:
          self.groupName = names[0]
          index = 0  
      
      else:
        self.groupName = None
        index = 0
     
      self.contGroupPulldown.setData(names, names, index)
 
      self.updateStructures()


  def updateStructures(self):
 
    nuc = self.mainApp.nuc
    codes, names = nuc.getStructureNames()
    
    if codes:
      if self.structure in codes:
        structure = self.structure
        index = codes.index(structure)

      elif nuc.structure.name.split('/')[-1] in codes:
        structure = nuc.structure.name.split('/')[-1]
        index = codes.index(structure)
        
      else:
        structure = codes[0]
        index = 0  
      
    else:
      structure = None
      index = 0
    
    names.append('<New>')
    codes.append(None)
    
    self.structPulldown.setData(names, codes, index)
    self.changeStructure(structure)

  
  def changeContactGroup(self, groupName):
  
    self.groupName = groupName
    
  
  def changeStructure(self, structureCode):
    
    if structureCode != self.structure:
    
      if structureCode is None:
        structGroup = self.mainApp.nuc.getStructureGroup()
        structureCode = structGroup.name.split('/')[-1]
      
      s = self.mainApp.nuc.structures
      
      #self.mainApp.nuc.setCurrentStructure(structureCode)  
      self.structure = structureCode
            
      self.loadParams()
      self.mainApp.structureOuterPanel.updateStrucSelectToolbar()
      self.updateStructures()
 
  
  def showEvent(self, event):
  
    self.updateContents()
  
    QtGui.QWidget.showEvent(self, event)
 

  def hideEvent(self, event):
  
    if self.mainApp.nuc:
      self.saveParams()
  
    QtGui.QWidget.hideEvent(self, event)
 
 
  def closeEvent(self, event):
  
    if self.mainApp.nuc:
      self.saveParams()
  
    QtGui.QWidget.closeEvent(self, event)
  
  
  def leaveEvent(self, event):
  
    if self.mainApp.nuc:
      self.saveParams()
  
  
  def loadParams(self):
    
    if not self.structure:
      return
    
    nuc = self.mainApp.nuc
    attrs = nuc._getCalculationGroup(self.structure).attrs
  
    self.numModelsEntry.set(attrs['numModels'])
    self.bboneRegCheck.set(bool(attrs['bboneReg']))
    self.bboneSpaceEntry.set(attrs['bboneSpace'])
    self.restrBinnedCheck.set(bool(attrs['restrBinned']))
    self.powerLawEntry.set(attrs['powerLaw'])
    self.seqScaleEntry.set(attrs['seqUnitScale']) 
    self.restrDistEntry.set(attrs['restrDist'])
    self.restrErrEntry.set(attrs['restrErr'] * 100.0)
    self.tempMaxEntry.set(attrs['tempMax'])
    self.tempMinEntry.set(attrs['tempMin'])
    self.tempStepsEntry.set(attrs['tempSteps'])
    self.dynStepsEntry.set(attrs['dynSteps'])
    self.hierProtocolCheck.set(bool(attrs['hierProtocol']))
    self.hierStartEntry.set(attrs['hierStart'])
    self.hierStepsEntry.set(attrs['hierSteps'])
    self.startStrucRadio.setIndex(attrs['startStruc'])
    self.randRadEntry.set(attrs['randRad'])
    self.randSeedEntry.set(attrs['randSeed'])
   
  
  def saveParams(self):
    
    if not self.structure:
      return
    
    nuc = self.mainApp.nuc
    attrs = nuc._getCalculationGroup(self.structure).attrs
  
    attrs['numModels'] = self.numModelsEntry.get()
    attrs['bboneReg'] = 1 if self.bboneRegCheck.get() else 0
    attrs['bboneSpace'] = self.bboneSpaceEntry.get()
    attrs['restrBinned'] = 1 if self.restrBinnedCheck.get() else 0
    attrs['powerLaw'] = self.powerLawEntry.get()
    attrs['seqUnitScale'] = self.seqScaleEntry.get()
    attrs['restrDist'] = self.restrDistEntry.get()
    attrs['restrErr'] = self.restrErrEntry.get() / 100.0
    attrs['tempMax'] = self.tempMaxEntry.get()
    attrs['tempMin'] = self.tempMinEntry.get()
    attrs['tempSteps'] = self.tempStepsEntry.get()
    attrs['dynSteps'] = self.dynStepsEntry.get()
    attrs['hierProtocol'] = 1 if self.hierProtocolCheck.get() else 0
    attrs['hierStart'] = self.hierStartEntry.get()
    attrs['hierSteps'] = self.hierStepsEntry.get()
    attrs['startStruc'] = self.startStrucRadio.getIndex()
    attrs['randRad'] = self.randRadEntry.get()
    attrs['randSeed'] = self.randSeedEntry.get()

    nuc.root.flush()
  
  
  def checkBinned(self, isBinned):
    
    if isBinned and not self.bboneRegCheck.get():
      self.bboneRegCheck.set(True)


  def checkSpacing(self, hasBackbone):
    
    if self.restrBinnedCheck.get() and not hasBackbone:
      self.restrBinnedCheck.set(False)
    
    
  def recalcRestraints(self):
  
    nuc = self.mainApp.nuc

    binned    = self.restrBinnedCheck.get()
    powLaw    = self.powerLawEntry.get()
    
    chromosomes = nuc.getDisplayedChromosomes()
    
    rDist = self.restrDistEntry.get()
    rErr = self.restrErrEntry.get() / 100.0
    
    lower = (1.0-rErr) * rDist
    upper = (1.0+rErr) * rDist
      
    if binned:
      bboneSep = self.bboneSpaceEntry.get() * 1000
      
    else:
      bboneSep = 0
    
    if self.groupName:
      groupName = self.groupName
    
    else:
      names = [x for x in self.mainApp.nuc.origContacts]
      if not names:
        return
      
      groupName = names[0]
        
    nuc.setRestraints(chromosomes, groupName, bboneSep=bboneSep,
                      binned=binned, scale=1.0, exponent=powLaw,
                      lower=lower, upper=upper, minCount=1, maxPopDist=5.0,
                      structure=self.structure)
                        
    
  def startCalculation(self):
  
    nuc = self.mainApp.nuc
    self.saveParams()
    
    self.startButton.setText('Running...')
    self.startButton.update()
    self.mainApp.tabbedPanel.setCurrentIndex(0)
    self.mainApp.update()
    
    structGroup = nuc.getStructureGroup(self.structure, None, True)
    structure = structGroup.name.split('/')[-1]
    numCpus = self.numCpusPulldown.get()
    job = nuc.calcStructure(self.groupName, numCpus,
                            self.updateFunc, structure=structure)
    
    self.startButton.setText('Go!')
    msg = 'Structure calculation complete.\nTime taken: %.2f seconds.'
    showInfo('Info', msg % job.timeElapsed(), parent=self)
    
    self.changeStructure(structure)


  def startCalculationChromo(self):
  
    nuc = self.mainApp.nuc
    self.saveParams()
    
    self.startButtonChromo.setText('Running...')
    self.startButtonChromo.update()
    self.mainApp.tabbedPanel.setCurrentIndex(0)
    self.mainApp.update()
   
    structGroup = nuc.getStructureGroup(self.structure, None, True)
    structure = structGroup.name.split('/')[-1]
    numCpus = self.numCpusPulldown.get()
    job = nuc.calcStructure(self.groupName, numCpus,
                            self.updateFunc, trans=False,
                            structure=structure)
    
    self.startButtonChromo.setText('Go!')
    msg = 'Structure calculation complete\nTime taken: %.2f seconds.'
    showInfo('Info', msg % job.timeElapsed(), parent=self)

    self.changeStructure(structure)


  def startCalculationPara(self):
  
    nuc = self.mainApp.nuc
    self.saveParams()
        
    self.startButton.disable()
    self.startButtonPara.disable()
    self.startButtonPara.setText('Running...')
    self.startButtonChromo.disable()
    self.startButtonPara.update()
    
    structGroup = nuc.getStructureGroup(self.structure, None)
    structure = structGroup.name.split('/')[-1]
    numCpus = self.numCpusPulldown.get()
    self.calcChromos = nuc.getDisplayedChromosomes()
    
    self.calcJob = nuc.calcStructure(self.groupName, numCpus,
                                     self.updateFunc, bgCalc=True,
                                     structure=structure)
    self.timer.start()
    
  
  def checkStatus(self):
  
    if self.calcJob and self.calcJob.isComplete():
      self.timer.stop()
      
      msg = 'Structure calculation complete!\nTime taken: %.2f seconds.'
      showInfo('Info', msg % self.calcJob.timeElapsed(), parent=self)
      
      structure = job.structure
      nuc = self.mainApp.nuc
      coords = self.calcJob.getResult()
      nuc.setAllCoords(coords, self.calcChromos, structure)
      
      self.calcJob.engine.stop()
      self.calcJob = None
 
      # superimpose ensemble
      if len(coords) > 1:
        nuc.modelAlign(chromosomes=self.calcChromos, structure=structure)
 
      nuc.save()
        
      self.changeStructure(structure)
      
      self.startButton.enable()
      self.startButtonPara.enable()
      self.startButtonPara.setText('Go!')
      self.startButtonChromo.enable()
       
      self.mainApp.structureGlWidget.construct(nuc)
      self.mainApp.structureGlWidget.update()
      
  
  def checkTemps(self, val):
  
    tMax = self.tempMaxEntry.get()
    tMin = self.tempMinEntry.get()
    
    if tMin > tMax:
      self.tempMaxEntry.set(tMin)
      self.tempMinEntry.set(tMax)
       
  
  def fullAnneal(self):
    
    self.restrBinnedCheck.set(False)
    self.bboneRegCheck.set(False)
    self.tempMaxEntry.set(5000.0)
    self.tempMinEntry.set(10.0)
    self.tempStepsEntry.set(1024)
    self.dynStepsEntry.set(64)
    self.hierProtocolCheck.set(True)
    self.hierStartEntry.set(2)
    self.hierStepsEntry.set(4)
    self.startStrucRadio.setIndex(0)
    
    
  def quickAnneal(self):
    
    self.numModelsEntry.set(1)
    self.bboneRegCheck.set(True)
    self.bboneSpaceEntry.set(500)
    self.restrBinnedCheck.set(True)
    self.tempMaxEntry.set(3000.0)
    self.tempMinEntry.set(10.0)
    self.tempStepsEntry.set(500)
    self.dynStepsEntry.set(2)
    self.hierProtocolCheck.set(True)
    self.hierStartEntry.set(4)
    self.hierStepsEntry.set(2)
    self.startStrucRadio.setIndex(0)
  
  
  def minimise(self):
    
    self.tempMaxEntry.set(100.0)
    self.tempMinEntry.set(10.0)
    self.tempStepsEntry.set(100)
    self.dynStepsEntry.set(100)
    self.hierProtocolCheck.set(False)
    self.startStrucRadio.setIndex(2)
    
