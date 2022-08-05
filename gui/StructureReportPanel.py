from PySide import QtCore, QtGui

from numpy import array

from gui.qtgui.Frame import Frame
from gui.qtgui.LabelFrame import LabelFrame
from gui.qtgui.Table import ObjectTable, Column
from gui.qtgui.ButtonList import ButtonList
from gui.qtgui.MessageDialog import showOkCancel
from gui.qtgui.Graph import GraphAxis, GraphDataSet, Graph
from gui.qtgui.Label import Label
from gui.qtgui.PulldownList import PulldownList


class StructureReportPanel(Frame):
  
  def __init__(self, parent, mainApp):
  
    Frame.__init__(self, parent)
    
    self.mainApp = mainApp
    self.structure = None
    self.model = None
    self.chromosome = None

    row = 0
    frame = LabelFrame(self, 'Structures', grid=(row,0))
   
    columns = [Column('ID', self.getStructId),
               Column('Name', self.getStructName,
                      setEditValue=self.setStrucName),
               Column('Num\nmodels', self.getStructNumModels),
               Column('Num\nparticles', self.getStructNumParticles, format='{:,}'),
               Column('Num\nrestraints', self.getStructNumRestraints, format='{:,}'),
               Column('Num\nchromos', self.getStructNumChromos),
               Column('Chromsomes', self.getStructChromos),
               ]
    
    self.structuresTable = ObjectTable(frame, columns, None,
                                       callback=self.selectStructure,
                                       sortCol=None, multiSelect=True,
                                       grid=(0,0))
    
    frame = LabelFrame(self, 'Models', grid=(row,1))
    
    columns = [Column('#', self.getModelId),
               Column('Displayed?', self.getModelSelected,
                      setEditValue=self.toggleModelDisplayed),
               Column('RMSD', self.getModelRmsd),
               Column('Long & short axes', self.getModelSize, format='%.3f x %.3f'),
               ]
    
    self.modelsTable = ObjectTable(frame, columns, None,
                                   callback=None, sortCol=None,
                                   multiSelect=True, grid=(0,0))
                                   
    texts = ['Display selected', 'Remove selected']
    callbacks = [self.selectModels, self.deleteModels]
    buttons = ButtonList(frame, texts, callbacks, grid=(1,0))
                                   
    row += 1
 
    frame = LabelFrame(self, 'Restraint violations', grid=(row,0), gridSpan=(1,2))
    
    xAxis = GraphAxis('Sequence position (Mb)',   labels=None, ticks=True)
    yAxis = GraphAxis('Restraint violation', labels=None,
                      grid=False, ticks=True, valRange=(0.0, 5.0))
    yAxis2 = GraphAxis('Structure RMSD',     labels=None,
                       grid=False, ticks=True, valRange=(0.0, 4.0))
    
    label = Label(frame, ' Chromosome: ', grid=(0,0))
    self.chromoPulldown = PulldownList(frame, callback=self.selectChromosome, grid=(0,1))
    
    self.reportGraph = Graph(frame, (xAxis, yAxis, yAxis2), [],
                             size=(400, 200), title='Chromosome',
                             fgColor=QtGui.QColor(255, 255, 255, 255),
                             bgColor=QtGui.QColor(0, 0, 0, 255),
                             mgColor=QtGui.QColor(64, 64, 64, 255),
                             grid=(1,0), gridSpan=(1,3)) 
    
    self.layout().setRowStretch(row, 2)
    frame.layout().setColumnStretch(2, 2)


  def selectChromosome(self, chromo):
    
    if chromo != self.chromosome:
      self.chromosome = chromo
      self._redrawGraphs()
    

  def updateChromosomes(self):
    
    chromo = self.chromosome
    nuc = self.mainApp.nuc
    coordsGroup = nuc._getCoordsGroup(self.structure)

    chromos = nuc.getChromosomes()
    chromos = [c for c in chromos if c in coordsGroup and coordsGroup[c].shape[1]]
       
    if nuc and chromos:
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
 
  
  def toggleModelDisplayed(self, i, val):
    
    nuc = self.mainApp.nuc
    models = set(nuc.structures[self.structure].attrs['displayModels'])
    
    if val:
      models.add(i)
    else:
      models.remove(i)
      
    nuc.structures[self.structure].attrs['displayModels'] = sorted(models)
  
  
  def deleteModels(self):
  
    selected = self.modelsTable.getSelectedObjects()
 
    if selected:
      n = len(selected)
 
      if n == 1:
        msg = 'Remove %d selected model?'
      else:
        msg = 'Remove %d selected models?'
 
      if not showOkCancel('Confirm', msg % n, parent=self):
        return

      nuc = self.mainApp.nuc
      nuc.removeModels(self.structure, selected)
      
      self.modelsTable.clearSelection()
      self.updateContents()
  
  
  def selectModels(self):
  
    selected = self.modelsTable.getSelectedObjects()
 
    if selected:
      nuc = self.mainApp.nuc
      nuc.structures[self.structure].attrs['displayModels'] = selected
      
      self.updateModels()
  
  
  def getModelRmsd(self, i):
    
    return float(self.rmsds[i])
    
  
  def getModelSize(self, i):
  
    nuc = self.mainApp.nuc
    return nuc.getModelSize(i, structure=self.structure)

  
  def getModelId(self, i):
  
    return i


  def getModelSelected(self, i):
    
    if self.structure:
      nuc = self.mainApp.nuc
      return i in nuc.structures[self.structure].attrs['displayModels']
  
    
  def setStrucName(self, structure, name):
  
    nuc = self.mainApp.nuc
    nuc.structures[structure].attrs['name'] = name
 
  
  def selectStructure(self, obj, row, col):
    
    if self.structure != obj:
      nuc = self.mainApp.nuc
      self.structure = obj
      
      if obj:
        self.rmsds = nuc.calcModelRmsds(structure=self.structure)[0]
      else:
        self.rmsds = []
      
      self.updateModels()
      self.updateChromosomes()
      self._redrawGraphs()
  
  
  def getStructId(self, code):
    
    return code


  def getStructName(self, code):
    
    nuc = self.mainApp.nuc
    return nuc.structures[code].attrs['name']
 
 
  def getStructNumModels(self, structure):
    
    nuc = self.mainApp.nuc
    return nuc.getNumModels(structure)
 
  
  def getStructNumParticles(self, structure):
  
    nuc = self.mainApp.nuc
    coordsGroup = nuc._getCoordsGroup(structure)
    
    n = 0
    for chromo in coordsGroup:
      shape = coordsGroup[chromo].shape
      
      if len(shape) != 3:
        continue
      
      n += shape[1]
    
    return n
    
    
  def getStructNumRestraints(self, structure):
  
    nuc = self.mainApp.nuc
    restGroup = nuc._getParticleGroup(structure)
  
    n = 0
    for chromoA in restGroup:
      subGroup = restGroup[chromoA]
 
      for chromoB in subGroup:
        n += subGroup[chromoB].shape[0]

    return n
   
  
  def getStructNumChromos(self, structure):
  
    nuc = self.mainApp.nuc
    coordsGroup = nuc._getCoordsGroup(structure)
    
    return len(coordsGroup.keys())
    
  
  def getStructChromos(self, structure):
  
    nuc = self.mainApp.nuc
    chromos = nuc.getChromosomes(structure)
    
    return ', '.join(chromos)
  
  
  def updateModels(self):
  
    if self.structure:
      nuc = self.mainApp.nuc
      models = list(range(nuc.getNumModels(self.structure)))
    else:
      models = []
    
    self.modelsTable.setObjects(models)
    self.modelsTable.update()
 
 
  def updateContents(self):
  
    if self.isVisible():
      nuc = self.mainApp.nuc
      codes = list(nuc.structures.keys())
      codes.sort(key=lambda a:'%5d' % int(a))
      
      self.structure = None
      self.structuresTable.setObjects(codes)
      self.structuresTable.update()
      
      if codes:
        self.selectStructure(codes[0], 0, 0)
      else:
        self.updateModels()
        self.updateChromosomes()
        self._redrawGraphs()


  def _redrawGraphs(self):
    
    # TBC : trans violations
    
    if not self.structure:
      return
    
    nuc = self.mainApp.nuc
    dataSets = []
    chrA = self.chromosome
    
    if chrA and nuc.getNumModels(self.structure):
      violDict = nuc.calcRestraintViolations([chrA], cis=True, trans=False, upperOnly=True,
                                             reportAll=False, structure=self.structure)
 
      rmsds, atomRmsds = nuc.calcModelRmsds(chromosomes=[chrA], structure=self.structure)
 
      if len(atomRmsds):
        particGroup = nuc._getParticleGroup(self.structure)      
        positions = particGroup[chrA]['positions']
        binSize = max(1, min(positions[-1]/100.0, 2e6))
        binData = {}
        xyVals = []
 
        for i, pos in enumerate(positions):
          bin = int(pos//binSize)
          rmsd = atomRmsds[i]
 
          if bin in binData:
            binData[bin].append(rmsd)
          else:
            binData[bin] = [rmsd,]
 
        for bin in sorted(binData):
          vals = array(binData[bin])
 
          x = bin*binSize/1e6
          y  = vals.mean()
          #dy = vals.std()
 
          xyVals.append((x,y))
 
        if xyVals:
          ds = GraphDataSet(xyVals, 'RMSD', '#808080', plotType='histogram',
                          symbol='circle', symbolSize=0.5, secondAxis=True)
          dataSets.append(ds)
 
      if violDict:
        xyVals = []
        for posA, posB, delta, target in violDict[(chrA, chrA)]:
          frac = delta/target
          xyVals.append((posA/1e6, frac))
 
        if xyVals:
          ds = GraphDataSet(xyVals, 'Violation', '#FF0000',  plotType='scatter',
                            symbol='circle', symbolSize=1.0)
          dataSets.append(ds)
      
    self.reportGraph.xAxis.valRange = (0.0, 1.0)
    self.reportGraph.yAxis.valRange = (0.0, 5.0)
    self.reportGraph.yAxis2.valRange = (0.0, 4.0)
    self.reportGraph.updateData(dataSets, title='Chromosome ' + (chrA or ''))
