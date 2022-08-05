from PySide import QtCore, QtGui

from gui.qtgui.ButtonList import ButtonList
from gui.qtgui.Frame import Frame
from gui.qtgui.LabelFrame import LabelFrame
from gui.qtgui.Label import Label
from gui.qtgui.MessageDialog import showOkCancel
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.Table import ObjectTable, Column

from numpy import array

Qt = QtCore.Qt

STR_CATS =  ['Structure models', 'Coords per model',
             'Structure restraints', 'Cis structure restraints',
             'Trans structure restraints',]

class OverviewPanel(Frame):
  
  def __init__(self, parent, mainApp):
  
    Frame.__init__(self, parent)
    
    self.mainApp = mainApp
    self.groupName = None
    self.structure = None
    self.restraintStats = {}
    self.contactStats = {}
    self.chromoStats = {}
    
    frame1 = LabelFrame(self, 'Contacts', grid=(0,0))

    columns = [Column('Group', self.getGroupName),
               Column('Bin\nsize', self.getBinSize, format='{:,}'),
               Column('Contacts',   self.getAllCount, format='{:,}'),
               Column('Cis',   self.getCisCount, format='{:,}'),
               Column('Trans', self.getTransCount, format='{:,}'),
               Column('%\nTrans', self.getTransPercent, format='{:.2f}'),
               Column('% Cis\nisolated', self.getCisIsolated, format='%.2f'),
               Column('% Trans\nisolated', self.getTransIsolated, format='%.2f'),
               ]
 
    self.numsTable = ObjectTable(frame1, columns, [],
                                 callback=self.selectGroup, sortCol=None,
                                 multiSelect=False, grid=(0,0), gridSpan=(1,2))
    
    texts = ['Revert to original', 'Remove contacts']
    callbacks = [self.revertContacts, self.removeContacts]
    self.contactButtons = ButtonList(frame1, texts, callbacks,
                                     grid=(1, 0), gridSpan=(1,2))
    
    frame2 = LabelFrame(self, 'Chromosomes', grid=(1,0))
    self.layout().setRowStretch(1, 2)
    
    label = Label(frame2, 'Structure:', grid=(0,0))
    
    self.structPulldown = PulldownList(frame2, callback=self.changeStructure,
                                       grid=(0,1), stretch=(0,1))
   
    columns = [Column('Name', self.getChromoName),
               Column('Points', self.getChromoPoints, format='{:,}'),
               Column('Region (Mb)', self.getChromoRegion),
               Column('Contacts', self.getChromoContacts, format='{:,}'),
               Column('obs/exp', self.getExptRestr, format='%.2f'),
               Column('Cis\nrestr.', self.getChromoCis, format='{:,}'),
               Column('Trans\nrestr.', self.getChromoTrans, format='{:,}'),
               ]
               
    self.chromoTable = ObjectTable(frame2, columns, [],
                                   callback=None, sortCol=0,
                                   multiSelect=True, grid=(1,0),
                                   gridSpan=(1,2))
                                   
    texts = ['Remove selected', ]
    callbacks = [self.removeChromosome, ]
    buttons = ButtonList(frame2, texts, callbacks, grid=(2,0), gridSpan=(1,2))

    """
    frame3 = LabelFrame(self, 'Structures', grid=(0,0))

    columns = [Column('Category', self.getOverviewCategory, alignment='l'),
               Column('Count', self.getOverviewValue, alignment='r'),
               Column('% Total', self.getOverviewPercent, alignment='r')]
               
    self.numsTable = ObjectTable(frame1, columns, range(len(CATS)),
                                 callback=None, sortCol=None,
                                 multiSelect=False, grid=(0,0))
    """
  
  def selectGroup(self, obj, row, col):
  
    self.groupName = obj
    
  
  def revertContacts(self):
  
    nuc = self.mainApp.nuc
    
    if self.groupName:
      group = nuc.getContactGroup(self.groupName)
      
      if group.parent == nuc.workContacts:
        msg = 'Really revert "%s" contacts to original values?' % self.groupName
        
        if showOkCancel('Confirm', msg, parent=self):
          nuc.revertContacts(self.groupName)
          self.updateContents()
        
  
  def removeContacts(self):
  
    nuc = self.mainApp.nuc
    
    if self.groupName:
      msg = 'Really remove ALL "%s" contacts?' % self.groupName
      
      if showOkCancel('Confirm', msg, parent=self):
        nuc.removeContacts(self.groupName)
        self.groupName = None
        self.updateContents()
    
    
  def removeChromosome(self):
  
    nuc = self.mainApp.nuc
    chromos = self.chromoTable.getSelectedObjects()
    
    if chromos:
      msg = 'Really remove %d chromosomes (%s) and all thier data?' % (len(chromos), ', '.join(chromos))
      
      if showOkCancel('Confirm', msg, parent=self):
        nuc.removeChromosomes(chromos)
        self.updateContents()
    
    
  def getChromoName(self, chromo):
  
    return chromo
  
  
  def getChromoPoints(self, chromo):
    
    nuc = self.mainApp.nuc
    particGroup = nuc._getParticleGroup(self.structure)
    
    if chromo in particGroup:
      return particGroup[chromo]['positions'].shape[0]
  
  
  def getChromoRegion(self, chromo):
    
    nuc = self.mainApp.nuc
    
    if nuc and chromo in nuc.chromosomes:
      start, end = nuc.getChromosomeLimits(chromo)
      return '%.1f-%.1f' % (float(start)/1e6,float(end)/1e6)
  
  
  def getChromoContacts(self, chromo):
  
    nuc = self.mainApp.nuc
    
    if nuc and self.groupName:
      return self.chromoStats[chromo][0]
      
    
  def getExptRestr(self, chromo):
  
    nuc = self.mainApp.nuc
    
    if nuc:
      cis, trans = self.restraintStats.get(chromo, (0,0))
      n = self.restraintStats[True] or 1
      particGroup = nuc._getParticleGroup(self.structure)
      
      pt = 0.0
      
      for c in particGroup:
        pt += particGroup[c]['positions'].shape[0]
      
      p = particGroup[chromo]['positions'].shape[0]
      
      obs = (cis+trans) / float(n)
      
      exp = p/(pt or 1.0)
      
      return obs / exp
      


  def getChromoCis(self, chromo):
  
    nuc = self.mainApp.nuc
    
    if nuc:
      cis, trans = self.restraintStats.get(chromo, (0,0))
      return cis
  
  
  def getChromoTrans(self, chromo):
  
    nuc = self.mainApp.nuc

    if nuc:
      cis, trans = self.restraintStats.get(chromo, (0,0))
      return trans
  
  
  def selectCell(self, obj, row, col):
  
    pass
  
  
  def showEvent(self, event):
    
    self.updateContents()
    
    
  def changeStructure(self, code):
    
    self.structure = code
    self.restraintStats = self.getRestraintStats()
    self.chromoTable.setObjects(self.mainApp.nuc.getChromosomes())
 
 
  def updateContents(self):
  
    if self.isVisible():
      nuc = self.mainApp.nuc
      
      enabled = False  
      if self.groupName:
        group = nuc.getContactGroup(self.groupName)
        
        if group and group.parent == nuc.workContacts:
          enabled = True
        
      self.contactButtons.buttons[0].setEnabled(enabled)
                 
      codes, names = nuc.getStructureNames()      
      
      if codes:
        if self.structure in codes:
          index = codes.index(self.structure)
        
        elif nuc.structure.name in codes:
          self.structure = nuc.structure.name
          index = codes.index(self.structure)
         
        else:
          self.structure = codes[0]
          index = 0
      
      else:
        index = 0
      
      self.structPulldown.setData(names, codes, index)
      
      self.restraintStats = self.getRestraintStats()
      self.contactStats, self.chromoStats = self.getContactStats()
      
      if nuc:
        self.chromoTable.setObjects(nuc.getChromosomes())
 
      else:
        self.chromoTable.setObjects([])
 
 
      self.chromoTable.update()
      
      names = nuc.getContactGroupNames()
      self.numsTable.setObjects(names)
      self.numsTable.update()
  
  
  def getGroupName(self, group):
  
    return group
    
  
  def getAllCount(self, group):
    
    n, c, t, ic, it = self.contactStats.get(group, [0,0,0,0,0])
    
    return n


  def getCisCount(self, group):
    
    n, c, t, ic, it = self.contactStats.get(group, [0,0,0,0,0])
    
    return c
  
  
  def getTransCount(self, group):
  
    n, c, t, ic, it = self.contactStats.get(group, [0,0,0,0,0])
    
    return t
  
  
  def getTransPercent(self, group):
  
    n, c, t, ic, it = self.contactStats.get(group, [0,0,0,0,0])
    
    if n:
      return 100.0 * t/float(n)
  
  
  def getCisIsolated(self, group):
  
    n, c, t, ic, it = self.contactStats.get(group, [0,0,0,0,0])
    
    if ic and c:
      return 100.0 * ic/float(c)
  
  
  def getTransIsolated(self, group):
  
    n, c, t, ic, it = self.contactStats.get(group, [0,0,0,0,0])
    
    if it and t:
      return 100.0 * it/float(t)
  
  
  def getBinSize(self, group):
      
    nuc = self.mainApp.nuc
    binSize = nuc.getContactsBinSize(group)
    
    return int(binSize)
    
    
  def getContactStats(self):
  
    countDict = {}
    chromoStats = {}
    nuc = self.mainApp.nuc
    
    for groupName in nuc.getContactGroupNames():
      n0 = 0
      c0 = 0
      t0 = 0
 
      for chromo in nuc.getChromosomes():
        c = nuc.getNumContacts(groupName, [chromo,], trans=False)
        t = nuc.getNumContacts(groupName, [chromo,], cis=False)
        chromoStats[chromo] = [c+t, c, t]
        n0 += c + t
        c0 += c
        t0 += t
      
      t0 /= 2 # Each Trans contact has been counted both ways (though only stored one way)
      
      if groupName in nuc.origContacts:
        isSingleCell = nuc.origContacts[groupName].attrs['isSingleCell']
      
      elif groupName in nuc.workContacts:
        isSingleCell = nuc.workContacts[groupName].attrs['isSingleCell']
      
      else:
        isSingleCell = False
      
      if isSingleCell:
        ic0 = nuc.getNumIsolatedContacts(groupName, trans=False)
        it0 = nuc.getNumIsolatedContacts(groupName, cis=False)
      
      else:
        ic0 = None
        it0 = None
      
      countDict[groupName] = [n0, c0, t0, ic0, it0]
        
    return countDict, chromoStats
    
    
  def getRestraintStats(self):
     
    countDict = {} 
    nuc = self.mainApp.nuc
    
    m = 0
    
    if self.structure:
      c = 0
      t = 0
 
      restraintGroup = nuc.structures[self.structure]['restraints']
 
      for chrA in restraintGroup:
        subGroup = restraintGroup[chrA]
 
        for chrB in subGroup:
          n = subGroup[chrB].shape[1]
          m += n
 
          if chrA == chrB:
            if chrA in countDict:
              countDict[chrA][0] += n
 
            else:
              countDict[chrA] = [n, 0]
 
          else:
            if chrA in countDict:
              countDict[chrA][1] += n
 
            else:
              countDict[chrA] = [0, n]
 
            if chrB in countDict:
              countDict[chrB][1] += n
 
            else:
              countDict[chrB] = [0, n]
    
    countDict[True] = m
    
    return countDict
