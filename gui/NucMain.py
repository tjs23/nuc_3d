from PySide import QtCore, QtGui
from os import path, remove
from shutil import copy2, move
from hashlib import sha256
from h5py import File
from random import seed, randint
from numpy import array
import sys, os, numpy, time, re

seed(time.time())
sys.path.append('../')

Qt = QtCore.Qt
QLeftButton = Qt.LeftButton
QMiddleButton = Qt.MiddleButton
QRightButton = Qt.RightButton
QPoint = QtCore.QPoint
QKeys = QtGui.QKeySequence
QAction = QtGui.QAction 

from gui.qtgui.Base import Icon
from gui.qtgui.Button import Button
from gui.qtgui.Entry import IntRangesEntry, IntEntry, Entry
from gui.qtgui.FileSelect import selectFile, selectFiles, FileType, FileDialog, selectSaveFile
from gui.qtgui.Frame import Frame, FlowFrame
from gui.qtgui.InputDialog import askString, askInteger, askFloat, askChoice
from gui.qtgui.Label import Label
from gui.qtgui.LabelFrame import LabelFrame
from gui.qtgui.Menu import Menu
from gui.qtgui.MessageDialog import showWarning, showYesNo, showMulti
from gui.qtgui.MessageDialog import showSaveDiscardCancel, showInfo, showOkCancel
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.Slider import FloatSlider, Slider
from gui.qtgui.SpinBox import FloatSpinBox, IntSpinBox
from gui.qtgui.Splitter import Splitter
from gui.qtgui.TabbedFrame import TabbedFrame
from gui.qtgui.Text import Console
from gui.qtgui.ToolBar import ToolBar
from gui.qtgui.Colors import ColorDialog, GradientEditor

from analyses import Graphs

from gui.ContactMapPanel import ContactMapPanel
from gui.DataTrackPanel import DataTrackPanel
from gui.FileSystemPanel import FileSystemPanel
from gui.GenomeBrowserPanel import GenomeBrowserPanel
from gui.ImagePanel import ImagePanel
from gui.InteractomePanel import InteractomeOuterPanel
from gui.OverviewPanel import OverviewPanel
from gui.LabelsPanel import LabelsPanel
from gui.RoiPanel import RoiPanel
from gui.StructureCalcPanel import StructureCalcPanel
from gui.StructureReportPanel import StructureReportPanel
from gui.StructurePanel import StructureOuterPanel

from formats.Util import getFileType, splitExtension
from formats.Util import STRUCTURE_FORMATS, CONTACT_FORMATS, DATA_TRACK_FORMATS, INTERACTIONS_FORMATS

from NucApi import DISPLAY_MODES, RESTRAINT_COLOR_MODES
from NucApi import DATA_TRACK_SYMBOLS, COLOR_MODES, DATA_TRACK_PEAK_TYPES
from NucApi import Nucleus, FILE_EXT, PROGRAM, EXTERNAL, INNATE, DERIVED, INTERACTIONS

from util.Io import checkRegularFile, pathExists
from util.Structure import getArcChromoCoords, getLinearChromoCoords

ICON_DIR =  path.join(path.dirname(path.dirname(__file__)), 'gui', 'icons')

ORD_A = ord('A')
SESSION_KEY = ''.join([chr(ORD_A+randint(0, 25)) for x in range(10)])
NEW_FILE_TEMP = '_Nuc3dTemp_%s.nuc.temp' % SESSION_KEY
MOVIE_FRAME_SIZES = ((1920,1080), # 1080p
                     (1280,1024),
                     (1280,720), # 720p
                     (1024,768),
                     (800,600))


class NucMain(QtGui.QMainWindow):

  def __init__(self, parent=None, filePath=None):
  
    QtGui.QMainWindow.__init__(self, parent)
    
    self.nuc = None
    self.dirPathNuc = '.'
    self.dirPathRef = '.'
    self.dirPathCont = self._getHomeDir() or '.'
    self.dirPathData = self.dirPathCont
    self.dirPathImage = self.dirPathCont
    self.dirPathMovie = self._getHomeDir() or '.'
    self.recentFiles = []
    
    self.coordImage = None # Eventually rplace with nuc.currentCoordImage etc.

    if filePath:
      filePath = path.abspath(filePath)
    
    self.inFilePath = filePath
    self.glLists = 0
    
    # Middle panels - pre linking
    
    self.structureOuterPanel = StructureOuterPanel(self, self, self.openFiles)
    self.structurePanel  = self.structureOuterPanel.innerPanel
    
    self.contactMapPanel = ContactMapPanel(self, self, self.openFiles)
    
    self.interactomeOuterPanel = InteractomeOuterPanel(self, self, self.openFiles)
    self.interactomePanel = self.interactomeOuterPanel.innerPanel
    
    self.structureReportPanel = StructureReportPanel(self, self)
    
    self.dataTrackPanel = DataTrackPanel(self, self)
    #self.imagePanel = ImagePanel(self)
    
    
    # Main Menu
       
    menuBar = self.menuBar()
    fileMenu = Menu(menuBar, '&File')
    contactMenu = Menu(menuBar, '&Contacts')
    trackMenu = Menu(menuBar, '&Data tracks')
    structMenu = Menu(menuBar, '&Structures')
    configMenu = Menu(menuBar, 'Con&figure')
    imageMenu = Menu(menuBar, '&Microscopy')
    helpMenu = Menu(menuBar, '&Help')
    
    
    # File Menu
    
    self.recentFilesMenu = Menu(fileMenu, 'Recent files',
                                setupFunc=self.setRecentFilesMenu)
    
    commands = (('&Open %s file' % FILE_EXT, self.openNucFile, QKeys.Open, 'open.png'),
                ('&Save %s file' % FILE_EXT, self.saveNucFile, QKeys.Save, 'save.png'),
                ('&Save as...',    self.saveNucFileAs, QKeys.SaveAs, None),
                ('&Set experiment reference file', self.setExperimentRef, QKeys("Ctrl+E"), None),
                ('&Set genome reference file', self.setGenomeRef, QKeys("Ctrl+G"), None),
                ('&Quit program',  self.close, QKeys("Ctrl+Q"), None))
    
    for text, cmd, keys, icon in commands:
      fileMenu.addItem(text, callback=cmd, shortcut=keys, icon=self.getIcon(icon))
    
    
    # Contacts Menu
    
    self.exportMatrixAction = None
    contactMenu.addItem('&Import contacts', self.importContacts)
    
    exportMenu = Menu(contactMenu, '&Export contacts',
                      setupFunc=self._setupExportContactsMenu)
 
    contactMenu.addItem('Export &image', self.exportContactImage)

    removeMenu = Menu(contactMenu, '&Remove')
    removeMenu.addItem('&Isolated', self.removeIsolatedContacts)
    removeMenu.addItem('&Structure violated', self.removeViolatedContacts)
    removeMenu.addItem('&Resolve ambiguous', self.resolveAmbigousContacts)
    
    createMenu = Menu(contactMenu, 'Create from...')
    createMenu.addItem('Structure &neighbours', self.generateRandomContacts)
    createMenu.addItem('Structure &violations', self.extractViolatedContacts)
    createMenu.addItem('&Merged contact groups', self.mergeContacts)
    createMenu.addItem('&Observed vs Expected', self.calcObsVsExpContacts)
     
    analysisMenu = Menu(contactMenu, 'Analyses')
    analysisMenu.addItem('&Calculate correlations', self.calcContactCorrelations)
    analysisMenu.addItem('&Pseudo-4C', self._graphPseudo4C)
    analysisMenu.addItem('&Sequence separation', self._graphContactSeqSep)
    analysisMenu.addItem('&Trans spacing', self._graphTransSpacing)
    analysisMenu.addItem('&Trans FFT', self._graphContactFourierTransform)

    contactMenu.addItem('&Normalise population data', self.normalisePopContacts)
     
    # Data Tracks Menu
    
    trackMenu.addItem('&Import data track', callback=self.importDataTrack)
    trackMenu.addItem('&Import genome features', callback=self.importGenomeFeature)
    trackMenu.addItem('&Import interactions', callback=self.importInteractions)

    exportMenu = Menu(trackMenu, 'Export data track',
                      setupFunc=self._setupExportDataMenu)

    exportIntMenu = Menu(trackMenu, 'Export interactions',
                         setupFunc=self._setupExportInteractionsMenu)

    refMenu = Menu(trackMenu, 'Reference data')
    refMenu.addItem('Add to &experimental reference', self.importExperimentRefData)
    refMenu.addItem('Add to &genome reference', self.importGenomeRefData)
      
    makeTrackMenu = Menu(trackMenu, 'Create data track from...')
    
    cisMenu   = Menu(makeTrackMenu, 'Contacts', setupFunc=self._setupCreateContactDataMenu)
    cisMenu   = Menu(makeTrackMenu, 'Cis contacts', setupFunc=self._setupCreateCisDataMenu)
    transMenu = Menu(makeTrackMenu, 'Trans contacts', setupFunc=self._setupCreateTransDataMenu)
    voidMenu  = Menu(makeTrackMenu, 'Void contact regions', setupFunc=self._setupCreateVoidDataMenu)
    densMenu  = Menu(makeTrackMenu, 'Data density', setupFunc=self._setupCreateDensityDataMenu)
        
    commands = (('Coordinate density', self.createCoordDensityDataTrack, None),
                ('Nucleus depth', self.createDepthDataTrack, None),
                ('Chromosome depth', self.createChromoDepthDataTrack, None),
                ('Tans interface distance', self.createTransDepthDataTrack, None),
                ('Radius of gyration', self.createRadGyrationDataTrack, None),
                ('Intermingled territories', self.createIntermingleDataTrack, None),
                ('Contact distance', self.createContactDistDataTrack, None),
                ('Recording', self.createRecordingDataTrack, None),
                ('Regional structure violations', self.createRegionContactDistDataTrack, None),
                ('Backbone distance', self.createBackboneDictDataTrack, None),
                ('Backbone angle',  self.createLinearityDataTrack, None))   
    
    for text, cmd, keys in commands:
      makeTrackMenu.addItem(text, callback=cmd, shortcut=keys)

    makeInteractMenu = Menu(trackMenu, 'Create interactions from...')    
    makeInteractMenu.addItem('Data track pairs', callback=self.createDataTrackPairInteractions, shortcut=None)
      
      
    # Structure Menu
    
    currentStructMenu = Menu(structMenu, 'Current structure',
                             setupFunc=self._setupCurrentStructureMenu)

    structMenu.addItem('&Import structure', self.importCoords)
    
    exportMenu = Menu(structMenu, 'Export structure',
                      setupFunc=self._setupExportStructureMenu)
    
    structMenu.addItem('Export &image', self.exportStructureImage)
   
    movieMenu = Menu(structMenu, 'Export movie')
    self.movieMenu = movieMenu

    movieMenu.addItem('Annealing', self.exportAnnealMovieImages)
    #movieMenu.addItem('Clipping', self.exportChromosomeMovieImages)
    movieMenu.addItem('Expansion', self.exportExpansionMovieImages)
    movieMenu.addItem('Rotation', self.exportRotationMovieImages)

    movieMenu.addSeparator()
    
    for i, size in enumerate(MOVIE_FRAME_SIZES):
      movieMenu.addItem('%d x %d' % size, object=size, checked=i==2, group=0)
    
    self.exportCurrentModelsAction = None
    analysisMenu = Menu(structMenu, 'Analyses')
    analysisMenu.addItem('&Contact distances', self._graphContactDistances)
    analysisMenu.addItem('&Violation', self._graphViolations)
    analysisMenu.addItem('&Chromosome intermingling', self._graphIntermingling)
    analysisMenu.addItem('&Hierarchical clustering', self._clusterStructures)
    analysisMenu.addItem('&Data track 3D clustering', self._clusterDataTrack3d)

    calcMenu = Menu(structMenu, 'Calculations')
    calcMenu.addItem('&Anneal structure',self.showStrucCalc)
    calcMenu.addItem('(Re)calculate &restraints',self.recalcRestraints)
    calcMenu.addItem('(Re)calculate &densities', self.recalcDensity)
    
    transformMenu = Menu(structMenu, 'Synthetic coords')
    transformMenu.addItem('Random walk', self.setCoordsRandomWalk)
    transformMenu.addItem('Linear stack', self.setCoordsLinearStack)
    transformMenu.addItem('Great circle', self.setCoordsGreatCircle)
      
    transformMenu = Menu(structMenu, 'Transform coords')
    transformMenu.addItem('Align models', self.alignModels)
    transformMenu.addItem('Center structure',self.centerStructure)
    transformMenu.addItem('Mirror structure',self.mirrorStructure)
      
    # Configure Menu
    
    panel = self.structureOuterPanel
    
    colorMenu = Menu(configMenu, '&Colour mode', icon=self.getIcon('paint.png'))
    colorMenu.addItem('Sequence postion', panel.colorSeq,
                      icon=self.getIcon('color-seq.png'),
                      tipText=None)
    colorMenu.addItem('Chromosome', panel.colorChromo,
                      icon=self.getIcon('color-chromo.png'),
                      tipText=None)
    colorMenu.addItem('Data tracks', panel.colorGenData,
                      icon=self.getIcon('color-data.png'),
                      tipText=None)
    #colorMenu.addItem('Faint', panel.colorFaint,
    #                  icon=self.getIcon('color-faint.png'),
    #                  tipText=None)
    #colorMenu.addItem('Region of interest', panel.displayRoi,
    #                  icon=self.getIcon('color-roi.png'),
    #                  tipText=None)
    colorMenu.addItem('Model number', panel.colorModel,
                      icon=self.getIcon('color-model.png'),
                      tipText=None)
    colorMenu.addItem('RMSD', panel.colorRmsd,
                      icon=self.getIcon('rmsd.png'),
                      tipText=None)
    
    colorOptMenu = Menu(configMenu, '&Colour options',
                        icon=self.getIcon('colors.png'),
                        setupFunc=panel._setupColorConfigMenu)
   
    dispMenu = Menu(configMenu, '&Render style', icon=self.getIcon('render.png'))
    dispMenu.addItem('Line', panel.displayLine,
                     icon=self.getIcon('display-line.png'),
                     tipText=None)
    dispMenu.addItem('Tube', panel.displayTube,
                     icon=self.getIcon('display-tube.png'),
                     tipText=None)
    dispMenu.addItem('Ball and stick',  panel.displayBallStick,
                     icon=self.getIcon('display-ball.png'),
                     tipText=None)
    self.restActionA = dispMenu.addItem('Show restraints', panel.toggleRestraints,
                                        checked=False,  tipText=None)
    Menu(configMenu, '&Render options', icon=self.getIcon('configure.png'),
         setupFunc=self._setupDispDetailsMenu)
    
    Menu(configMenu, '&Mouse motion',
         setupFunc=self._setupMouseMotionMenu, icon=self.getIcon('mouse.png'))
    
      
    # Image Menu

    imageMenu.addItem('&Import point cloud', self.importImageCoords)
    imageMenu.addItem('&Center point clouds', self.centreImageCoords)
    
    coordImageMenu = Menu(imageMenu, 'Point cloud image',
                          setupFunc=self._setupCoordImageMenu)
      
    # Help Menu
    
    commands = (('&About %s' % PROGRAM, self._notImplemented, None),
                ('&Documentation', self._notImplemented, None),
                ('&Tutorials',  self._notImplemented, None),
                ('&Reset session settings', self._resetQtSettings, None))
    
    for text, cmd, keys in commands:
      helpMenu.addItem(text, callback=cmd, shortcut=keys)
    
      

    # Toolbars
    
    # File toolbar
    
    icons = ['new.png', 'open.png', 'save.png']
    icons = [self.getIcon(i) for i in icons]
    texts = ['&New %s' % FILE_EXT, '&Load %s' % FILE_EXT, '&Save %s' % FILE_EXT,]
    shortcuts = [QKeys.New, QKeys.Open, QKeys.Save]
    funcs = [self.newNucFile, self.openNucFile,
             self.saveNucFile]
    
    self.fileToolbar = ToolBar(self, 'File toolbar', funcs, icons, texts, 
                               shortcuts, objName='fileToolbar',
                               iconSize=32, areas='tblr')
    self.addToolBar(Qt.LeftToolBarArea, self.fileToolbar)
    
    # Side panels toolbar
    
    icons = ['file-system.png', 'struc-calc.png', 'stats.png',
             'color-roi.png', 'labels.png', 'command_line.png']
    icons = [self.getIcon(i) for i in icons]
    texts = ['Open file browser', 'Structure calculations', 'Stats overview',
             'Regions of interest', 'Positional labels', 'Python console']
    shortcuts = [None] * len(texts)
    funcs = [self.toggleFileBrowser, self.toggleStrucCalc,
             self.toggleStatsOverview, self.toggleRoiTable,
             self.toggleLabelsTable, self.toggleConsole]
    
    self.panelToolbar = ToolBar(self, 'Side panels toolbar', funcs, icons, texts, 
                                shortcuts, objName='panelToolbar', areas='l',
                                movable=False, iconSize=32)
    self.addToolBar(Qt.LeftToolBarArea, self.panelToolbar)
        
    
     # Chromosome toolbar
    
    self.chromoActions = []
    bg = self.chromoButtonGroup = QtGui.QButtonGroup(self)
    bg.setExclusive(False)
    bg.connect(bg, QtCore.SIGNAL('buttonClicked(int)'),  self.toggleChromosomes)
    self.chromoToolbar = ToolBar(self, 'Chromosome toolbar', [], None, None, 
                                 objName='chromoToolbar', areas='tblr', iconSize=(32,32))
    self.addToolBar(Qt.RightToolBarArea, self.chromoToolbar)
                                  
    callbacks = [self.showNoChromos, self.showAllChromos]
    icons = ['chromo-none.png', 'chromo-all.png']
    icons = [self.getIcon(i) for i in icons]
    texts = ['Show no chromosomes', 'Show all chromosomes']                              
    self.chromoToolbar.setActions(callbacks, icons, texts, shortcuts=None)
    
    configChromosButton = Button(self.chromoToolbar, ' ',
                                 icon=self.getIcon('configure.png'),
                                 iconSize=32)
    configChromosMenu = Menu(configChromosButton, 'Configure chromosome options')
    
    configChromosMenu.addItem('Invert selection', callback=self.chromosInvertSelected)
    configChromosMenu.addItem('Show unselected as faint', callback=self.chromosFaintUnselected)
    
    colorChromosMenu = Menu(configChromosMenu, 'Reset chromosome colours')
    
    colorChromosMenu.addItem('Gradient', callback=self.chromosResetColors)
    colorChromosMenu.addItem('Cycle',    callback=self.chromosResetCyclingColors)
    
    configChromosButton.setMenu(configChromosMenu)
    self.chromoToolbar.addWidget(configChromosButton)

   
    """
    
    gDataTypes = [EXTERNAL, DERIVED, INNATE, INTERACTIONS]
    titles = ['Experimental', 'Derived', 'Genome', 'Interactions']
    
    self.gDataButtonGroups = {}
    self.gDataActions = {}
    self.gDataToolbars = {}
    for i, typ in enumerate(gDataTypes):
      bg =  QtGui.QButtonGroup(self)
      bg.setExclusive(False)
      bg.connect(bg, QtCore.SIGNAL('buttonClicked(int)'), lambda x, t=typ:self.toggleDataTrack(x, t))
      self.gDataButtonGroups[typ] = bg
      self.gDataActions[typ] = []
      name = titles[i]
      tb = ToolBar(self, '%s data tracks toolbar' % name, [], None, None, 
                   objName='dataTracksToolbar%s' % name, areas='tbl')
      tb.addWidget(Label(tb, name+':'))
      self.gDataToolbars[typ] = tb
      self.addToolBar(Qt.BottomToolBarArea, tb)
    """
     

    # Display toolbar
    # Main Layout
    
    middle = Frame(self)
    
    self.splitter = Splitter(middle, grid=(0,0), stretch=(1,1))
    self.splitter.setObjectName('Splitter')
    
    
    # Data tracks toolbar
    
    self.dataTrackButtonGroup = QtGui.QButtonGroup(self)
    self.dataTrackButtonGroup.setExclusive(False)
    self.dataTrackButtonGroup.connect(self.dataTrackButtonGroup, QtCore.SIGNAL('buttonClicked(int)'), lambda x:self.toggleDataTrack(x))
    
    #tb = ToolBar(self, 'Data tracks toolbar', [], None, None, 
    #             objName='dataTracksToolbar', areas='tbl')
    #tb.addWidget(Label(tb, name+':'))
    
    self.dataTrackFrame = FlowFrame(middle, grid=(1,0))
    
    self.setCentralWidget(middle)
    
    # Left panels
    
    self.leftFrame = QtGui.QWidget()
    self.fileSystemPanel = FileSystemPanel(self.leftFrame, self.openFiles)
    self.strucCalcPanel = StructureCalcPanel(self.leftFrame, self,
                                             self.structureOuterPanel.refreshStructure)
    self.overviewPanel = OverviewPanel(self.leftFrame, self)
    self.roiPanel = RoiPanel(self.leftFrame, self)
    self.labelsPanel = LabelsPanel(self.leftFrame, self)
   
    msg = 'Nuc API interactive console\nQt window available as "gui"\nHDF5 file available as "root"'
    
    self.consolePanel = Console(self.leftFrame, msg, self.nuc, closeFunc=self.toggleConsole)

    self.splitter.addWidget(self.leftFrame)
    self.leftFrame.hide()
    self.fileSystemPanel.hide()
    self.strucCalcPanel.hide()
    
    
    # Middle panels
    
    
    texts = ['Structure display', 'Contact map', 'Interactome',
             'Structure report', 'Data track report'] # , 'Images']
    icons = None
    widgets = [self.structureOuterPanel, self.contactMapPanel,
               self.interactomeOuterPanel, self.structureReportPanel,
               self.dataTrackPanel] #, self.imagePanel]
    self.tabbedPanel = TabbedFrame(self, texts, icons,
                                   self.selectTab, widgets)

    self.splitter.addWidget(self.tabbedPanel)
    
    self.defaultState = self.saveState()
    self.defaultGeometry = self.saveGeometry()
    
    # Finally
    self._readQtSettings()


  def _isTempFile(self, filePath):
  
    return filePath.endswith('.nuc.temp')
  
  
  def _getTempFilePath(self, filePath):
    
    
    if self._isTempFile(filePath):
      return filePath
      
    else:
      root, fex = path.splitext(filePath)
      ending = '_%s.nuc.temp' % SESSION_KEY
    
      return '%s%s' % (root, ending)
    
  
  def _getPersistantFilePath(self, tempPath):
  
    ending = '_%s.nuc.temp' % SESSION_KEY
    
    if self._isTempFile(tempPath):
      return tempPath[:-len(ending)] + '.nuc'
      
    else:
      return tempPath
  
  
  def _isFileInUse(self, filePath):


    dirName, fileName = path.split(filePath)
    
    tempFileName = self._getTempFilePath(fileName)
    permFileName = self._getPersistantFilePath(fileName)
    
    for fileName2 in os.listdir(dirName):
      if fileName2.endswith('.nuc.temp'):
        if fileName2 != tempFileName:
          if self._getPersistantFilePath(fileName2) == permFileName:
            return fileName2
    
    return False
    
  
  def _setChromoColor(self, chromo):
    
    colorObj = ColorDialog(self).getColor(self.nuc.getChromoColor(chromo).rgbHex)
            
    if colorObj:
      self.nuc.setChromoColor(chromo, colorObj.getRgbF())
      self.updateContents()
  
  
  def _setChromoColorMode(self, obj):
  
    chromo, mode = obj
    if self.nuc.setChromoDisplayParams(chromo, colorMode=mode):
      self.updateContents()


  def _setChromoDisplayMode(self, obj):
    
    chromo, mode = obj 
    if self.nuc.setChromoDisplayParams(chromo, displayMode=mode):
      self.updateContents()
  
  
  def _setChromoRoiAll(self, chromo):
  
    start, end = self.nuc.getChromosomeLimits(chromo)
    
    if self.nuc.setChromoRegionOfInterest(chromo, start, end):
      self.updateContents()
      

  def _setChromoRoiStart(self, val, chromo):
  
    start, end = self.nuc.getChromoRegionOfInterest(chromo)
    
    if self.nuc.setChromoRegionOfInterest(chromo, int(1e6*val), end):
      self.updateContents()


  def _setChromoRoiEnd(self, val, chromo):
  
    start, end = self.nuc.getChromoRegionOfInterest(chromo)
    
    if self.nuc.setChromoRegionOfInterest(chromo, start, int(1e6*val)):
      self.updateContents()
  
  
  def _setChromoRoiPosition(self, frac, chromo, start_entry, end_entry):
     
    start, end = self.nuc.getChromoRegionOfInterest(chromo)
    delta = end - start
    
    p1, p2 = self.nuc.getChromosomeLimits(chromo)
    val = p1 + (p2-p1) * (frac/100.0)
    
    start = max(val-delta/2, p1)
    end   = min(val+delta/2, p2)
    
    if self.nuc.setChromoRegionOfInterest(chromo, start, end):
      start_entry.set(start/1e6, doCallback=False)
      end_entry.set(end/1e6, doCallback=False)
      self.updateContents()
  
  
  def _set_roi_slider(self, chromo, slider, start=None, end=None):
    
    first, last = self.nuc.getChromosomeLimits(chromo)
    p1, p2 = self.nuc.getChromoRegionOfInterest(chromo)
    
    if start is None:
      start = p1
    else:
      start = start.get()
      
    if end is None:
      end = p2  
    else:
      end = end.get()
    
    frac = 50 * (start+end) / (last-first)
    
    slider.set(frac, doCallback=False)
    
  
  def _setChromoHighlight(self, selected):
    
    setDisplay =  self.nuc.setChromoDisplayParams
    
    for chromo in self.nuc.getChromosomes():
      if chromo == selected:
        # Chromo coloured tubes
        setDisplay(chromo, isShown=True, useLabels=False, colorMode=1, displayMode=2)
      else:
        # Faint lines
        setDisplay(chromo, isShown=True, useLabels=False, colorMode=5, displayMode=1)
        
    self.updateContents()
  

  def _setDataColor(self, obj):
    
    typ, code = obj
    
    prev = self.nuc.getDataTrackColor(typ, code).rgbHex
    colorObj = ColorDialog(self).getColor(prev)
            
    if colorObj:
      if self.nuc.setDataTrackColor(typ, code,  colorObj.getRgbF()):
        self.updateContents()


  def _setInteractionsColor(self, code):

    prev = self.nuc.getInteractionsColor(code).rgbHex
    colorObj = ColorDialog(self).getColor(prev)
            
    if colorObj:
      if self.nuc.setInteractionsColor(code,  colorObj.getRgbF()):
        self.updateContents()
  
 
  def _setDataSymbol(self, obj):
  
    typ, code, index = obj
    
    if self.nuc.setDataTrackSymbol(typ, code,  index):
      self.updateContents()


  def _setDataTrackType(self, obj):
  
    typ, code, index = obj
    
    if self.nuc.setDataTrackPeakType(typ, code, index):
      self.updateContents()
  
      
  def _setDataScale(self, val, typ, code):
   
    if self.nuc.setDataTrackScale(typ, code, val):
      self.updateContents()
   
   
  def _setDataMinValue(self, val, typ, code):

    if self.nuc.setDataTrackThresholds(typ, code, minVal=val):
      self.updateContents()


  def _setDataMaxValue(self, val, typ, code):

    if self.nuc.setDataTrackThresholds(typ, code, maxVal=val):
      self.updateContents()
 
  def _setDataScale(self, val, typ, code):
   
    if self.nuc.setDataTrackScale(typ, code, val):
      self.updateContents()
  
  def _setDataTrackDetails(self, val, typ, code):
   
    if self.nuc.setDataTrackDetails(typ, code, val):
      self.updateContents()

  def _setInteractionsDetails(self, val, code):
   
    if self.nuc.setInteractionsDetails(code, val):
      self.updateContents()
  
  def _setInteractionsLineWidth(self, val, code):
  
    if self.nuc.setInteractionsStyle(code, lineWidth=val):
      self.updateContents()
  
  
  def _setInteractionsMinValue(self, val, code):
  
    if self.nuc.setInteractionsThresholds(code, minVal=val):
      self.updateContents()
  
  
  def _setupDataTrackContextMenu(self, pos, code):
     
     menu = Menu(self, 'Data track context menu')
     
     if code in self.nuc.interactions:
       menu.addItem('Interactions "%s"' % code, callback=None, widget=None).setEnabled(False)
       menu.addSeparator()
 
       menu.addItem('Set display color', callback=self._setInteractionsColor,
                    icon=self.getIcon('colors.png'), object=code)
 
       # Could set the style here...
       
       # dataLayer.attrs['options'] = (0, showText, style, null)
       # dataLayer.attrs['display'] = (r, g, b, a, lineWidth, minVal, maxVal)
 
       menu.addSeparator()
       attrs = self.nuc.interactions[code].attrs
       
       val = attrs['details']
       entry = Entry(menu, val,  grid=None,
                     callback=lambda v, c=code:self._setInteractionsDetails(v, c))
       menu.addItem('Details', widget=entry)
              
       menu.addSeparator()

       val = attrs['display'][4] # lineWidth
       entry = FloatSpinBox(menu, val, 0.0, step=0.1, grid=None,
                            callback=lambda v, c=code:self._setInteractionsLineWidth(v, c))
       menu.addItem('Line width', widget=entry)
 
       val = attrs['display'][5] # min value
       entry = FloatSpinBox(menu, val, 0.0, step=0.1, grid=None,
                            callback=lambda v, c=code: self._setInteractionsMinValue(v, c))
       entry.setDecimals(2)
       menu.addItem('Min. value', widget=entry)
 
       menu.addSeparator()
       menu.addItem('Export interactions', callback=self.exportInteractions, object=code)

       menu.addItem('Remove interactions', callback=self.removeInteractions, object=code)
     
     else:
       for source in EXTERNAL, INNATE, DERIVED:
         if code in self.nuc.dataTracks[source]:
           break
       
       currentSymbol = self.nuc.getDataTrackSymbol(source, code)
     
       obj = (source, code)
       menu.addItem('Data track "%s"' % code, callback=None, widget=None).setEnabled(False)
       menu.addSeparator()
 
       menu.addItem('Set display color', callback=self._setDataColor,
                    icon=self.getIcon('colors.png'), object=obj)
                    
       symbolMenu = Menu(menu, 'Structure symbol')
       for i, text in enumerate(DATA_TRACK_SYMBOLS):
         symbolMenu.addItem(text, group=1, checked=i == currentSymbol,
                            object=(source, code, i), callback=self._setDataSymbol)
 
       currentTrackType = self.nuc.getDataTrackPeakType(source, code)
       trackTypeMenu = Menu(menu, 'Data track type')
       for i, text in enumerate(DATA_TRACK_PEAK_TYPES):
         trackTypeMenu.addItem(text, group=1, checked=i == currentTrackType,
                               object=(source, code, i), callback=self._setDataTrackType)
 
       # dataLayer.attrs['options'] = (0, showText, shape, 0, 0, 0, 0, 0)
       # dataLayer.attrs['display'] = (r, g, b, 0.0, scale, threshold, 0.0, 0.0)
       attrs = self.nuc.dataTracks[source][code].attrs
 
       menu.addSeparator()
       
       if 'details' in attrs:
         val = attrs['details']
       else:
         val = code
         attrs['details'] = val
       
       entry = Entry(menu, val,  grid=None,
                     callback=lambda v, s=source, c=code:self._setDataTrackDetails(v, s, c))
       menu.addItem('Details', widget=entry)
              
       menu.addSeparator()
 
       val = attrs['display'][4] # Scale
       entry = FloatSpinBox(menu, val, 0.0, step=0.1, grid=None,
                            callback=lambda v, s=source, c=code:self._setDataScale(v, s, c))
       entry.setDecimals(3)
       menu.addItem('Size scale', widget=entry)
 
       val = attrs['display'][5] # min value
       entry = FloatSpinBox(menu, val, 0.0, step=0.1, grid=None,
                            callback=lambda v, s=source, c=code: self._setDataMinValue(v, s, c))
       entry.setDecimals(3)
       menu.addItem('Min. value', widget=entry)
 
       #val = attrs['display'][6] # max value
       #entry = FloatSpinBox(menu, val, 0.0, step=0.1, grid=None,
       #                     callback=lambda v, s=source, c=code: self._setDataMaxValue(v, s, c))
       #entry.setDecimals(2)
       #menu.addItem('Max. value', widget=entry)
 
       menu.addSeparator()
       menu.addItem('Export data', callback=self.exportDataTrack,
                    object=obj)

       menu.addItem('Remove dataset', callback=self.removeDataTrack,
                    object=obj)
 
     menu.popup(pos)

   
  def contextMenuEvent(self, event):

    widget = self.childAt(event.pos())
  
    if self.nuc and widget:
      if isinstance(widget, QtGui.QPushButton):
        buttonGroup = widget.group()
        
        if buttonGroup is self.chromoButtonGroup:
          name = widget.text()
          name = name.strip()
          menu = Menu(self, 'Chromosome context menu')
 
          menu.addItem('', widget=LabelFrame(menu, 'Chromosome %s' % name))
          menu.addSeparator()
 
          menu.addItem('Set ID color', callback=self._setChromoColor,
                       icon=self.getIcon('colors.png'), object=name)
 
          currentMode = self.nuc.chromosomes[name].attrs['display'][2]
          colorMenu = Menu(menu, 'Color mode')
          for i, text in enumerate(COLOR_MODES):
            colorMenu.addItem(text, group=1, checked=i == currentMode,
                              object=(name, i), callback=self._setChromoColorMode)
 
          currentMode = self.nuc.chromosomes[name].attrs['display'][3]
          dispMenu = Menu(menu, 'Render style')
          for i, text in enumerate(DISPLAY_MODES):
             dispMenu.addItem(text, group=1, checked=i == currentMode,
                              object=(name, i), callback=self._setChromoDisplayMode)
 
          menu.addItem('Highlight Chr. %s' % name, object=name, callback=self._setChromoHighlight)      
        
          menu.popup(self.mapToGlobal(event.pos()))
          return
 
        elif buttonGroup is self.dataTrackButtonGroup:
          code = widget.obj
          self._setupDataTrackContextMenu(self.mapToGlobal(event.pos()), code)
 
          return
          
      elif isinstance(widget, GenomeBrowserPanel):
        gPos = self.mapToGlobal(event.pos())
        pos = widget.mapFromGlobal(gPos)
        key = widget._getScreenDataTrack(pos.x(), pos.y())
        
        if key is not None:
          typ, code = key
          self._setupDataTrackContextMenu(gPos, code)
          return
                     
    QtGui.QMainWindow.contextMenuEvent(self, event)
  
  
  def showEvent(self, event):
    
    QtGui.QMainWindow.showEvent(self, event)
    
    if not self.nuc:
      if self.inFilePath:
        self.openNucFile(self.inFilePath)
      
      """
      elif self.recentFiles:
        msg = 'Select %s file' % PROGRAM
        objs = self.recentFiles + [None]
        texts = self.recentFiles + ['New file']
 
        filePath = showMulti('Query', msg, texts, objs, self)
 
        if filePath:
          self.openNucFile(filePath)
      """
 
      if not self.nuc:
        self.newNucFile()
  
  
  def _setMouseMotionFunc(self, opt):
  
    isRotation, moveFunc = opt
    
    for panel in [self.structurePanel, self.interactomePanel]:
    
      if isRotation:
        setattr(panel, moveFunc, panel._mouseRotate)
      else:
        setattr(panel, moveFunc, panel._mouseTranslate)


  def _setupCoordImageMenu(self, menu):
  
    menu.clear()
  
    codes = self.nuc.getCoordImages()
    
    menu.addItem('None', self.setCoordImage, False, group=0, checked=not self.coordImage)
    
    for code in sorted(codes):
      menu.addItem(code, self.setCoordImage, code, group=0, checked=self.coordImage == code)
  

  def _setupCurrentStructureMenu(self, menu):
  
    menu.clear()
  
    nuc = self.nuc
    
    for code, name in zip(*nuc.getStructureNames()):
      menu.addItem(name, self.setCurrentStructure, code)
  
  
  def _setupCreateCisDataMenu(self, menu):
  
    menu.clear()
  
    nuc = self.nuc
    for groupName in nuc.getContactGroupNames():
      menu.addItem(groupName, self.createCisDataTrack, groupName)


  def _setupCreateContactDataMenu(self, menu):
  
    menu.clear()
  
    nuc = self.nuc
    for groupName in nuc.getContactGroupNames():
      menu.addItem(groupName, self.createContactDataTrack, groupName)
     
  
  def _setupCreateTransDataMenu(self, menu):
  
    menu.clear()
  
    nuc = self.nuc
    for groupName in nuc.getContactGroupNames():
      menu.addItem(groupName, self.createTransDataTrack, groupName)
  
  
  def _setupCreateVoidDataMenu(self, menu):
  
    menu.clear()
  
    nuc = self.nuc
    for groupName in nuc.getContactGroupNames():
      menu.addItem(groupName, self.createVoidDataTrack, groupName)
  
  
  def _setupMouseMotionMenu(self, menu):
  
    menu.clear()
    
    panel = self.structurePanel
    lbr = panel.mouseLeftMoveFunc == panel._mouseRotate
    mbr = panel.mouseMiddleMoveFunc == panel._mouseRotate
    rbr = panel.mouseRightMoveFunc == panel._mouseRotate
    
    menu.addItem('Left button rotates', self._setMouseMotionFunc,
                 (True, 'mouseLeftMoveFunc'),
                 checked=lbr, group=0)
    menu.addItem('Left button moves',  self._setMouseMotionFunc,
                 (False, 'mouseLeftMoveFunc'),
                 checked=not lbr, group=0) 
     
    menu.addItem('Middle button rotates', self._setMouseMotionFunc,
                 (True, 'mouseMiddleMoveFunc'),
                 checked=mbr, group=1)
    menu.addItem('Middle button moves', self._setMouseMotionFunc,
                 (False, 'mouseMiddleMoveFunc'),
                 checked=not mbr, group=1)
    
    menu.addItem('Right button rotates', self._setMouseMotionFunc,
                 (True, 'mouseRightMoveFunc'),
                 checked=rbr, group=2)
    menu.addItem('Right button moves', self._setMouseMotionFunc,
                 (False, 'mouseRightMoveFunc'),
                 checked=not rbr, group=2)
  
  
  def setGradientColors(self, colors, name,):
    
    rgbs = [[int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)] for c in colors]
    rgbs = (numpy.array(rgbs) / 255.0).tolist()
    
    group = self.nuc.display['colorSchemes']
    schemeAttrs = self.nuc._setAttr(group, name, rgbs)
    
    self.updateContents()    
    

  def _setDensityRadius(self, value):
        
    nuc = self.nuc
    nuc.setDisplaySizes(density=value)


  def setSphereDetail(self, value):
    
    if self.nuc.setDisplayDetailLevels(ballDetail=value):
      self.updateContents()


  def setLineSmooth(self, value):
  
    if self.nuc.setDisplayDetailLevels(lineSmooth=value):
      self.updateContents()


  def setTubeSmooth(self, value):
  
    if self.nuc.setDisplayDetailLevels(tubeSmooth=value):
      self.updateContents()
  
  
  def setTubeDetail(self, value):
  
    if self.nuc.setDisplayDetailLevels(tubeDetail=value):
      self.updateContents()

            
  def setSphereSize(self, radius):
  
    if self.nuc.setDisplaySizes(ball=radius):
      self.updateContents()


  def setStickSize(self, frac):
       
    if self.nuc.setDisplaySizes(stick=frac):
      self.updateContents()


  def setLineWidth(self, value):

    if self.nuc.setDisplaySizes(line=value):
      self.updateContents()
    

  def setTubeSize(self, radius):
  
    if self.nuc.setDisplaySizes(tube=radius):
      self.updateContents()
      
  def setRestraintLinesCis(self, isSelected):
        
    if self.nuc.setRestraintsDisplayed(cis=isSelected):
      self.updateContents()
  
  
  def setRestraintLinesTrans(self, isSelected):
    
    if self.nuc.setRestraintsDisplayed(trans=isSelected):
      self.updateContents()
  
  
  def setTextLabels(self, isSelected):
    
    if self.nuc.setTextDisplayed(isSelected):
      self.updateContents()
  
    
  def setChromoLabels(self, isSelected):
  
    if self.nuc.setChromoTextDisplayed(isSelected):
      self.updateContents()    
    
  
  def _setupModelMenu(self, menu):
  
    menu.clear()
    nuc = self.nuc
    
    if nuc:
      nModels = nuc.getNumModels()
      models = range(nModels)
      selected = set(nuc.getDisplayedModels())
      
      for m in models:
        menu.addItem('Model %d' % (m+1), checked=m in selected, 
                     callback=self._setDisplayModels,
                     object=[m,], group=1)        
      
      menu.addItem('All models', 
                   callback=self._setDisplayModels,
                   object=list(models))        
     
      menu.addSeparator()
      value = [m+1 for m in selected]
      rangeEntry = IntRangesEntry(menu, value, minValue=1, maxValue=nModels,
                                  callback=self._setDisplayModelsRange)
      
      rangeEntry.setMaximumWidth(100)                    
                                    
      menu.addItem('Selection:', callback=self._setDisplayModelsRange,widget=rangeEntry)
          
          
  def _setupDispDetailsMenu(self, menu):
    
    menu.clear()
    nuc = self.nuc
    
    if nuc:
      attrs = nuc.display.attrs

      sphDetail = attrs['detailLevel'][0]
      lineSmooth = attrs['detailLevel'][1]
      tubeSmooth = attrs['detailLevel'][2]
      tubeDetail = attrs['detailLevel'][3]
 
      sphereSize = attrs['sizes'][0]
      stickSize = attrs['sizes'][1]
      lineWidth = attrs['sizes'][2]
      tubeSize = attrs['sizes'][3]   
    
      showCis = attrs['options'][2]
      showTrans = attrs['options'][3]
      showChromoLabels = attrs['options'][4]
      showLabels = attrs['options'][5]
      
      sphSizeEntry = FloatSpinBox(menu, sphereSize, minValue=0.001, maxValue=10.0,
                                  step=0.5, callback=self.setSphereSize,
                                  multiplier=1.1)

      stickSizeEntry = FloatSpinBox(menu, stickSize, minValue=0.0, maxValue=1.0,
                                    step=0.25, callback=self.setStickSize)

      sphDetailEntry = IntSpinBox(menu, sphDetail, minValue=0, maxValue=4, step=1,
                                  callback=self.setSphereDetail)
 
      lineWidthEntry = FloatSpinBox(menu, lineWidth, minValue=0.5, maxValue=5.0,
                                    step=0.5, callback=self.setLineWidth)

      lineSmoothEntry = IntSpinBox(menu, lineSmooth, minValue=0, maxValue=4, step=1,
                                   callback=self.setLineSmooth)

      tubeSizeEntry = FloatSpinBox(menu, tubeSize, minValue=0.001, maxValue=10.0,
                                   step=0.5, callback=self.setTubeSize, multiplier=1.1)

      tubeSmoothEntry = IntSpinBox(menu, tubeSmooth, minValue=0, maxValue=4, step=1,
                                   callback=self.setTubeSmooth)

      tubeDetailEntry = IntSpinBox(menu, tubeDetail, minValue=0, maxValue=4, step=1,
                                   callback=self.setTubeDetail)

      menu.addItem('Cis restraints', checked=bool(showCis),
                   callback=lambda a:self.setRestraintLinesCis(a.isChecked()))
      menu.addItem('Trans restraints', checked=bool(showTrans), 
                   callback=lambda a:self.setRestraintLinesTrans(a.isChecked()))
      
                   
      menu.addItem('Chromosome labels', checked=bool(showChromoLabels), 
                   callback=lambda a:self.setChromoLabels(a.isChecked()))
      menu.addItem('Position labels', checked=bool(showLabels), 
                   callback=lambda a:self.setTextLabels(a.isChecked()))
      menu.addSeparator()

      menu.addItem('Ball radius:', widget=sphSizeEntry)
      menu.addItem('Stick width:', widget=stickSizeEntry)
      menu.addItem('Sphere detail:', widget=sphDetailEntry)
      menu.addSeparator()
      
      menu.addItem('Line width:', widget=lineWidthEntry)
      menu.addItem('Line smooth:', widget=lineSmoothEntry)
      menu.addSeparator()
      
      menu.addItem('Tube width:', widget=tubeSizeEntry)
      menu.addItem('Tube smooth:', widget=tubeSmoothEntry)
      menu.addItem('Tube detail:', widget=tubeDetailEntry)
  
  
  def _setTitle(self, title):
    
    self.setWindowTitle('%s: %s' % (PROGRAM, title) )
  
  
  def _getHomeDir(self):
  
    return os.environ.get('HOME') or os.environ.get('HOMEPATH')


  def _resetQtSettings(self):
    
    settings = QtCore.QSettings()
    settings.clear()

    self.dirPathNuc = None
    self.dirPathRef = None
    self.dirPathCont = None
    self.dirPathImage = None
    self.dirPathData = None
    self.recentFiles = [] 
       
    self.restoreState(self.defaultState)
    self.restoreGeometry(self.defaultGeometry)
    self.structurePanel.bgColor = [0.0, 0.0, 0.0, 1.0]
    self.interactomePanel.bgColor = [0.0, 0.0, 0.0, 1.0]

    
  def _readQtSettings(self):
    
    settings = QtCore.QSettings()

    state = settings.value("state", None)
    geometry = settings.value("geometry", None)
    
    dirPathNuc = settings.value("dirPathNuc", None)
    dirPathRef = settings.value("dirPathRef", None)
    dirPathCont = settings.value("dirPathCont", None)
    dirPathImage = settings.value("dirPathImage", None)
    dirPathData = settings.value("dirPathData", None)
    strucBgColor = settings.value("strucBgColor", None)
    recentFiles = settings.value("recentFiles", None)
    mouseMoveFuncs = settings.value("mouseMoveFuncs", None)

    if state:
      self.restoreState(state)
    
    if geometry:
      self.restoreGeometry(geometry)

    if dirPathNuc and pathExists(dirPathNuc) and path.isdir(dirPathNuc):
      self.dirPathNuc = dirPathNuc
      self.fileSystemPanel.openDir(dirPathNuc)
    
    if dirPathRef and pathExists(dirPathRef) and path.isdir(dirPathRef):
      self.dirPathRef = dirPathRef
    
    if dirPathCont and pathExists(dirPathCont) and path.isdir(dirPathCont):
      self.dirPathCont = dirPathCont
    
    if dirPathImage and pathExists(dirPathImage) and path.isdir(dirPathImage):
      self.dirPathImage = dirPathImage

    if dirPathData and pathExists(dirPathData) and path.isdir(dirPathData):
      self.dirPathData = dirPathData
    
    if strucBgColor:
      rgb = [float(x) for x in strucBgColor]
      self.structurePanel.bgColor = rgb
      self.interactomePanel.bgColor = rgb
    
    if recentFiles:
      if not isinstance(recentFiles, list):
        recentFiles = [recentFiles,]
      
      recentFiles = [f for f in recentFiles if pathExists(f) and f.endswith(FILE_EXT)]
        
      self.recentFiles = recentFiles

    if mouseMoveFuncs:
      lbr, mbr, rbr = mouseMoveFuncs
      self._setMouseMotionFunc((int(lbr), 'mouseLeftMoveFunc'))
      self._setMouseMotionFunc((int(mbr), 'mouseMiddleMoveFunc'))
      self._setMouseMotionFunc((int(rbr), 'mouseRightMoveFunc'))
    
    
  def _writeQtSettings(self):
    
    settings = QtCore.QSettings()
    
    settings.setValue("geometry", self.saveGeometry())
    settings.setValue("state", self.saveState())
    
    settings.setValue("dirPathNuc", self.dirPathNuc)
    settings.setValue("dirPathRef", self.dirPathRef)
    settings.setValue("dirPathCont", self.dirPathCont)
    settings.setValue("dirPathImage", self.dirPathImage)
    settings.setValue("dirPathData", self.dirPathData)
    settings.setValue("strucBgColor", self.structurePanel.bgColor)
    
    settings.setValue("recentFiles", self.recentFiles)
  
    panel = self.structurePanel
    lbr = 1 if panel.mouseLeftMoveFunc == panel._mouseRotate else 0
    mbr = 1 if panel.mouseMiddleMoveFunc == panel._mouseRotate else 0
    rbr = 1 if panel.mouseRightMoveFunc == panel._mouseRotate else 0
    mouseMoveFuncs = [lbr, mbr, rbr]
    
    settings.setValue("mouseMoveFuncs", mouseMoveFuncs)
  
  
  def _showSidePanel(self, selected):
  
    panels = [self.overviewPanel, self.roiPanel, self.labelsPanel, self.strucCalcPanel,
              self.fileSystemPanel, self.consolePanel]
    
    if self.leftFrame.isVisible():
      if self.leftFrame.width() < 10:
        sizes = self.splitter.sizes()
        sizes[0] = 250
        sizes[1] -= 250
        self.splitter.setSizes(sizes)
    
    for panel in panels:
      if panel is selected:
        panel.show()
      else:
        panel.hide()

    self.leftFrame.show()
  
  
  def _toggleSidePanel(self, selected):
  
    panels = [self.overviewPanel, self.roiPanel, self.labelsPanel, self.strucCalcPanel,
              self.fileSystemPanel, self.consolePanel]
    
    if self.leftFrame.isVisible():
      if self.leftFrame.width() < 10:
        sizes = self.splitter.sizes()
        sizes[0] = 250
        sizes[1] -= 250
        self.splitter.setSizes(sizes)
       
      elif selected.isVisible():
        for panel in panels:
          panel.hide()
        self.leftFrame.hide()
        return
    
    for panel in panels:
      if panel is selected:
        panel.show()
      else:
        panel.hide()

    self.leftFrame.show()
 
  
  def toggleFileBrowser(self):
    
    self._toggleSidePanel(self.fileSystemPanel)

 
  def toggleStrucCalc(self):
  
    self._toggleSidePanel(self.strucCalcPanel)
 
 
  def  toggleStatsOverview(self):
  
    self._toggleSidePanel(self.overviewPanel)
  
  
  def  toggleLabelsTable(self):
  
    self._toggleSidePanel(self.labelsPanel)


  def  toggleRoiTable(self):
  
    self._toggleSidePanel(self.roiPanel)


  def  toggleConsole(self):
  
    self._toggleSidePanel(self.consolePanel)
  
  
  def showStrucCalc(self):
  
    self._showSidePanel(self.strucCalcPanel)
    
  
  def setRecentFilesMenu(self, menu):
  
    menu.clear()
    
    for filePath in self.recentFiles:
      menu.addItem(filePath, callback=self.openNucFile,
                   object=filePath)
 
  
  def selectTab(self, index):
    
    panels = self.tabbedPanel.widgets
    panels[index].updateContents()
    
    if index in (0,2):
      self.nuc.notifiers.add('/display/attrs/models')
    
    if index < 4:
      self.chromoToolbar.show()
    else:
      self.chromoToolbar.hide()
  
     
  def updateContents(self):
    
    if self.nuc:
      tab = self.tabbedPanel.currentIndex()
      
      self.selectTab(tab)
      self.updateChromosomeToolbar()
      self.updateDataTrackToolbars()
      
      leftPanels = [self.labelsPanel,
                    self.roiPanel,
                    self.overviewPanel,
                    self.strucCalcPanel]
                    
      for panel in leftPanels:
        if panel.isVisible():
          panel.updateContents()

    
  def _checkCopyFile(self, source, detination, actionMsg='copying a file'):
     
    try:
      temp = File(source, 'r')
    
    except Exception as err:
      msg = 'The HDF5 file %s was found to be unreadable and likely corrupt when %s.' % (source, actionMsg)
      msg += ' Further action will be aborted. A backup file may exist as a .temp file.'
      msg += ' Corrupt .temp files may be removed if the .nuc file is OK\n\n'  
      msg += ' Original Python error: %s' % err
      showWarning('Serious error', msg, self)
      return
    
    try:
      copy2(source, detination)
     
    except Exception as err:
      msg = 'A file system error occurred when %s.' % actionMsg
      msg += ' The program will now exit to avoid further problems.'
      msg += ' Check file system write permission and free space.'
      msg += ' Previous data may be present as ".nuc.temp" files.\n\n'
      msg += ' Original Python error: %s' % err
      showWarning('Serious error', msg, self)
      sys.exit(0)
    
    return True
    
    
  def _warning(self, msg):
        
    showWarning('*Warning*', msg, self)
 
 
  def _checkSetNucFile(self, filePath):
    
    
    if self.nuc:
      if self.nuc.root.filename == self._getTempFilePath(filePath):
        return
        
      elif self.nuc.root.filename == filePath:
        return
      
    valid, msg = checkRegularFile(filePath)
    if not valid:
      showWarning('Failed to open file', msg, self)
      return
    
    self.setCursor(QtCore.Qt.WaitCursor)
    tempPath = self._getTempFilePath(filePath)
    otherTemp = self._isFileInUse(filePath)
       
    if otherTemp and checkRegularFile(otherTemp)[0]:
      #lastMod = time.time() -  os.stat(otherTemp).st_atime
      
      msg = 'A temp file "%s" exists. This is could be the result of an unclean shutdown' % otherTemp
      msg += ' or the file could be open in another active session.\n\n Continue to open the file?\n\n'
      
      texts = ['No','Use temp data', 'Yes and remove temp', 'Yes']
      objects = [0,1,2,3]        
      response = showMulti('Query', msg, texts, objects, self)
      
      if response == 0:
        self.unsetCursor()
        return           
      
      elif response == 1:
        if not self._checkCopyFile(otherTemp, tempPath, 'restoring a temp file'):
          self.unsetCursor()
          return     
      
      elif response == 2:
        remove(otherTemp)
        if not self._checkCopyFile(filePath, tempPath, 'creating an active temp file'):
          self.unsetCursor()
          return
      
      else:
        if not self._checkCopyFile(filePath, tempPath, 'creating an active temp file'):
          self.unsetCursor()
          return
         
    else:
      if not self._checkCopyFile(filePath, tempPath, 'creating an active temp file'):
       self.unsetCursor()
       return
            
    msg = 'Save previous %s data?' % PROGRAM
    if self.nuc and self.nucModified() and showYesNo('Confirm', msg, self):
      self.saveNucFile()
      
    elif self.nuc:
      fileName = self.nuc.root.filename
        
      if self._isTempFile(fileName) and pathExists(self._getPersistantFilePath(fileName)): # Discard
        remove(fileName)

      elif fileName == NEW_FILE_TEMP and pathExists(NEW_FILE_TEMP):
        remove(fileName)
    
    self.nuc = Nucleus(tempPath)
    self.nuc._warning = self._warning
    self.nuc.gui = self
    self.nuc.save()
    
    filePath = self._getPersistantFilePath(tempPath)
    self._setTitle(filePath)
    
    if not pathExists(filePath):
      if not self._checkCopyFile(tempPath, filePath, 'saving a new file'):
        self.unsetCursor()
        return
    
    if self.nuc.experimentRef and checkRegularFile(self.nuc.experimentRef)[0]:
      dirName, fileName = path.split(self.nuc.experimentRef)
      self.dirPathRef = dirName

    if self.nuc.genomeRef and checkRegularFile(self.nuc.genomeRef)[0]:
      dirName, fileName = path.split(self.nuc.genomeRef)
      self.dirPathRef = dirName
      
    dirName, fileName = path.split(filePath)
    self.dirPathNuc = dirName
  
    #self.structurePanel.viewCoords = list(self.nuc.display.attrs['view3d'][:3])
    #self.structurePanel.rotation = self.nuc.getGlobalTransform()
    
    self.consolePanel.setContext(self.nuc)
    self.structurePanel.chromoGlLists = {}
    self.structurePanel.restraintGlList = None
    self.structurePanel.glLists = []
    self.coordImage = None
    
    self.updateContents()
    self.unsetCursor()
    
    return True


  def _getFileHash(self, filePath, blockSize=65536):
    
    hashObj = sha256()
    fileObj = open(filePath, 'rb')    
    data = fileObj.read(blockSize)
    
    while data:
      hashObj.update(data)
      data = fileObj.read(blockSize)
    
    return hashObj.digest()
  
  
  def nucModified(self):
    
    if self.nuc:
      return self.nuc.isModified

    else:
      return False
  
  
  def _saveNucFile(self, nuc):
  
    nuc.setStructureView(self.structurePanel.view)
    nuc.setGlobalTransform(self.structurePanel.rotation)
    nuc.save()
    self._writeQtSettings()
  
  
  def saveNucFile(self):
  
    if self.nuc:
      fileName = self.nuc.root.filename
      
      if fileName == NEW_FILE_TEMP:
        self.saveNucFileAs()
      
      elif self._isTempFile(fileName):
        self._saveNucFile(self.nuc)
        destFile = self._getPersistantFilePath(fileName)
        if self._checkCopyFile(fileName, destFile, 'saving data to file'):
          self._setRecentFile(destFile)
        
      else:
        self._saveNucFile(self.nuc)
        self._setRecentFile(fileName)
  
  
  def saveNucFileAs(self):
  
    msg = 'Set %s save file name and location' % PROGRAM
    fileTypes = [FileType(PROGRAM, ['*%s' % FILE_EXT])]
    filePath = selectSaveFile(self, msg, self.dirPathNuc, fileTypes)
  
    if filePath:
      nuc = self.nuc
      root, fex = path.splitext(filePath)
      self.dirPathNuc = path.dirname(filePath)
      
      if fex != FILE_EXT:
        filePath =  root + FILE_EXT
      
      self._saveNucFile(nuc)
      self._setTitle(filePath)

      oldPath = nuc.root.filename
      
      if self._isTempFile(filePath):
        newPath = filePath
      else:
        newPath = self._getTempFilePath(filePath)
      
      if self._isTempFile(oldPath): # Old temp newcomes new temp
        move(oldPath, newPath) 
        
      else: # Old file left, and copied to new temp
        if not self._checkCopyFile(oldPath, newPath, 'saving data to a new file'):
          return
      
      mainPath = self._getPersistantFilePath(newPath)
      if self._checkCopyFile(newPath, mainPath, 'saving data to a new file'): # Make non-temp file
        self._setRecentFile(mainPath) 
           
        nuc.__init__(newPath, nuc.version, nuc.experimentRef, nuc.genomeRef)
         
         
  def newNucFile(self):
    
    msg = 'Save previous %s data?' % PROGRAM
    if self.nuc and self.nucModified() and showYesNo('Confirm', msg, self):
      self.saveNucFile()
    
    if pathExists(NEW_FILE_TEMP): # Unclean shutdown
      remove(NEW_FILE_TEMP)
        
    self.nuc = Nucleus(NEW_FILE_TEMP)
    self.nuc._warning = self._warning
    self.nuc.gui = self
    self.nuc.save()
    
    self._setTitle('*NEW*')
    
    self.consolePanel.setContext(self.nuc)  
    self.structurePanel.chromoGlLists = {}
    self.structurePanel.restraintGlList = None
    self.structurePanel.glLists = []

    self.updateContents()
    #only update at the beginning, not to lose the slider settings
     
     
  def _setRecentFile(self, filePath):
  
    if filePath in self.recentFiles:
      self.recentFiles.remove(filePath)
 
    self.recentFiles.insert(0, filePath)

    if len(self.recentFiles) > 10:
      self.recentFiles = self.recentFiles[:10]
  
  
  def openFiles(self, filePaths, replace=None):

    if not isinstance(filePaths, (list, tuple, set)):
      filePaths = [filePaths,]
    
    for filePath in filePaths:
      dataType, format = getFileType(filePath)
       
      if dataType == 'structure':
        self.importCoords(filePath, format)
 
      elif dataType == 'dataTrack':
        self.importDataTrack(filePath, format)
 
      elif dataType == 'contacts':
        self.importContacts(filePath, format)

      elif dataType == 'interactions':
        self.importInteractions(filePath, format)
 
      elif dataType == 'nuc':
        self.openNucFile(filePath, replace)
 
      else:
        self._warning('Type of file "%s" not understood' % filePath)
        break

      
  def openNucFile(self, filePath=None, replace=None):
  
    if not filePath:
      msg = 'Select %s file' % FILE_EXT
      fileTypes = [FileType(PROGRAM, ['*%s' % FILE_EXT])]
      filePath = selectFile(self, msg, fileTypes=fileTypes,
                            directory=self.dirPathNuc)
    
    if filePath:
      valid, msg = checkRegularFile(filePath)
      
      if valid:
        isNew = self._checkSetNucFile(path.abspath(filePath))
        
        if isNew:
          self.fileSystemPanel.openDir(self.dirPathNuc)
        
        self._setRecentFile(filePath)
        
      else:
        showWarning('Failed to open file', msg, self)


  def setExperimentRef(self, filePath=None):
    
    if not self.nuc:
      return
    
    if not filePath:
      msg = 'Select %s file for experiment reference' % FILE_EXT
      fileTypes = [FileType(PROGRAM, ['*%s' % FILE_EXT])]
                                
      dialog = FileDialog(self, msg, self.dirPathRef, doSave=True,
                          default='tissue_cell.nuc', fileTypes=fileTypes,
                          warnOverwrite=False)
               
      filePath = dialog.getFile()
   
    if filePath:
      valid, msg = checkRegularFile(filePath)
      self.dirPathRef = path.dirname(filePath)
      
      if valid:
        if self.nuc.setExperimentRef(filePath):
          self.updateContents()
        
      else:
        try:
          rNuc = Nucleus(filePath, experimentRef=False, genomeRef=False, mode='a')
          rNuc.save()
          
          if self.nuc.setExperimentRef(filePath): 
            self.updateContents()
 
        except Exception() as err:
          msg = 'Referece .nuc file "%s" could not be created.\nPython error:\n%s'
          self._warning(msg % (filePath, err))
      
      
  def setGenomeRef(self, filePath=None):
    
    if not self.nuc:
      return
    
    if not filePath:
      msg = 'Select %s file for genome reference' % FILE_EXT
      fileTypes = [FileType(PROGRAM, ['*%s' % FILE_EXT])]
      
      dialog = FileDialog(self, msg, self.dirPathRef, doSave=True,
                          default='genome_build.nuc', fileTypes=fileTypes,
                          warnOverwrite=False)
               
      filePath = dialog.getFile()
       
    if filePath:
      valid, msg = checkRegularFile(filePath)
      self.dirPathRef = path.dirname(filePath)
      
      if valid:
        if self.nuc.setGenomeRef(filePath):
          self.updateContents()
        
      else:
        try:
          rNuc = Nucleus(filePath, experimentRef=False, genomeRef=False, mode='a')
          rNuc.save()
          
          if self.nuc.setGenomeRef(filePath): 
            self.updateContents()
 
        except Exception() as err:
          msg = 'Referece .nuc file "%s" could not be created.\nPython error:\n%s'
          self._warning(msg % (filePath, err))
  
  
  def getIcon(self, fileName):
    
    if fileName:
      filePath = path.join(ICON_DIR, fileName)
      
    else:
      filePath = None
    
    return QtGui.QIcon(filePath) 
     
       
  def _notImplemented(self, *args):
  
    showWarning('Cannot proceed', 'Function not implemented', self)
  

  def closeEvent(self, event=None): # Overwrite
    
    if self.nuc:
      if self.nucModified():
        msg = 'Save current data to file before exit?'
        fileName = self.nuc.root.filename
 
        answer = showSaveDiscardCancel('Query', msg, self)
 
        if answer is None: # Cancelled
          if event:
            event.ignore()
 
          return
 
        elif answer: # Yes, save
          self._saveNucFile(self.nuc)
          self.nuc.root.close()

          if self._isTempFile(fileName):
            destFile = self._getPersistantFilePath(fileName)
            if self._checkCopyFile(fileName, destFile, 'saving data to file'):
              remove(fileName)
 
        elif self._isTempFile(fileName) and pathExists(self._getPersistantFilePath(fileName)): # No, discard
          if pathExists(fileName):
            remove(fileName)
          else:
            print('Warning temprary file "%s" deleted prematurely' % fileName)  

        elif fileName == NEW_FILE_TEMP and pathExists(NEW_FILE_TEMP):
          remove(fileName)
      
      else:
        fileName = self.nuc.root.filename
        
        if self._isTempFile(fileName) and pathExists(self._getPersistantFilePath(fileName)): # Discard
          remove(fileName)
        
        elif fileName == NEW_FILE_TEMP and pathExists(NEW_FILE_TEMP):
          remove(fileName)
        
      self._writeQtSettings()
      
      if event:
        event.accept()
 
    elif event:
      self._writeQtSettings()
      event.accept()
      
      
  def removeUnrestrainedPoints(self):
    
    
    restraintDict = self.nuc.getRestraints(usePositions=True)
    posDict = {}
       
    for chromoA, chromoB in restraintDict:
      restraints = restraintDict[(chromoA, chromoB)]
      posA = restraints[:,0]
      posB = restraints[:,1]
      
      if chromoA in posDict:
        posDict[chromoA].update(posA)
      else:
        posDict[chromoA] = set(posA)
      
      if chromoB in posDict:
        posDict[chromoB].update(posB)
      else:
        posDict[chromoB] = set(posB)
        
    for chromo in posDict:
      posDict[chromo] = sorted(posDict[chromo])
      
    # Restraint indices need to be reset...
      
    self.nuc.addChromosomes(posDict, interpolateCoords=True)
    self.updateContents()
    
 
  def setCoordsRandomWalk(self):
    
    chromos = self.nuc.getDisplayedChromosomes()
    
    if not chromos:
      return
    
    chromoStr = ', '.join(chromos)
    msg = 'Really randomise structures for chromosomes: %s ?' % (chromoStr,)
    
    if not showOkCancel('Confirm', msg, parent=self):
      return

    models = list(range(self.nuc.getNumModels()))
    self.nuc.setRandomCoords(models, chromos, randWalk=True, 
                             randSeed=int(time.time()))       
    
    self.updateContents()


  def setCoordsLinearStack(self):
    
    chromos = self.nuc.getDisplayedChromosomes()
    
    if not chromos:
      return
    
    chromoStr = ', '.join(chromos)
    msg = 'Really set structures for chromosomes: %s ?' % (chromoStr,)
    
    if not showOkCancel('Confirm', msg, parent=self):
      return
    
    sortChromos = []
    for chromo in self.nuc.getDisplayedChromosomes():
      start, end = self.nuc.getChromosomeLimits(chromo)
      delta = end - start
      if not delta:
        continue
      
      sortChromos.append(chromo)
    
    nChromos = len(sortChromos) 
    for i, chromo in enumerate(sortChromos):
      particGroup = self.nuc._getParticleGroup() # Will use current structure
      seqPos = array(particGroup[chromo]['positions'])
      coords, centre = getLinearChromoCoords(100e6, seqPos, i, nChromos, 5.0)
      for model in range(self.nuc.getNumModels()):
        self.nuc.setModelCoords(coords, model, chromosomes=[chromo,])
    
    self.updateContents()
        
      
  def setCoordsGreatCircle(self):
    
    chromos = self.nuc.getDisplayedChromosomes()
    
    if not chromos:
      return
    
    chromoStr = ', '.join(chromos)
    msg = 'Really set structures for chromosomes: %s ?' % (chromoStr,)
    
    if not showOkCancel('Confirm', msg, parent=self):
      return

    offsets = {}
    totalSeqLen = 0
    sortChromos = []
    
    for chromo in self.nuc.getDisplayedChromosomes():
      start, end = self.nuc.getChromosomeLimits(chromo)
      delta = end - start
      if not delta:
        continue
      
      delta  *= 1.05
      sortChromos.append(chromo)
      offsets[chromo] = totalSeqLen
      totalSeqLen += delta   
    
    for i, chromo in enumerate(sortChromos):
      particGroup = self.nuc._getParticleGroup() # Will use current structure
      seqPos = array(particGroup[chromo]['positions'])
      coords, centre = getArcChromoCoords(seqPos, offsets[chromo], totalSeqLen, 10.0)
    
      for model in range(self.nuc.getNumModels()):
        self.nuc.setModelCoords(coords, model, chromosomes=[chromo,])

    self.updateContents()

             
  def importExperimentRefData(self):
  
    if self.nuc:
      rNuc = self.nuc.getExperimentRefNuc()
      if not rNuc:
        msg = 'Experiment reference file not set. Set now?'
        if showOkCancel('Confirm', msg, parent=self):
          self.setExperimentRef()
        else:
          return
     
      rNuc = self.nuc.getExperimentRefNuc()
      if not rNuc:
        return
          
      msg = 'Select experimental data to load'
      fileTypes = []
      filePath = selectFile(self, msg, self.dirPathData, fileTypes)
 
      if filePath:
        self.dirPathData = path.dirname(filePath)
        fileName, fex = path.splitext(path.split(filePath)[1])
        code = askString('Query', 'Enter data name:', fileName, parent=self)
        
        if code:
          rNuc.importDataTrack(filePath, EXTERNAL, code, format=None)
          rNuc.save()
          self.nuc.updateProxyDataTracks(rNuc, EXTERNAL)
          self.updateContents()
        
          
  def importGenomeRefData(self):
  
    if self.nuc:
      rNuc = self.nuc.getGenomeRefNuc()
      if not rNuc:
        msg = 'Genome reference file not set. Set now?'
        if showOkCancel('Confirm', msg, parent=self):
          self.setGenomeRef()
        else:
          return
     
      rNuc = self.nuc.getGenomeRefNuc()
      if not rNuc:
        return
          
      msg = 'Select genome data to load'
      fileTypes = []
      filePath = selectFile(self, msg, self.dirPathData, fileTypes)
 
      if filePath:
        self.dirPathData = path.dirname(filePath)
        fileName, fex = path.splitext(path.split(filePath)[1])
        code = askString('Query', 'Enter data name:', fileName, parent=self)
        
        if code:
          rNuc.importDataTrack(filePath, INNATE, code, format=None)
          rNuc.save()
          self.nuc.updateProxyDataTracks(rNuc, INNATE)
          self.updateContents()
  
  
  def _getUserInputWord(self, msg, default=''):
  
    word = askString('Query', msg, default, parent=self) or ''
    word = word.strip()
    
    while re.search('\W', word):
      word = askString('Query', msg+'\n(letters, numbers and underscore only)', default, parent=self) or ''
      word = word.strip()
      
      if not word:
        return
    
    return word
        
  
  def importGenomeFeature(self, filePath=None):
  
    if not filePath:
      fileTypes = [FileType('GFF/GTF', ['*.gff*','*.gtf*'])]
      msg = 'Select geneome feature file to load'
      filePath = selectFile(self, msg, self.dirPathData, fileTypes)
    
    if filePath:
      self.dirPathData = path.dirname(filePath)
      
      from formats.GFF import getFeatures
      from gui.qtgui.MessageDialog import MessageDialog
      
      self.setCursor(QtCore.Qt.WaitCursor)
      featureDict = getFeatures(filePath)
      self.unsetCursor()
      
      if not featureDict:
        return
      
      features = sorted(featureDict.keys())
      feature = features[0]
      
      if len(features) > 1:
        msg = 'What type of feature should be imported?'
        dialog = MessageDialog('Query', 'Select', msg,
                               QtGui.QMessageBox.Question, self)

        texts = ['{} ({:,})'.format(f, featureDict[f]) for f in features]
        pulldown = PulldownList(dialog, texts, features, index=0)
        dialog.layout().addWidget(pulldown, 2, 2)
                
        dialog.setStandardButtons(QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel)
        dialog.setDefaultButton(QtGui.QMessageBox.Ok)
 
        resp = dialog.exec_()
        
        if resp != QtGui.QMessageBox.Ok:
          return
          
        feature = pulldown.get()
      
      fileName = splitExtension(path.split(filePath)[1])[0]
      code = self._getUserInputWord('Enter data name:', feature)
      
      if code:
        self.setCursor(QtCore.Qt.WaitCursor)
 
        counts = self.nuc.importDataTrack(filePath, 'GFF', INNATE, code, None, feature=feature)
 
        self.updateContents()
        self.unsetCursor()
 
        if counts:
          nChromos, nRegions = counts
          msg = 'Imported %d features for %d chromosomes'
          showInfo('Info', msg % (nRegions, nChromos))
      
                            
  def importDataTrack(self, filePath=None, format=None):
 
    if not filePath:
    
      allTypes = []
      for fmt in DATA_TRACK_FORMATS:
        allTypes += DATA_TRACK_FORMATS[fmt]
    
      fileTypes = [FileType('Supported formats', allTypes)]
      for fmt, exts in DATA_TRACK_FORMATS.iteritems():
        fileTypes.append( FileType(fmt, exts) )
 
      msg = 'Select data track files to load'
      filePaths = selectFiles(self, msg, self.dirPathData, fileTypes)
    
    for filePath in filePaths:
      self.dirPathData = path.dirname(filePath)
      fileRoot, fileExt = splitExtension(filePath)
      fileExt = '*' + fileExt
      
      if not format:
        for format in DATA_TRACK_FORMATS:
          if fileExt in DATA_TRACK_FORMATS[format]:
            break
 
        else:
          self._warning('File extension "%s" not supported' % fileExt[1:])
          return
      
      if format == 'Nuc3D':
        self.nuc.importDataTrack(filePath, format, None, None, None)
        self.updateContents()
        self.unsetCursor()
        
      else:
       
        fileName = path.split(fileRoot)[1]
        code = self._getUserInputWord('Enter data name:', fileName)
 
        if code:
          binSize = askInteger('Query', 'Data bin size (kb):', 0, parent=self) or 0
          binSize *= 1000
 
          self.setCursor(QtCore.Qt.WaitCursor)
 
          msg = 'What category is the imported data?'
          texts = ['External experiment',
                   'Innate to genome',
                   'Contact/structure derived']
          objects = [EXTERNAL, INNATE, DERIVED]
          source = showMulti('Query', msg, texts, objects, self)
 
          #try:
          counts = self.nuc.importDataTrack(filePath, format, source, code, binSize or None)

          self.updateContents()
          self.unsetCursor()
 
          if counts:
            nChromos, nRegions = counts
            msg = 'Imported %d values for %d chromosomes'
            showInfo('Info', msg % (nRegions, nChromos))

        
        
  def importInteractions(self, filePath=None, format=None):
 
    if not filePath:
    
      allTypes = []
      for fmt in INTERACTIONS_FORMATS:
        allTypes += INTERACTIONS_FORMATS[fmt]
    
      fileTypes = [FileType('Supported formats', allTypes)]
      for fmt, exts in INTERACTIONS_FORMATS.iteritems():
        fileTypes.append( FileType(fmt, exts) )
 
      msg = 'Select interactions file to load'
      filePath = selectFile(self, msg, self.dirPathData, fileTypes)
    
    if filePath:
      self.dirPathData = path.dirname(filePath)
      fileRoot, fileExt = splitExtension(filePath)
      fileExt = '*' + fileExt
      
      if not format:
        for format in INTERACTIONS_FORMATS:
          if fileExt in INTERACTIONS_FORMATS[format]:
            break
 
        else:
          self._warning('File extension "%s" not supported' % fileExt[1:])
          return
          
          
      if format == 'Nuc3D':
        self.nuc.importInteractions(filePath, format, None)
        self.updateContents()
        self.unsetCursor()
      
      else:
        fileName = path.split(fileRoot)[1]
        code = self._getUserInputWord('Enter data name:', fileName)
 
        if code:
 
          self.setCursor(QtCore.Qt.WaitCursor)
 
          counts = self.nuc.importInteractions(filePath, format, code)

          self.updateContents()
          self.unsetCursor()
 
          if counts:
            nChromos, nRegions = counts
            msg = 'Imported %d interactions for %d chromosomes'
            showInfo('Info', msg % (nRegions, nChromos))
 
  
  def _convertTextToRegion(self, text):
    
    def err(msg=None):
      if not msg:
        msg = 'Could not interpret chromosome region "%s"\n' % text
        msg += 'Region specifications should be <chr>:<middle> or <chr>:(start>-<end>'
      showWarning('Error', msg, self)    
    
    chrPos = text.split(':')
    
    if len(chrPos) != 2:
      err()
      return
    
    chromo, pos = chrPos
    chromo = chromo.strip().upper()
    pos = pos.strip()
    pos = pos.replace(' ', '')
    pos = pos.replace(',', '')
    
    if chromo not in self.nuc.getDisplayedChromosomes():
      err('Chromosome "%s" not currently displayed' % chromo)
    
    pos = pos.split('-')
    
    if len(pos) > 2:
      err()
      return
      
    for i, p in enumerate(pos):
      
      try:
        p = int(p)
        pos[i] = p
      
      except ValueError:
        err('Postion "%s" could not be interpreted as a number' % p)
        return
    
    if len(pos) == 2:
      pos = sorted(pos)  
      width = pos[1] - pos[0]
      pos = sum(pos)/2
       
    else:
      width = None
      pos = pos[0]
    
    return (chromo, pos, width)
  
  
  def _isFileContactMatrix(self, filePath, format):
  
    if format in ('SAM', 'PFE'):
      return False
    
    elif format == 'TVS':
      return True
    
    elif format == 'HDF5':
      from h5py import File
      root = File(filePath, mode='r')
      isMatrix = 'contactMatrix' in root
      root.close()
      return isMatrix
   
    elif format == 'JSON':
      fileObj = open(filePath, 'r')
      data = fileObj.read(128)
      isMatrix = 'contactMatrix' in data
      fileObj.close()
      return isMatrix
    
    elif format == 'NDArray':
     from numpy import load
     dataDict = load(filePath)
     for key in dataDict:
       if key.startswith('contactMatrix/'):
         isMatrix = True
         break
     
     else:
       isMatrix = False
     return isMatrix

  
  def importImageCoords(self, filePath=None):
  
    from numpy import loadtxt, load
  
    if not filePath:
    
      #allTypes = []
      #for exts in CONTACT_FORMATS.values():
      #  allTypes += exts
    
      #fileTypes = [FileType('SAM/BAM', CONTACT_FORMATS['SAM']),
      #             FileType('Paired Fragment End', CONTACT_FORMATS['PFE']),
      #             FileType('Supported formats', allTypes)]
      
      #for format, exts in CONTACT_FORMATS.iteritems():
      #  if format in ('SAM', 'PFE'):
      #    continue
          
      #  fileTypes.append( FileType(format, exts) )
     
      fileTypes = [FileType('All files', ['*.txt','*.npy','*.csv']),]
      msg = 'Select image point coords data file to load'
      filePath = selectFile(self, msg, self.dirPathImage, fileTypes)
  
    if filePath:
      self.dirPathImage, fileName = path.split(filePath)
      defaultName = path.splitext(fileName)[0]
      code = self._getUserInputWord('Image name:', defaultName)
      
      if not code:
        return
      
      if filePath.endswith('.npy'):
        data = load(filePath)
      
      else:
        with open(filePath) as file_obj:
          data = []
          
          for line in file_obj:
            row = line[:-1].split(',')
            row = [float(x) for x in row]
            data.append(row)
          
          data = array(data)
          
          while data.max() < 10000:
            data *= 10.0
          
          print data
          
      self.nuc.setImageCoords(code, data)
  
  
  def centreImageCoords(self):
  
    self.nuc.centreImageCoords()
    
  
  def importContacts(self, filePath=None, format='SAM'):

    if filePath:
      filePaths = [filePath]
    
    else:
      allTypes = []
      for exts in CONTACT_FORMATS.values():
        allTypes += exts
    
      fileTypes = [FileType('Nuc Chromo Contact', CONTACT_FORMATS['NCC']),
                   FileType('SAM/BAM', CONTACT_FORMATS['SAM']),
                   FileType('Paired Fragment End', CONTACT_FORMATS['PFE']),
                   FileType('Supported formats', allTypes)]
      
      for fmt, exts in CONTACT_FORMATS.iteritems():
        if fmt in ('SAM', 'PFE'):
          continue
          
        fileTypes.append( FileType(fmt, exts) )
 
      msg = 'Select contact data file to load'
      filePaths = selectFiles(self, msg, self.dirPathCont, fileTypes) or []
     
    
    for filePath in filePaths:
    
      self.dirPathCont = path.dirname(filePath)
      fileRoot, fileExt = splitExtension(filePath)
      fileName = path.basename(fileRoot)
      defaultName = re.sub( r'\W+', '_', fileName)
      defaultName = re.sub(r'_+', '_', defaultName)
      
      fileExt = '*' + fileExt
 
      for format in CONTACT_FORMATS:
        if fileExt in CONTACT_FORMATS[format]:
          break
 
      else:
        self._warning('File extension "%s" not supported' % fileExt[1:])
        return
      
      self.setCursor(QtCore.Qt.WaitCursor)
      
      if format == 'Nuc3D':
        self.nuc.importContacts(filePath, format, None, None,
                                isSingleCell=None, updateChromos=True)

        self.updateContents()
        self.unsetCursor()
        
      else:
          
        isMatrix = self._isFileContactMatrix(filePath, format)
 
        if isMatrix:
          isSingleCell = False
          binSize = None

        else:
          isSingleCell = showYesNo('Query', 'Is the data single-cell?', parent=self)
 
          if isSingleCell:
            binSize = None
 
          else:
            binSize = askInteger('Query', 'Bin size (kb)',
                                 initialValue=50, minValue=1,
                                 maxValue=1000000, parent=self) or 50
            binSize *= 1000
 
 
        groupName = self._getUserInputWord('Contact group name:', defaultName)
 
        if not groupName:
          return
 
        msg = 'Contact dataset "%s" already exists. Overwite existing data?'
        if self.nuc.getContactGroup(groupName) and not showYesNo('Confirm', msg % groupName, parent=self):
          return
      
        self.nuc.importContacts(filePath, format, groupName, binSize,
                                isSingleCell=isSingleCell, updateChromos=True)
      
      if filePaths:
        self.updateContents()
        self.unsetCursor()
                             
  
  def mirrorStructure(self):
    
    nuc = self.nuc
    if nuc:
      chromos = nuc.getDisplayedChromosomes()
      nuc.modelMirror(models=None, chromosomes=chromos)
      self.updateContents()


  def centerStructure(self):
    
    nuc = self.nuc
    
    if nuc:
      chromos = nuc.getDisplayedChromosomes()
      nuc.modelCentre(models=None, chromosomes=chromos)
      self.updateContents()
      
      
  def removeIsolatedContacts(self):
     
    nuc = self.nuc
    
    if nuc:
      threshold = askInteger('Query', 'Input threshold separation',
                             initialValue=int(1e6), minValue=(1e3),
                             maxValue=int(1e8), parent=self)
                             
      if threshold:
        contactGroups = nuc.getSelectedContactGroups()
        
        n = 0
        for group in contactGroups:
          n += nuc.removeIsolatedContacts(group, threshold=threshold)
            
        self.updateContents()
        showInfo('Info', 'Removed %d isolated contacts' % (n))
      
      
  def removeViolatedContacts(self):
  
    if self.nuc:
      threshold = askFloat('Query', 'Distance threshold',
                            initialValue=5.0, minValue=1e-3,
                            maxValue=1e4, parent=self)

      if threshold:
        chromos = self.nuc.getDisplayedChromosomes()
        contactGroups = self.nuc.getSelectedContactGroups()
        n = 0
        for group in contactGroups:
          n += self.nuc.removeViolatedContacts(group, chromos, threshold)
        
        self.updateContents()
        showInfo('Info', 'Removed %d violated contacts' % (n))


  def resolveAmbigousContacts(self):
  
    if self.nuc:
      ideal_dist = askFloat('Query', 'Ideal distance',
                            initialValue=2.0, minValue=1e-3,
                            maxValue=1e4, parent=self)

      if ideal_dist:
        chromos = self.nuc.getDisplayedChromosomes()
        contactGroups = self.nuc.getSelectedContactGroups()
        n = 0
        for group in contactGroups:
          n += self.nuc.resolveAmbiguousContacts(group, chromos, ideal_dist)
        
        self.updateContents()
        showInfo('Info', 'Resolved %d ambigous contacts' % (n))


  def extractViolatedContacts(self):
  
    if self.nuc:
      msg = 'Contact group name? (letters, numbers and underscore only)'
      groupName = askString('Query', msg, 'violated', self)
      
      if not groupName:
        return 
      
      threshold = askFloat('Query', 'Threshold (relative to mean distance)',
                            initialValue=4.0, minValue=1e-3,
                            maxValue=1e4, parent=self)

      if threshold:
        chromos = self.nuc.getDisplayedChromosomes()
        contactGroups = self.nuc.getSelectedContactGroups()
        
        n = 0
        for source in contactGroups:
          n += self.nuc.getViolatedContacts(source, groupName, chromos, threshold)
        
        self.updateContents()
        showInfo('Info', 'Identified %d violated contacts' % (n))


  def mergeContacts(self):
  
    if self.nuc:
      contactGroups = self.nuc.getSelectedContactGroups()
      
      if len(contactGroups) < 2:
        msg = 'Must have at least two contact groups selected/displayed'
        showWarning('Warning', msg, self)
        return
               
      msg = 'Merged contact group name? (letters, numbers and underscore only)'

      group_a = contactGroups[0]
      group_merge = askString('Query', msg, group_a+'_merge', self)
      
      if not group_merge:
        return 
      
       
      for i, group_b in enumerate(contactGroups[1:]):
        if i == len(contactGroups)-2:
          name = group_merge
        else:
          name = 'merge_temp_%d' % i
        
        self.nuc.mergeContacts(group_a, group_b, name=name, remove_frac=None)
        group_a = name # Merge third group etc into merge of first two
        
      self.updateContents()  
        
  
  def calcObsVsExpContacts(self):
  
    contactGroups = self.nuc.getSelectedContactGroups()
    
    if contactGroups:
      for group in contactGroups:
        self.nuc.makeObsVsExpContactMap(group)
 
      self.updateContents()

 
  def generateRandomContacts(self):
 
    if self.nuc:
      neighbours = askInteger('Query', 'Number of closest neighbours to consider:',
                              initialValue=10, minValue=1,
                              maxValue=100, parent=self)
      if not neighbours:
        return
     
      chromos = self.nuc.getDisplayedChromosomes()
      nRestraints = askInteger('Query', 'Number of contacts to generate:',
                               initialValue=len(chromos)*5000, minValue=1,
                               maxValue=int(1e8), parent=self)

      if not nRestraints:
        return
      
      groupName = ' '
      while re.search('\W', groupName):
        msg = 'Contact group name? (letters, numbers and underscore only)'
        groupName = askString('Query', msg, 'structural', self)
        
        if not groupName:
          return 
      
      models = self.nuc.getDisplayedModels() or [0,]
      n = 0
      
      for i, model in enumerate(models):
        replace = i == 0 
        n += self.nuc.setRandomContacts(groupName, model, nRestraints, chromos,
                                        neighbours, replace)
      
      self.updateContents()
      showInfo('Info', 'Generated %d contacts' % (n))


  def calcContactCorrelations(self):
  
    nuc = self.nuc

    contactGroups = nuc.getSelectedContactGroups()
    if contactGroups:
      groupName = contactGroups[0]
    
    else:
      msg = 'No contact groups currently displayed'
      nuc._warning(msg)
      return # Nothing selected  
      
    chromosomes = nuc.getDisplayedChromosomes()
    
    nuc.plotContactCovariance(groupName, chromosomes)
    

  def normalisePopContacts(self):
  
    nuc = self.nuc

    contactGroups = nuc.getSelectedContactGroups()
    if contactGroups:
      groupName = contactGroups[0]
    
    else:
      msg = 'No contact groups currently displayed'
      nuc._warning(msg)
      return # Nothing selected  

    nuc.normaliseContacts(groupName)
  
  def _setupExportContactsMenu(self, menu):
  
    menu.clear()
    nuc = self.nuc
    
    for name, group in nuc.getContactGroups():
      menu.addItem(name, object=name,
                   callback=self.exportContacts)
                   
    menu.addSeparator()
     
    checked = bool(self.exportMatrixAction and self.exportMatrixAction.isChecked())
    
    menu.addItem('&Sparse pairs', object=0, checked=not checked, group=0)
    action = menu.addItem('&Full matrix', object=1, checked=checked, group=0)
    self.exportMatrixAction = action
      

  def _setupExportStructureMenu(self, menu):

    menu.clear()
    nuc = self.nuc

    codes, names = nuc.getStructureNames()
    
    for code, name in zip(codes, names):
      menu.addItem('%s:%s' % (code, name), object=code,
                   callback=self.exportCoords)
      
    menu.addSeparator()
     
    checked = bool(self.exportCurrentModelsAction and self.exportCurrentModelsAction.isChecked())
    
    action = menu.addItem('&Displayed model(s)', object=0, checked=checked, group=0)
    self.exportCurrentModelsAction = action
    menu.addItem('&All models', object=1, checked=not checked, group=0)


  def _setupExportDataMenu(self, menu):

    menu.clear()
    nuc = self.nuc
    
    codes = []
    counts = {}
 
    for typ in nuc.dataTracks:
      for code in nuc.dataTracks[typ]:
        if 'options' in nuc.dataTracks[typ][code].attrs:
          codes.append((typ, code))
          if code in counts:
            counts[code] += 1
          else:
            counts[code] = 1
 
    codes.sort()
    
    for i, (typ, code) in enumerate(codes):
      if counts[code] > 1:
        name = '%s (%s)' % (code, typ)
      else:
        name = code
    
      if (i>0) and (typ != codes[i-1][0]):
        menu.addSeparator()
    
      color = nuc.getDataTrackColor(typ, code).rgbHex
      icon = Icon(color=color)
      menu.addItem(name, icon=icon, object=(typ, code),
                   callback=self.exportDataTrack)
      
  def _setupCreateDensityDataMenu(self, menu):

    menu.clear()
    nuc = self.nuc
    
    codes = []
    counts = {}
 
    for typ in nuc.dataTracks:
      for code in nuc.dataTracks[typ]:
        if 'options' in nuc.dataTracks[typ][code].attrs:
          codes.append((typ, code))
          if code in counts:
            counts[code] += 1
          else:
            counts[code] = 1
 
    codes.sort()
    
    for i, (typ, code) in enumerate(codes):
      if counts[code] > 1:
        name = '%s (%s)' % (code, typ)
      else:
        name = code
    
      if (i>0) and (typ != codes[i-1][0]):
        menu.addSeparator()
    
      color = nuc.getDataTrackColor(typ, code).rgbHex
      icon = Icon(color=color)
      menu.addItem(name, icon=icon, object=(typ, code),
                   callback=self.createSpatialDensityDataTrack)


  def _setupExportInteractionsMenu(self, menu):

    menu.clear()
    nuc = self.nuc
    codes = sorted(nuc.interactions.keys())

    for code in codes:
      color = nuc.getInteractionsColor(code).rgbHex
      icon = Icon(color=color)
      menu.addItem(code, icon=icon, object=code,
                   callback=self.exportInteractions)
 
 
  def importCoords(self, filePath=None, format=None):

    if not filePath:
    
      allTypes = []
      for fmt in ('TSV', 'JSON', 'HDF5', 'NDArray', 'Nuc3D'): # Not PDB or Excel
        allTypes += STRUCTURE_FORMATS[fmt]
    
      fileTypes = [FileType('Supported formats', allTypes)]
      for fmt, exts in STRUCTURE_FORMATS.iteritems():
        fileTypes.append( FileType(fmt, exts) )
 
      msg = 'Select coordinate data file to load'
      filePath = selectFile(self, msg, self.dirPathData, fileTypes)
    
    if filePath:
      fileRoot, fileExt = splitExtension(filePath)
      fileExt = '*' + fileExt
      
      if not format:
        for format in STRUCTURE_FORMATS:
          if fileExt in STRUCTURE_FORMATS[format]:
            break
 
        else:
          self._warning('File extension "%s" not supported' % fileExt[1:])
          return
    
      self.dirPathData = path.dirname(filePath)
      self.setCursor(QtCore.Qt.WaitCursor)
      
      if format == 'Nuc3D':
        self.nuc.importCoords(filePath, format, None)
        self.updateContents()
        self.unsetCursor()
        
      else:
        structure = self.nuc.getNextStructureCode()
        counts = self.nuc.importCoords(filePath, format, structure)
        nChromos, nModels, nCoords = counts

        self.nuc.setCurrentStructure(structure)
        self.updateContents()
        self.unsetCursor()
 
        msg = 'Imported %d coordinates for %d model(s) of %d chromosomes'
        showInfo('Info', msg % (nCoords, nModels, nChromos))

     
  def exportContacts(self, group):
  
    msg = 'Select contacts export filename'
    
    options = (('Paired Fragment End', 'PFE'),
               ('HDF5 binary', 'HDF5'),
               ('NumPy binary','NDArray'),
               ('Nuc JSON', 'JSON'),
               ('Text matrix','TSV'))
      
    isMatrix = self.exportMatrixAction.isChecked()
    
    fileTypes = []
    for label, key in options:
      if isMatrix:
        if (key == 'PFE'):
          continue
     
      else:
        if (key == 'TSV'):
          continue
        
      fileTypes.append( FileType(label, CONTACT_FORMATS[key]) )
    
    fileName = 'Contacts_%s_Export' % (group, )
    
    dialog = FileDialog(self, msg, directory=self.dirPathData,
                        doSave=True, default=fileName,
                        fileTypes=fileTypes)
  
    filePath = dialog.getFile()
    
    if filePath:
      fileRoot, fileExt = splitExtension(filePath)
      filePat = '*' + fileExt.lower()
      
      for format in CONTACT_FORMATS:
        if filePat in CONTACT_FORMATS[format]:
          break
      
      else: # Rename
        selected = dialog.selectedNameFilter()
        
        for label, key in options:
          if selected.startswith(label):
            format = key
            fileExt = CONTACT_FORMATS[format][0][1:]
            filePath = fileRoot + fileExt
            break
      
      self.dirPathData = path.dirname(filePath)
         
      if pathExists(filePath):
        if not showOkCancel('Confirm', 'Overwrite "%s"?' % fileName, parent=self):
          return
      
      self.setCursor(QtCore.Qt.WaitCursor)
      
      if isMatrix:
        binSize = self.nuc.getContactsBinSize(group) or 1
        
        if binSize < 2:
          binSize = askInteger('Query', 'Bin size (kb)',
                               initialValue=50, minValue=1,
                               maxValue=1000000, parent=self) or 50
          binSize *= 1000
        
      else:
        binSize = None
      
      counts = self.nuc.exportContacts(filePath, format, group, binSize)
      
      self.unsetCursor()  
      
      if counts:
        showInfo('Info', 'Exported %d contacts' % counts)

  
  def exportCoords(self, structure):
  
    msg = 'Select data export file name'
    
    options = (('Nuc Text', 'TSV'),
               ('Nuc JSON', 'JSON'),
               ('HDF5 binary', 'HDF5'),
               ('NumPy binary','NDArray'),
               ('Protein Data Bank', 'PDB'),
               ('Excel spreadsheet', 'XLSX'))
    
    fileTypes = []
    for label, key in options:
      fileTypes.append( FileType(label, STRUCTURE_FORMATS[key]) )
    
    fileName = 'Structure_%s_Export' % (structure, )
    
    dialog = FileDialog(self, msg, directory=self.dirPathData,
                        doSave=True, default=fileName,
                        fileTypes=fileTypes)
  
    filePath = dialog.getFile()
    
    if filePath:
      fileRoot, fileExt = splitExtension(filePath)
      filePat = '*' + fileExt.lower()
      
      for format in STRUCTURE_FORMATS:
        if filePat in STRUCTURE_FORMATS[format]:
          break
      
      else: # Rename
        selected = dialog.selectedNameFilter()
        
        for label, key in options:
          if selected.startswith(label):
            format = key
            fileExt = STRUCTURE_FORMATS[format][0][1:]
            filePath = fileRoot + fileExt
            break
      
      self.dirPathData = path.dirname(filePath)
         
      if pathExists(filePath):
        if not showOkCancel('Confirm', 'Overwrite "%s"?' % fileName, parent=self):
          return
      
      if self.exportCurrentModelsAction.isChecked():
        models = self.nuc.getDisplayedModels(structure)
      else:
        models = None # All available
      
      self.setCursor(QtCore.Qt.WaitCursor)
      
      counts = self.nuc.exportCoords(filePath, format, structure=structure, models=models)
      
      self.unsetCursor()  
      
      if counts:
        nChromos, nModels, nCoords = counts
        showInfo('Info', 'Exported %d coordinates for %d model(s) of %d chromosomes' % (nCoords, nModels, nChromos))


  def exportDataTrack(self, obj):
    
    typ, code = obj
    
    msg = 'Select data track export file name'
    
    options = (('BED track','BED'),
               ('Wiggle track','WIG'),
               ('BedGraph track','bedGraph'),
               ('ENCODE BroadPeak','broadPeak'),
               ('Nuc Text','TSV'),
               ('HSF5 binary ','HDF5'),
               ('Nuc JSON','JSON'),
               ('NumPy binary','NDArray'))
                 
    fileTypes = []
    for label, key in options:
      fileTypes.append( FileType(label, DATA_TRACK_FORMATS[key]) )
    
    fileName = code + DATA_TRACK_FORMATS['BED'][0][1:]

    dialog = FileDialog(self, msg, directory=self.dirPathData,
                        doSave=True, default=fileName,
                        fileTypes=fileTypes)
  
    filePath = dialog.getFile()
   
    if filePath:
      fileRoot, fileExt = splitExtension(filePath)
      filePat = '*' + fileExt.lower()
      
      if filePat == '*.broadpeak':
        filePat = '*.broadPeak'
        
      if filePat == '*.broadpeak.gz':
        filePat = '*.broadPeak.gz'
      
      for format in DATA_TRACK_FORMATS:
        if filePat in DATA_TRACK_FORMATS[format]:
          break
      
      else: # Rename
        selected = dialog.selectedNameFilter()
        
        for label, key in options:
          if selected.startswith(label):
            format = key
            fileExt = DATA_TRACK_FORMATS[format][0][1:]
            filePath = fileRoot + fileExt
            break
      
      self.dirPathData = path.dirname(filePath)
         
      if pathExists(filePath):
        if not showOkCancel('Confirm', 'Overwrite "%s"?' % fileName, parent=self):
          return

      self.setCursor(QtCore.Qt.WaitCursor)
      
      counts = self.nuc.exportDataTrack(typ, code, filePath, format)
      
      self.unsetCursor()  
      
      if counts:
        nChromos, nRegions = counts
        showInfo('Info', 'Exported %d values for %d chromosomes' % (nRegions, nChromos))


  def exportInteractions(self, code):
    
    msg = 'Select interactions export file name'
    
    options = (('Nuc Text','TSV'),
               ('HSF5 binary ','HDF5'),
               ('Nuc JSON','JSON'),
               ('NumPy binary','NDArray'))
                 
    fileTypes = []
    for label, key in options:
      fileTypes.append( FileType(label, INTERACTIONS_FORMATS[key]) )
    
    fileName = code + '.tsv'

    dialog = FileDialog(self, msg, directory=self.dirPathData,
                        doSave=True, default=fileName,
                        fileTypes=fileTypes)
  
    filePath = dialog.getFile()
   
    if filePath:
      fileRoot, fileExt = splitExtension(filePath)
      filePat = '*' + fileExt.lower()
      
      for format in INTERACTIONS_FORMATS:
        if filePat in INTERACTIONS_FORMATS[format]:
          break
      
      else: # Rename
        selected = dialog.selectedNameFilter()
        
        for label, key in options:
          if selected.startswith(label):
            format = key
            fileExt = INTERACTIONS_FORMATS[format][0][1:]
            filePath = fileRoot + fileExt
            break
      
      self.dirPathData = path.dirname(filePath)
         
      if pathExists(filePath):
        if not showOkCancel('Confirm', 'Overwrite "%s"?' % fileName, parent=self):
          return

      self.setCursor(QtCore.Qt.WaitCursor)
      
      counts = self.nuc.exportInteractions(code, filePath, format)
      
      self.unsetCursor()  
      
      if counts:
        nChromos, nRegions = counts
        showInfo('Info', 'Exported %d interactions for %d chromosomes' % (nRegions, nChromos))
       
       
  def removeDataTrack(self, obj):
    
    typ, code = obj
    nuc = self.nuc
    msg = 'Really remove "%s" dataset?' % (code,) 
    
    if showOkCancel('Confirm', msg, parent=self):
      nuc.removeDataTrack(typ, code)
      self.updateContents()


  def removeInteractions(self, code):
    
    nuc = self.nuc
    msg = 'Really remove "%s" interactions?' % (code,) 
    
    if showOkCancel('Confirm', msg, parent=self):
      nuc.removeInteractions(code)
      self.updateContents()
  
  
  def _exportWidgetImage(self, filePath, widget):
  
    if self.nuc:
    
      if not filePath:
        msg = 'Image file'
        fileTypes = [FileType('PNG', ['*.png']),]
        filePath = selectSaveFile(self, msg, fileTypes=fileTypes)
        self.update()
        time.sleep(0.2)
    
      if filePath:
        widget.update()
        pixmap = QtGui.QPixmap.grabWindow(widget.winId())
        pixmap.save(filePath, 'PNG')

       
  def exportContactImage(self, filePath=None):
    
    self._exportWidgetImage(filePath, self.contactMapPanel.contactMap)
  
  
  def exportStructureImage(self, filePath=None):
  
    self._exportWidgetImage(filePath, self.structurePanel)
  
  
  def saveRotationMovFrames(self, filePath=None, angle=360.0, axis=[0,1,0],
                            nframes=360, fps=30, codec='mpeg4'):
    """
    Create movie by rotating about a given axis through a given angle,
    spliting into specified number of frames.
    Defaults to 30 frames per second.
    """
    import subprocess, os
    
    dAngle = angle/float(nframes)
    widget = self.structurePanel
    
    w = widget.width()
    h = widget.height()
    
    rootName, ext = os.path.splitext(filePath)
    tempPath = rootName + '_%s_raw' % time.time() + ext
    
    # TBD: movie writer class
    
    cmd =['avconv', '-y', # No propt for overwrite
          '-s', '%dx%d' % (w, h),
          '-r', str(fps),
          '-frames:d', str(nframes),
          '-t', '%.3f' % (nframes/fps), 
          '-an', # Disable audio, -stats 
          '-codec:v', 'rawvideo',
          '-f', 'rawvideo', # force i/o format (container not codec)
          '-pix_fmt', 'rgba',
          '-i', '-', tempPath]
   
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
 
    for i in range(0, nframes+1):
      # Rotate structure display
      widget.rotateView(dAngle, axis)
       
      widget.update()
    
      qimage = widget.grabFrameBuffer(True)

      proc.stdin.write(qimage.bits())
   
    proc.stdin.close()
    proc.wait()
 
    cmd =['avconv', '-y', '-i', tempPath, '-codec:v', codec, filePath]
    proc = subprocess.Popen(cmd)
    proc.wait()
    
    os.unlink(tempPath)
    
    
  def saveExpansionMovFrames(self, filePath=None, angle=360.0, axis=[0,1,0], expand_from=1.0, expand_to=2.5,
                             nrot_frames=350, nexp_frames=40, fps=30, codec='mpeg4'):
    """
    Expand the nucleus by shifting the chromosome centre-of-masses axially away/towards the origin,
    and save nframes images, starting numbering from start_frame.
    Hint: assume 30 frames per second for a good resolution movie.
    """
    
    import subprocess, os
    
    nframes = 2 * (nrot_frames + nexp_frames)
    dAngle = angle/float(nrot_frames)
    widget = self.structurePanel
    
    w = widget.width()
    h = widget.height()
    
    rootName, ext = os.path.splitext(filePath)
    tempPath = rootName + '_%s_raw' % time.time() + ext
    
    # TBD: movie writer class
    
    cmd =['ffmpeg', '-y', # No propt for overwrite
          '-s', '%dx%d' % (w, h),
          '-r', str(fps),
          '-t', '%.3f' % (nframes/fps), 
          '-an', # Disable audio, -stats 
          '-codec:v', 'rawvideo',
          '-f', 'rawvideo', # force i/o format (container not codec)
          '-pix_fmt', 'bgra',
          '-i', '-',
          '-frames:d', str(nframes),
          tempPath,
          ]
   
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    
    print(cmd)
    print(nframes)
        
    expstep = (expand_to - expand_from) / float(nexp_frames)
    
    for i in range(0, nexp_frames): # Expand
      widget.parent().setExpansion(expand_from + i * expstep )
      widget.update()      
      qimage = widget.grabFrameBuffer(True)
      print(111, i, nexp_frames)
      proc.stdin.write(qimage.bits())
       
    for i in range(0, nrot_frames): # Rotate
      widget.rotateView(dAngle, axis)
      widget.update()      
      qimage = widget.grabFrameBuffer(True)
      print(222, i, nrot_frames)
      proc.stdin.write(qimage.bits())

    expstep = (expand_from - expand_to) / float(nexp_frames)
    
    for i in range(0, nexp_frames): # Contract
      widget.parent().setExpansion(expand_to + i * expstep )
      widget.update()      
      qimage = widget.grabFrameBuffer(True)
      print(333, i, nexp_frames)
      proc.stdin.write(qimage.bits())

    for i in range(0, nrot_frames): # Rotate
      widget.rotateView(dAngle, axis)
      widget.update()
      qimage = widget.grabFrameBuffer(True)
      print(444, i, nrot_frames)
      proc.stdin.write(qimage.bits())   

    proc.stdin.close()
    proc.wait()
 
    cmd =['ffmpeg', '-y', '-i', tempPath, '-codec:v', codec, filePath]
    proc = subprocess.Popen(cmd)
    proc.wait()
    
    os.unlink(tempPath)


  def backbone_trace_movie_images(self, filePath, chromo=None, section_length=20, section_centre_from=0.0,
                                  section_centre_to=1.0, nframes=101, start_frame=0):
    """
    Trace the backbone of a chromosome using a section of length section_length
    with its centre sweeping from section_centre_from to section_centre_to,
    and save nframes images, starting numbering from start_frame.
    Hint: assume 30 frames per second for a good resolution movie.
    """

    sectionstep = (section_centre_to - section_centre_from) / float(nframes - 1)
    frame = start_frame

    #chromosome
    if chromo in self.nuc.getChromosomes():
      #print visible chromosomes for lots of display section centres
      if self.nuc.getChromoDisplayParams(chromo)[0] and self.nuc.getChromoDisplayParams(chromo)[2] != 5: #isShown and not Faint
        for i in range(0, nframes):
          section_centre = section_centre_from + i * sectionstep
          self.nuc.setChromoSectionDisplayParams(chromo, section_centre=section_centre, section_numres=section_length)
          #redraw chromosomes
          self.structureOuterPanel.refreshStructure()
          frame = self.writeMovieFrame(self.structurePanel, filePath, frame)

    return frame
  
  
  def _selectMoveFrameFileRoot(self, default=None):
                          
    msg = 'Select move frame image root file name'
    fileTypes = [FileType('PNG', ['*.png']),]
    filePath = selectSaveFile(self, msg, self.dirPathMovie,
                              fileTypes, default)

    if not filePath: # Cancelled
      return
    
    # Remove any file extension 
    fileRoot = path.splitext(filePath)[0]
    
    # Remember last movie frame directory
    self.dirPathMovie = path.dirname(filePath)
    
    # Brief pause so file dialog is gone
    # by the time any images are grabbed
    self.update()
    time.sleep(0.2)
    
    return fileRoot
  

  def writeMovieFrame(self, widget, fileRoot, frame):
    """
    Save a single numbered image of a Qt widget
    """
     
    filePath = '%s_%03d.png' % (fileRoot, frame)
    
    widget.update()
    
    img = widget.grabFrameBuffer(True)
    img.save(filePath, 'PNG')

    return frame+1


  def adjustFrameSize(inputFunc): # Self here?
    """
    Decorator to wrap the movie frame creating functions
    Temporarily sets the widget size and hence output frame size
    """    
    
    def wrappedFunc(self, *args, **kw):
      
      actionGroup = self.movieMenu._groupDict[0]
      size = actionGroup.checkedAction().data()

      if size:
        prevGeom = self.structurePanel.geometry()
        self.structurePanel.setGeometry(0,0,size[0],size[1])
      
      result = inputFunc(self, *args, **kw)
      
      if size:
        self.structurePanel.setGeometry(prevGeom)
      
      return result
    
    return wrappedFunc


  @adjustFrameSize
  def exportChromosomeMovieImages(self, filePath=None):
    #sweep along chromosome section centre, and save each image from 0.00 to 1.00 every 0.01, saved as image_%3d.png
  
    if self.nuc:
    
      if not filePath:
        filePath = self._selectMoveFrameFileRoot()
        
      if filePath:  
        for chromo in self.nuc.getChromosomes():
          self.backbone_trace_movie_images(filePath+chromo, chromo=chromo, section_length=20,
                                           section_centre_from=0.0, section_centre_to=1.0,
                                           nframes=101, start_frame=0)

  @adjustFrameSize
  def exportRotationMovieImages(self, filePath = None):
    #rotate the molecule around the vertical axis, and save images at every 1 degree change.
    
    if self.nuc:
      
      if not filePath:
        msg = 'Select move file name'
        fileTypes = [FileType('MOV', ['*.mov']),] # Add AVI etc in future
        filePath = selectSaveFile(self, msg, self.dirPathMovie,
                                  fileTypes, 'StructureRotation.mov')
 
        if not filePath: # Cancelled
          return
 
        # Remove any file extension
        fileRoot = path.splitext(filePath)[0]
 
        # Remember last movie frame directory
        self.dirPathMovie = path.dirname(filePath)
 
        # Brief pause so file dialog is gone
        # by the time any images are grabbed
        self.update()
        time.sleep(0.2)
     
      if filePath:
        self.saveRotationMovFrames(filePath, angle=360.0, axis=[0,1,0], nframes=360)

  
  @adjustFrameSize
  def exportExpansionMovieImages(self, filePath=None):
    #sweep along expansion parameter, and save each image from 1.0 to 5.00 every 0.1, saved as image_%3d.png
    
    if self.nuc:
    
      if not filePath:
        msg = 'Select move file name'
        fileTypes = [FileType('MOV', ['*.mov']),] # Add AVI etc in future
        filePath = selectSaveFile(self, msg, self.dirPathMovie,
                                  fileTypes, 'StructureExpansion.mov')
 
        if not filePath: # Cancelled
          return
 
        # Remove any file extension
        fileRoot = path.splitext(filePath)[0]
 
        # Remember last movie frame directory
        self.dirPathMovie = path.dirname(filePath)
 
        # Brief pause so file dialog is gone
        # by the time any images are grabbed
        self.update()
        time.sleep(0.2)
    
      if filePath:
        self.saveExpansionMovFrames(filePath)
  
                                             
  @adjustFrameSize
  def exportAnnealMovieImages(self, filePath=None):
    """
    Export structure annealing movie frames
    Annealing uses current GUI settings
    """
    
    if self.nuc:
     
      msg = 'Recalculate structure for movie creation using current settings?'
      if not showOkCancel('Confirm', msg, parent=self):
        return
    
      if not filePath:
        filePath = self._selectMoveFrameFileRoot(default='AnnealFrame')
        
      if filePath:
        self._currentMovieFrame = [filePath, 0]
        self._movieRefreshStructure() # Starting frame
        
        panel = self.strucCalcPanel
        panel.updateFunc = self._movieRefreshStructure
        panel.startCalculation()
        panel.updateFunc = self.structureOuterPanel.refreshStructure
    
    
  def _movieRefreshStructure(self):
    """
    Temporary structure display refresh function
    which also saves a single movie frame.
    Can be used as callback while structure calculation occurs
    """
    
    fileRoot, frame = self._currentMovieFrame
    self.structureOuterPanel.refreshStructure()
    self.writeMovieFrame(self.structurePanel, fileRoot, frame)
    self._currentMovieFrame[1] += 1
  
  def createRecordingDataTrack(self):
     
     source = DERIVED
     code = 'Recording'
    
     self.nuc.setDataTrack(code, source, {}, {}, annoDict=None, stranded=None,
                           modelDict=None, color=(0.0, 1.0, 0.0), scale=1.0, threshold=0.0,
                           showText=True, shape=0)
     
     self.nuc.setDataTrackPeakType(source, code, 2)
                    
                    
  def createCoordDensityDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      nModels = nuc.getNumModels()
      chromos = nuc.getDisplayedChromosomes()
      
      if nModels and chromos:
        major, minor = nuc.getModelSize(0, chromos)
        nuc.calcDensity(0, chromos, major/5.0)
        
        for chromo in chromos:
          chrColorMode = nuc.getChromoDisplayParams(chromo)[2]

          if chrColorMode == 2:
            self.updateContents()
            break


  def createSpatialDensityDataTrack(self, obj):
  
    nuc = self.nuc
    
    if nuc:
      typ, code = obj

      self.setCursor(QtCore.Qt.WaitCursor)
      nuc.calcDataTrackSpatialDensity(code, typ)
    
      self.updateContents()
      self.unsetCursor()          


  def createRegionContactDistDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      thresh = askFloat('Query', 'Threshold distance?', 4.0)
      if thresh is None:
        return 
             
      self.setCursor(QtCore.Qt.WaitCursor)
      for group_name in nuc.getSelectedContactGroups():
        nuc.calcRegionViolDataTrack(group_name, threshold=thresh)
    
      self.updateContents()
      self.unsetCursor()          


  def createDepthDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      nModels = nuc.getNumModels()
      chromos = nuc.getDisplayedChromosomes()
      #models = list(range(nModels))
      
      if nModels and chromos:
        nuc.calcDepths(0, chromos, 5.0, 100)
      
      self.updateContents()  


  def createTransDepthDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      nModels = nuc.getNumModels()
      chromos = nuc.getDisplayedChromosomes()
      #models = list(range(nModels))
      
      if nModels and chromos:
        nuc.calcDepths(0, chromos, 2.50, 200, transInterface=True)
      
      self.updateContents()  


  def createRadGyrationDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      chromos = nuc.getDisplayedChromosomes()
      
      window = askInteger('Query', 'Window width (particles)?', 11)
      if not window:
        return
      
      if chromos:
        nuc.calcRadGyration(chromos, window)
      
      self.updateContents()  


  def createContactDistDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      nuc.calcContactDistances()
      
      self.updateContents()  


  def createContactDirectionDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      nuc.calcContactDirectionality()
      
      self.updateContents()  


  def createIntermingleDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      chromos = nuc.getDisplayedChromosomes()
      
      rad = askFloat('Query', 'Particle radius?', 2.0)
      if not rad:
        return
      
      if chromos:
        nuc.calcChromoIntermingling(chromos, rad)
      
      self.updateContents()  


  def createChromoDepthDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      nModels = nuc.getNumModels()
      chromos = nuc.getDisplayedChromosomes()
      #models = list(range(nModels))
      
      if nModels and chromos:
        # defaults to currently displayed structure 
        nuc.calcDepths(0, chromos, 5.0, 100,
                       separateChromos=True)
      
      self.updateContents()  
  
  
  def setCoordImage(self, code):
    
    self.coordImage = code
    self.updateContents()
  

  def setCurrentStructure(self, code):
    
    nuc = self.nuc
    
    if nuc.structure.name.split('/')[-1] != code:
      nuc.setCurrentStructure(code)
      self.updateContents()
  
  
  def createContactDataTrack(self, groupName):
    
    nuc = self.nuc
    
    if nuc:
      nuc.calcContactDataTrack(groupName, cis=True, trans=True, binSize=100000)
      self.updateContents()
      

  def createCisDataTrack(self, groupName):
    
    nuc = self.nuc
    
    if nuc:
      nuc.calcContactDataTrack(groupName, cis=True, trans=False, binSize=100000)
      self.updateContents()
   
   
  def createTransDataTrack(self, groupName):
  
    nuc = self.nuc
    
    if nuc:
      nuc.calcContactDataTrack(groupName, cis=False, trans=True, binSize=100000)
      self.updateContents()


  def createVoidDataTrack(self, groupName):
  
    nuc = self.nuc
    
    if nuc:
      nuc.calcContactVoidRegions(groupName)
      self.updateContents()


  def createBackboneDictDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      nuc.calcBackboneDistances()
      self.updateContents()
  
  
  def createLinearityDataTrack(self):
  
    nuc = self.nuc
    
    if nuc:
      regDict = {}
      valDict = {}
      coordsGroup = nuc._getCoordsGroup()
      particGroup = nuc._getParticleGroup()
      
      for chrA in coordsGroup:
        data = array(coordsGroup[chrA])
        nModels, nCoords, nDims = data.shape
        
        if not nModels: 
          continue

        if nCoords < 3:
          continue
        
        vals = []
        
        for i in range(nModels):
          coords = data[i]
          core = coords[1:-1]
          
          d1 = numpy.clip(coords[:-2] - core, -1e99, 1e99)
          d2 = numpy.clip(core - coords[2:],  -1e99, 1e99)
          
          s1 = numpy.sqrt((d1*d1).sum(axis=1))
          s2 = numpy.sqrt((d2*d2).sum(axis=1))
          
          s = numpy.clip(s1*s2, 1e-99, 1e99)
          proj = (d1*d2).sum(axis=1) / s
         
          vals.append(numpy.arccos(proj))
        
        vals = numpy.array(vals).mean(axis=0)
        pos  = numpy.array(particGroup[chrA]['positions'])
        
        regDict[chrA] = numpy.vstack([pos[:-2], pos[2:]]).T
        norm = vals-vals.min()
        norm /= norm.max() or 1.0
        valDict[chrA] = numpy.vstack([vals, norm]).T
      
      nuc.setDataTrack('bboneAngle', DERIVED, regDict,
                        valDict, color=(0.5, 0.4, 1.0))
      self.updateContents()
  
  def createDataTrackPairInteractions(self):
  
    codes = []
    
    for source in (EXTERNAL, INNATE, DERIVED):
      codes += self.nuc.getDataTrackCodes(source)
    
    if codes:
      code = askChoice('Query', 'Select data track', codes, self)
      
      if code:
        self.nuc.calcDataTrackPairIteractions(code)
        self.updateContents()  

        
  def alignModels(self):
     
     structure = self.nuc.structure.name.split('/')[-1]
     chromos = self.nuc.getDisplayedChromosomes(structure)
     
     self.nuc.modelAlign(chromosomes=chromos, structure=structure)
     self.structurePanel.update()
     self.updateContents()
 
  
  def recalcRestraints(self, ):
  
    self.strucCalcPanel.recalcRestraints()
    self.updateContents()


  def recalcDensity(self, ):
  
    self.strucCalcPanel.recalcDensity()
    self.updateContents()
  
  
  def annealHighRes(self, ):
  
    self.strucCalcPanel.fullAnneal()
    self.strucCalcPanel.startCalculation()
   
    
  def annealLowRes(self, ):

    self.strucCalcPanel.quickAnneal()
    self.strucCalcPanel.startCalculation()
   
    
  def annealMinimise(self, ):
  
    self.strucCalcPanel.minimise()
    self.strucCalcPanel.startCalculation()
  
  
  def showAllChromos(self):
  
    if self.nuc:
      allChromos = self.nuc.getChromosomes()
      for chromo in allChromos:
        self.nuc.setChromoDisplayed(chromo, True)
    
      self.updateContents()


  def showNoChromos(self):
  
    if self.nuc:
      allChromos = self.nuc.getChromosomes()
      for chromo in allChromos:
        self.nuc.setChromoDisplayed(chromo, False)
    
      self.updateContents()
  
  
  def showOneChromo(self):
  
    if self.nuc:
      allChromos = self.nuc.getChromosomes()
      if not allChromos:
        return
      
      selected = self.nuc.getDisplayedChromosomes()
      if selected:
        keep = selected[0]
        
      else:
        keep = allChromos[0]
      
      for chromo in allChromos:
        if chromo == keep:
          self.nuc.setChromoDisplayed(chromo, True)
        else:  
          self.nuc.setChromoDisplayed(chromo, False)
    
      self.updateContents()
  
  
  def chromosInvertSelected(self):
  
    allChromos = self.nuc.getChromosomes()
    selected = self.nuc.getDisplayedChromosomes()
   
    for chromo in allChromos:
      if chromo in selected:
        self.nuc.setChromoDisplayed(chromo, False)
      else:  
        self.nuc.setChromoDisplayed(chromo, True)
  
    self.updateContents()
  
   
  def chromosFaintUnselected(self):
  
    setDisplay =  self.nuc.setChromoDisplayParams
    selected = self.nuc.getDisplayedChromosomes()
    
    for chromo in self.nuc.getChromosomes():
      if chromo not in selected:
        # Faint lines
        setDisplay(chromo, isShown=True, useLabels=False, colorMode=5, displayMode=1)
        
    self.updateContents()  
  
  
  def chromosResetColors(self):
    
    if self.nuc:
      self.nuc.resetChromoColors()
      self.updateContents()


  def chromosResetCyclingColors(self):
    
    if self.nuc:
      self.nuc.resetChromoColors(False)
      self.updateContents()
  

  def toggleChromosomes(self, toggled):
    
    if self.nuc:
      buttons = self.chromoButtonGroup.buttons()
      selection = set([b.text().strip() for b in buttons if b.isChecked()])
      
      if selection != set(self.nuc.getDisplayedChromosomes()):
        chromoGroup = self.nuc.chromosomes
        
        for chromo in self.nuc.getChromosomes():
          self.nuc.setChromoDisplayed(chromo, chromo in selection)
          
      self.updateContents()
      
      
  def updateChromosomeToolbar(self):
  
    if self.nuc:
      chromos = self.nuc.getChromosomes()
      buttons = self.chromoButtonGroup.buttons()
      selected = set(self.nuc.getDisplayedChromosomes())
      
      for i, chromo in enumerate(chromos):
        text = '%-2s' % chromo
        isSelected = chromo in selected
        color = self.nuc.getChromoColor(chromo).rgbHex
        icon = Icon(None, color, 12)
        tipText = 'Hide/show chromosome %s' % chromo
        
        if i < len(self.chromoActions):
          button = buttons[i]
          button.setChecked(isSelected)
          button.setText(text)
          button.setIcon(icon)
          button.setToolTip(tipText)
          self.chromoActions[i].setVisible(True)
        
        else:
          button = Button(self.chromoToolbar, text, icon=icon,
                          tipText=tipText, grid=None)
          button.setCheckable(True)
          button.setChecked(isSelected)
 
          self.chromoButtonGroup.addButton(button)
          self.chromoButtonGroup.setId(button, i)
          action = self.chromoToolbar.addWidget(button)
          self.chromoActions.append(action)
      
      a = len(self.chromoActions)
      c = len(chromos)
      
      if a > c:
        for i in range(c, a):
          self.chromoActions[i].setVisible(False)
  
        
  def toggleDataTrack(self, toggled):
    
    nuc = self.nuc
     
    if nuc:
      buttons = self.dataTrackButtonGroup.buttons()
      selection = set([b.obj for b in buttons if b.isChecked()])
      
      interactions = set(nuc.getDisplayedInteractions())
      data_tracks = set(nuc.getDisplayedDataTracks())
      
      if selection != interactions | data_tracks:
        for code in self.nuc.interactions:
          isShown = code in selection
          nuc.setInteractionsDisplayed(code, isShown)
 
        for typ in EXTERNAL, INNATE, DERIVED:
          for code in self.nuc.dataTracks[typ]:
            isShown = code in selection
            nuc.setDataTrackDisplayed(typ, code, isShown)
        
        self.updateContents()
      
      
  def updateDataTrackToolbars(self):
    
    nuc = self.nuc
    
    if self.nuc:
      selected = set([x[1] for x in nuc.getDisplayedDataTracks()])
      selectedInter = set(nuc.getDisplayedInteractions())
 
      butGroup = self.dataTrackButtonGroup
 
      frame = self.dataTrackFrame
      codes = sorted([c for c in nuc.interactions])
      group_data = [(c.lower(), c, nuc.interactions[c], True) for c in codes]
      
      for typ in   EXTERNAL, INNATE, DERIVED:
        codes = [c for c in nuc.dataTracks[typ]]
        group_data += [(c.lower(), c, nuc.dataTracks[typ][c], False) for c in codes]
      
      group_data.sort()

      i = 0
      for sk, code, group, isInter in group_data:
        if 'display' not in group.attrs:
          continue
 
        if not len(group):
          continue
 
        color = numpy.array(group.attrs['display'][:3])
        color = '#%02X%02X%02X' % tuple(255*color)
        icon = Icon(None, color, 12)
 
        if isInter:
          isSelected = code in selectedInter
          tipText = 'Hide/show interactions "%s"' % (code)
          obj = code
 
        else:
          isSelected = code in selected
          tipText = 'Hide/show data track "%s"' % (code)
          obj = code
 
        if i < len(butGroup.buttons()):
          button = butGroup.buttons()[i]
          button.setChecked(isSelected)
          button.setText(code)
          button.setIcon(icon)
          button.setToolTip(tipText)
          button.obj = obj
          button.setVisible(True)
 
        else:
          button = Button(self.dataTrackFrame, code, icon=icon,
                          tipText=tipText, grid=None)
          button.setCheckable(True)
          button.setChecked(isSelected)
          button.obj = obj
          
          #font = button.font()
          #font.setPointSize(8)
          #button.setFont(font)
 
          butGroup.addButton(button)
          butGroup.setId(button, i)
 
        i += 1
     
      buttons = butGroup.buttons()
      n = len(buttons)
      if n > i:
        for j in range(i, n):
          buttons[j].setVisible(False)
            

  def _setColorMode(self, value):
    
    refresh = False
    
    nuc = self.nuc
    for chromo in nuc.getChromosomes():
      if nuc.setChromoDisplayParams(chromo, colorMode=value):
        refresh = True
    
    if refresh:
      self.updateContents()
 
 
  def _setDisplayMode(self, value):
    
    refresh = False
    
    nuc = self.nuc
    for chromo in nuc.getChromosomes():
      if nuc.setChromoDisplayParams(chromo, displayMode=value):
        refresh = True
    
    if refresh:
      self.updateContents()
 

  def toggleRestraints(self, *args):
    
    if self.nuc:
      attrs = self.nuc.display.attrs
      opts = list(attrs['options'])      
      showCis    = opts[2]
      showTrans  = opts[3]   
      
      if showCis or showTrans:
        self.nuc.setRestraintsDisplayed(False, False)
        self.restActionA.setChecked(False)
        self.restActionB.setChecked(False)
      
      else:
        self.nuc.setRestraintsDisplayed(True, True)
        self.restActionA.setChecked(True)
        self.restActionB.setChecked(True)
       
      self.updateContents()
      
  def _graphPseudo4C(self):
     
    groupNames = self.nuc.getSelectedContactGroups()
    
    if groupNames:
      chromo = askString('Query', 'Chromosome:', '1', parent=self)
      
      if chromo not in self.nuc.getChromosomes():
        msg = 'Chromosome "{}" not present'.format(chromo)
        showWarning('Warning', msg, self)
        return
      
      pos =  askInteger('Query', 'Chromosome position (bp):', int(92e6))
      if not pos:
        return
      
      Graphs.graphPseudo4c(self.nuc, chromo, pos, groupNames)
     
     
  def _graphContactFourierTransform(self):
  
    for groupName in self.nuc.getSelectedContactGroups():
      Graphs.graphContactFourierTransform(self.nuc,  groupName=groupName)


  def _graphContactSeqSep(self):
    
    for groupName in self.nuc.getSelectedContactGroups():
      Graphs.graphContactSeqSep(self.nuc, groupName)
  
  
  def _graphViolations(self):
   
    Graphs.graphViolations(self.nuc)
    
    
  def _graphContactDistances(self):
        
    for groupName in self.nuc.getSelectedContactGroups():
      Graphs.graphContactDistances(self.nuc, groupName, structure=self.nuc.structure.name.split('/')[-1])


  def _graphIntermingling(self):
  
    chromos = self.nuc.getDisplayedChromosomes()
    if not chromos:
      msg = 'No chromosomes selected'
      showWarning('Warning', msg, self)
      return
    
    Graphs.graphChromoIntermingling(self.nuc, chromos)
    
    
  def _graphTransSpacing(self): 
    
    for groupName in self.nuc.getSelectedContactGroups():
      Graphs.graphTransSpacing(self.nuc, groupName)


  def _clusterDataTrack3d(self):
    
    from analyses import SpatialClustering
    
    keys = self.nuc.getDisplayedDataTracks()
    
    if not keys:
      return
      
    """
    if len(keys) == 1:
      typ, code = keys[0]
    
    else: 
      texts = [k[1] for k in keys]
      msg = 'Select data track to analyse'
      typ, code = showMulti('Select data track', msg, texts, objects=keys, parent=self)
    """
    SpatialClustering.dbScanClusterDataTrack(self.nuc)


  def _clusterStructures(self): 
    
    if self.nuc.structure:
      linkage = self.nuc.clusterStructures([self.nuc.structure.name.split('/')[-1]])
     
      from matplotlib import pyplot
      from scipy.cluster import hierarchy
      
      pyplot.figure()
      pyplot.title('Structure hierarchy')
      graph = hierarchy.dendrogram(linkage)
      pyplot.show()
      
