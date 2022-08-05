from sys import path
from time import time
from math import ceil
from numpy import array, float32, int32, ones, zeros, hstack
from numpy import dstack, eye, log, uint8, uint32, empty, clip, concatenate

from gui.GenomeBrowserPanel import GenomeBrowserPanel

from gui.qtgui.Base import Icon
from gui.qtgui.ToolBar import ToolBar
from gui.qtgui.Button import Button
from gui.qtgui.Menu import Menu
from gui.qtgui.SpinBox import FloatSpinBox
from gui.qtgui.Slider import FloatSlider
from gui.qtgui.LabelFrame import LabelFrame
from gui.qtgui.Entry import IntEntry
from gui.qtgui.Label import Label
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.Colors import ColorDialog, GradientEditor

from PySide import QtCore, QtGui, QtOpenGL
from cUtil.apiUtil import binContacts
from cUtil.drawing import addContactMapPoints, addContactMapBlocks
from cUtil.drawing import addContactMapMesh, addContactMapRegions, getContactMapRegions
from OpenGL import GL, GLU, GLUT

LeftButton = QtCore.Qt.LeftButton
MiddleButton = QtCore.Qt.MiddleButton
RightButton = QtCore.Qt.RightButton
QPoint = QtCore.QPoint
Qt = QtCore.Qt

class ContactMapPanel(QtGui.QWidget):

  def __init__(self, parent, mainApp, openFunc, **kw):
    
    QtGui.QWidget.__init__(self, parent)
   
    layout = QtGui.QGridLayout(self)
    layout.setRowStretch(0, 0)
    layout.setRowStretch(1, 2)
    layout.setRowStretch(2, 0)
    layout.setColumnStretch(0, 0)
    layout.setColumnStretch(1, 2)
    
    layout.setSpacing(0)
    layout.setContentsMargins(0,0,0,0) 

    self.xRegionText = ''
    self.yRegionText = ''
    self.mainApp = mainApp

    # Contacts toolbar
     
    icons = ['zoom-in.png', 'zoom-out.png',
             'zoom-center.png', 'move-center.png']
    
    icons = [mainApp.getIcon(i) for i in icons]
    
    texts = ['Zoom in','Zoom out','Reset zoom','Center']
    
    funcs = [self.zoomIn, self.zoomOut, self.zoomReset, self.moveCenter]
   
    self.contactToolbar = ToolBar(self, 'Contacts toolbar', funcs, icons, texts, 
                                  objName='contactToolbar', areas='tblr', iconSize=32)
    
    self.contactGroupButton = Button(self.contactToolbar, '', iconSize=32,
                                      icon=mainApp.getIcon('contacts.png'),
                                      tipText='Set contact dataset')
    self.contactGroupMenu = Menu(self.contactGroupButton, 'Select contact dataset to display', 
                                 setupFunc=self._setupContactGroupMenu) 
                                                                 
                                                        
    self.contactGroupButton.setMenu(self.contactGroupMenu)
    self.contactToolbar.addWidget(self.contactGroupButton)

    colorConfigButton = Button(self.contactToolbar, '', icon=mainApp.getIcon('colors.png'),
                               tipText='Colour settings', iconSize=32)
    colorConfigMenu = Menu(colorConfigButton, 'Colours configure',
                           setupFunc=self._setupContactColorConfigMenu)
    colorConfigButton.setMenu(colorConfigMenu)
    self.contactToolbar.addWidget(colorConfigButton)
   
    contactConfigButton = Button(self.contactToolbar, '',
                                 tipText='Map brightness (gamma) correction',
                                 icon=mainApp.getIcon('brightness.png'), iconSize=32)
    contactConfigMenu = Menu(contactConfigButton, 'Gamma adjustment')
    self.gammaWidget =  FloatSlider(contactConfigMenu, startVal=1.0, endVal=-1.0, value=0.0,
                                    direction='h', step=0.05, bigStep=None, callback=self._gammaUpdateProxy)                      
    contactConfigMenu.addItem(u'log\u2081\u2080(\u03B3)', widget=self.gammaWidget)
   
    contactConfigButton.setMenu(contactConfigMenu)
    self.contactToolbar.addWidget(contactConfigButton)

    contactConfigButton = Button(self.contactToolbar, '  ',
                                 tipText='Configure contact map and data tracks',
                                 icon=mainApp.getIcon('configure.png'), iconSize=32)
    contactConfigMenu = Menu(contactConfigButton, 'Data track width',
                             setupFunc=self._setupContactOptionsMenu)
    contactConfigButton.setMenu(contactConfigMenu)
    self.contactToolbar.addWidget(contactConfigButton)

    # Coords
      
    self.contactRegionButton = Button(self.contactToolbar, '', grid=(0,1))
    self.contactRegionButton.setMaximumSize(250, 44)
    contactRegionMenu = Menu(self.contactRegionButton, 'Contact map region',
                             setupFunc=self._setupContactRegionMenu)
    self.contactRegionButton.setMenu(contactRegionMenu)
    self.contactToolbar.addWidget(self.contactRegionButton)
    
    self.pixmapCache = {}
    self.genomeBrowserX = GenomeBrowserPanel(self, mainApp, self.showContactXRegion)
    self.genomeBrowserY = GenomeBrowserPanel(self, mainApp, self.showContactYRegion, isVertical=True)
    self.contactMap     = ContactMapWidget(self, mainApp, openFunc) 
    
    corner = QtGui.QLabel(' ', self,)
    corner.setStyleSheet("background-color: rgb(64, 64, 64);")
    corner.setAutoFillBackground(True)
    
    layout.addWidget(self.contactToolbar, 0, 0, 1, 2)
    layout.addWidget(self.genomeBrowserY, 1, 0)
    layout.addWidget(self.contactMap, 1, 1)
    layout.addWidget(corner, 2, 0)
    layout.addWidget(self.genomeBrowserX, 2, 1)
    
    self.setLayout( layout )
    #self.setStyleSheet("color: white; background-color: black; selection-color: yellow; selection-background-color: rbg(0,0,64);")
    self.setAutoFillBackground(True)
    self.setAcceptDrops(True) 
  
  
  def _gammaUpdateProxy(self, value):
  
    self.contactMap.draw()

    
  def showContactXRegion(self, xChromo1, xSeq1, xChromo2, xSeq2):
    
    xSeqB1 = xSeq1 % 1000
    xSeqK1 = ((xSeq1 % 1000000) - xSeqB1) / 1000
    xSeqM1 = (xSeq1 - xSeqK1) / 1000000
    
    xSeqB2 = xSeq2 % 1000
    xSeqK2 = ((xSeq2 % 1000000) - xSeqB2) / 1000
    xSeqM2 = (xSeq2 - xSeqK2) / 1000000
    
    #if xChromo1 == xChromo2:
    #  data = (xChromo1, xSeqM1, xSeqK1, xSeqB1, xSeqM2, xSeqK2, xSeqB2)
    #  text = '%s:%3.3d %3.3d %3.3d - %3.3d %3.3d %3.3d' % data
    
    data = (xChromo1, xSeqM1, xSeqK1, xSeqB1, xChromo2, xSeqM2, xSeqK2, xSeqB2)
    text = '%s:%3.3d %3.3d %3.3d - %s:%3.3d %3.3d %3.3d' % data
    
    self.xRegionText = text 
    self.contactRegionButton.setText('%s\n%s' % (self.xRegionText, self.yRegionText))   
    

  def showContactYRegion(self, yChromo1, ySeq1, yChromo2, ySeq2):

    ySeqB1 = ySeq1 % 1000
    ySeqK1 = ((ySeq1 % 1000000) - ySeqB1) / 1000
    ySeqM1 = (ySeq1 - ySeqK1) / 1000000

    ySeqB2 = ySeq2 % 1000
    ySeqK2 = ((ySeq2 % 1000000) - ySeqB2) / 1000
    ySeqM2 = (ySeq2 - ySeqK2) / 1000000
    
    #if yChromo1 == yChromo2:
    #  data = (yChromo1, ySeqM1, ySeqK1, ySeqB1, ySeqM2, ySeqK2, ySeqB2)
    #  text = '%s:%3.3d %3.3d %3.3d - %3.3d %3.3d %3.3d' % data
    
    data = (yChromo1, ySeqM1, ySeqK1, ySeqB1, yChromo2, ySeqM2, ySeqK2, ySeqB2)
    text = '%s:%3.3d %3.3d %3.3d - %s:%3.3d %3.3d %3.3d' % data

    self.yRegionText = text 
    self.contactRegionButton.setText('%s\n%s' % (self.xRegionText, self.yRegionText))   


  def _toggleContactGroup(self, obj):
  
    groupName, isSelected, diagonalSide = obj
    self.mainApp.nuc.setChromoGroupSelected(groupName, not isSelected)
    
    if diagonalSide:
      self.mainApp.nuc.getContactGroup(groupName).attrs['diagonalSide'] = diagonalSide
    
    self.updateContents()
       
  
  def _setupContactGroupMenu(self, menu):
    
    menu.clear()
    nuc = self.mainApp.nuc
    
    if nuc:
      
      selected = set(nuc.getSelectedContactGroups())
      
      if self.contactMap.splitMap:
        label = Label(menu, ' <b>Left:</b>')
        menu.addItem('', widget=label)
        for name in nuc.getContactGroupNames():
          isSelected = (name in selected) and nuc.getContactsDiagonalSide(name) > 0
          menu.addItem(name, checked=isSelected,
                       callback=self._toggleContactGroup,
                       object=(name, isSelected, 1))

        label = Label(menu, ' <b>Right:</b>')
        menu.addItem('', widget=label)
        for name in nuc.getContactGroupNames():
          isSelected = (name in selected) and nuc.getContactsDiagonalSide(name) < 0
          menu.addItem(name, checked=isSelected,
                       callback=self._toggleContactGroup,
                       object=(name, isSelected, -1))
      
      else:
        for name in nuc.getContactGroupNames():
          isSelected = (name in selected)
          menu.addItem(name, checked=isSelected,
                       callback=self._toggleContactGroup,
                       object=(name, isSelected, None))


  def zoomIn(self):

    self.zoom(0.66666666667)
    
    
  def zoomOut(self):
  
    self.zoom(1.5)
    
    
  def zoomReset(self):
  
    self.resetZoom()
  
  
  def moveCenter(self):
  
    self.resetView()
  
  
  def toggleFlipY(self, action):
  
    self.contactMap.flipY = action.isChecked()
    self.updateContents()


  def toggleSplitMap(self, action):
  
    self.contactMap.splitMap = action.isChecked()
    self.updateContents()
   

  def _setupContactOptionsMenu(self, menu):
  
    menu.clear()
    nuc = self.mainApp.nuc
    
    if nuc:
      trackWidth = self.genomeBrowserX.trackWidth
      widget =  FloatSpinBox(menu, trackWidth, minValue=10.0, maxValue=200.0,
                             step=5.0, callback=self._setDataTrackWidth)                       
      menu.addItem('Genome track width:', widget=widget)
      menu.addItem('Switch direction', callback=self.toggleFlipY,
                   checked=self.contactMap.flipY)
      menu.addItem('Split along diagonal', callback=self.toggleSplitMap,
                   checked=self.contactMap.splitMap)


  def _setupContactRegionMenu(self, menu):
  
    menu.clear()
    nuc = self.mainApp.nuc
    
    if nuc:
      chromosomes = nuc.getDisplayedChromosomes()
      if not chromosomes:
        return
      
      seqLen = 0
      chromoRegions = {}
      for chromo in chromosomes:
        start, end = nuc.getChromosomeLimits(chromo)
        size = end-start
        if not size:
          continue
        
        chromoRegions[chromo] = (seqLen, seqLen+size)
        seqLen += size
        
      fx = sum(self.genomeBrowserX.viewRegion)/2.0
      fy = sum(self.genomeBrowserY.viewRegion)/2.0
      px = int(fx * seqLen)
      py = int(fy * seqLen) 
      bx = 0
      by = 0
      cx = chromosomes[0]
      cy = chromosomes[0]
      
      for chromo in chromoRegions:
        start, end = chromoRegions[chromo]
        
        if start <= px < end:
          bx = px-start
          cx = chromo
        
        if start <= py < end:
          by = py-start
          cy = chromo
      
      ix = chromosomes.index(cx)
      iy = chromosomes.index(cy)

      chrNames = ['Chr %s' % x for x in chromosomes]
      frame = LabelFrame(menu, 'Set view centre:')
      menu.addItem('', widget=frame)
      
      Label(frame, 'Horiz:', grid=(0,0))
      self._contactRegSelectChrX = PulldownList(frame, chrNames, chromosomes,
                                                callback=self._checkContactRegion,
                                                index=ix, grid=(0,1))
      self._contactRegSelectPosX = IntEntry(frame, bx, grid=(0,2),
                                            callback=self._setContactRegion)                    

      Label(frame, ' Vert:', grid=(1,0))
      self._contactRegSelectChrY= PulldownList(frame, chrNames, chromosomes,
                                               callback=self._checkContactRegion,
                                               index=iy, grid=(1,1))                       
      self._contactRegSelectPosY = IntEntry(frame, by, grid=(1,2),
                                            callback=self._setContactRegion)                         

  def _setDensityRadius(self, value):
        
    nuc = self.mainApp.nuc
    nuc.setDisplaySizes(density=value)
    

  def _setContactGroupColor(self, groupName):
   
    nuc = self.mainApp.nuc
    group = nuc.getContactGroup(groupName)
    prev  = nuc.getContactGroupColor(groupName).rgb
    
    colorObj = ColorDialog(self).getColor(prev)
            
    if colorObj:
      if nuc.setContactGroupColor(groupName,  colorObj.getRgbF()):
        self.updateContents()


  def _setupContactColorConfigMenu(self, menu):
  
    menu.clear()
    nuc = self.mainApp.nuc
    
    if nuc:    
      
      for groupName, group in nuc.getContactGroups():
      
        color = nuc.getContactGroupColor(groupName).rgbHex
        icon = Icon(None, color, 12)
        
        menu.addItem(groupName, icon=icon, object=groupName,
                     callback=self._setContactGroupColor)
  
  
  def _setDataTrackWidth(self, value):

    self.genomeBrowserX.trackWidth = value
    self.genomeBrowserY.trackWidth = value
    self.updateContents()
    

  def _checkContactRegion(self, *args):
    
    nuc = self.mainApp.nuc
    
    xChromo = self._contactRegSelectChrX.get()
    yChromo = self._contactRegSelectChrY.get()
    xPos = self._contactRegSelectPosX.get()
    yPos = self._contactRegSelectPosY.get()
    
    xStart, xEnd = nuc.getChromosomeLimits(xChromo)
    yStart, yEnd = nuc.getChromosomeLimits(yChromo)
    
    if not (xStart <= xPos < xEnd):
      xPos = (xStart+xEnd)/2
      self._contactRegSelectPosX.set(xPos, doCallback=False)
    
    if not (yStart <= yPos < yEnd):
      yPos = (yStart+yEnd)/2
      self._contactRegSelectPosY.set(yPos, doCallback=False)

      
  def _setContactRegion(self, *args):
    
    xChromo = self._contactRegSelectChrX.get()
    yChromo = self._contactRegSelectChrY.get()
    xPos = self._contactRegSelectPosX.get()
    yPos = self._contactRegSelectPosY.get()
      
    if xChromo and yChromo:
      self.setRegion(xChromo, yChromo, xPos, yPos)
    
    else:
      self.update()
    
    
  def setRegion(self, xChromo, yChromo, xPos, yPos, xWidth=None, yWidth=None):
  
    self.contactMap.setRegion(xChromo, yChromo, xPos, yPos, xWidth, yWidth)
  
  
  def translate(self, dx, dy):
  
    self.contactMap.translate(dx, dy)
    
 
  def resetView(self):
  
    self.contactMap.resetView()
    

  def resetZoom(self):
  
    self.contactMap.resetZoom()
  
  
  def zoom(self, val):
  
    self.contactMap.zoom(val)
  
  
  def updateContents(self):
  
    self.contactMap.draw()
    self.genomeBrowserX.update()
    self.genomeBrowserY.update()
    
  
class ContactMapWidget(QtOpenGL.QGLWidget):

  def __init__(self, parent, mainApp, openFunc, glList=512, **kw):
    
    QtOpenGL.QGLWidget.__init__(self, parent)
    
    self.mainApp = mainApp
    self.openFunc = openFunc
    self.glList = glList
    self.genomeBrowserX = parent.genomeBrowserX
    self.genomeBrowserY = parent.genomeBrowserY
    self.gammaWidget = parent.gammaWidget
    
    self.select_start = None
    self.viewSize = [640, 480]
    self.viewRegion = (0.0, 1.0, 0.0, 1.0) # Fraction of genome
    self.prevPos = QPoint()
    self.boxWidth = 1
    self.flipY = False
    self.splitMap = False
   
    self.blank = zeros((800,800,3), dtype=uint8)
    self.matrix = self.blank   
    
    self.setStyleSheet("background-color: rgb(0, 0, 0);")
    self.setAutoFillBackground(True)
    self.setMouseTracking(True)

    self.prevPos = QPoint()


  def setRegion(self, xChromo, yChromo, xPos, yPos, xWidth=None, yWidth=None):
    
    #print "R", xChromo, yChromo, xPos, yPos, xWidth, yWidth
    
    x1, x2, y1, y2 = self.viewRegion
    
    nuc = self.mainApp.nuc
    chromosomes = nuc.getDisplayedChromosomes()
    
    if not chromosomes:
      return
    
    x0 = 0.0
    y0 = 0.0
      
    seqLen = 0
    for chromo in chromosomes:
      start, end = nuc.getChromosomeLimits(chromo)
      size = (end-start)
      if not size:
        continue
      
      if chromo == xChromo:
        x0 = seqLen + xPos
      
      if chromo == yChromo:
        y0 = seqLen + yPos
      
      seqLen += size
    
    seqLen =float(seqLen)
    
    x0 /= seqLen
    y0 /= seqLen
    
    ar = (y2-y1)/(x2-x1)
    
    if xWidth:
      xWidth /= seqLen
    else:
      xWidth = x2-x1
    
    if yWidth:
      yWidth /= seqLen
    else:
      yWidth = y2-y1
    
    xWidth /= 2.0
    yWidth /= 2.0
    
    xWidth = max(1e3/seqLen, min(x0, 1.0-x0, xWidth))
    yWidth = max(1e3/seqLen, min(y0, 1.0-y0, yWidth))
    
    w = min(xWidth, yWidth)
    
    x1 = x0 - w 
    x2 = x0 + w
    y1 = y0 - w * ar
    y2 = y0 + w * ar
        
    self._updateView(x1, x2, y1, y2)
    
  
  def initializeGL(self):
  
    #self.draw()
    GL.glClearColor(0.2,0.2,0.2,1.0)
    
    GL.glEnable(GL.GL_TEXTURE_2D)
    texture = GL.glGenTextures(1)
    GL.glBindTexture(GL.GL_TEXTURE_2D, texture)
    GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
  
    self.initGlParams()
  
  
  def initGlParams(self):
     
    GL.glDisable(GL.GL_DEPTH_TEST)
    GL.glDisable(GL.GL_BLEND)
    GL.glDisable(GL.GL_CLIP_PLANE0)
    GL.glDisable(GL.GL_CULL_FACE)
    GL.glDisable(GL.GL_MULTISAMPLE)
    GL.glDisable(GL.GL_LIGHTING)
    GL.glDisable(GL.GL_LIGHT1)
    GL.glDisable(GL.GL_COLOR_MATERIAL)
    
    GL.glEnable(GL.GL_TEXTURE_2D)
    
    GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
    GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST)
    GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_REPEAT)
    GL.glTexParameterf(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_REPEAT)
    
    GL.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE) 


  def resizeGL(self, width, height):
    
    ar = width/float(height)
 
    x1, x2, y1, y2 = self.viewRegion

    dy = (y2-y1)
    dx = dy * ar
    
    x2 = min(1.0, x1 + dx)
    
    self.viewRegion = x1, x2, y1, y2
    self.viewSize = width, height

    
    # Whole thing gets smaller, no clipping
    GL.glViewport(0, 0, width, height)
    self.draw()
  
  
  def paintGL(self):
    
    self.initGlParams()
    
    w, h = self.viewSize
    #x1, x2, y1, y2 = self.viewRegion
    #sy, sx, sz = self.matrix.shape
    
    GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
    
    GL.glMatrixMode(GL.GL_PROJECTION)
    
    GL.glLoadIdentity()
    GL.glOrtho(0, w, h, 0, 0, 1)
   
    GL.glMatrixMode(GL.GL_MODELVIEW)

    GL.glLoadIdentity()
    GL.glDepthMask(False)
    #GL.glScalef(1.0, 1.0, 1.0)
    #GL.glTranslatef(0.0, 0.0, 0.0)
    GL.glCallList(self.glList)
  
    
  def dragEnterEvent(self, event):
    
    if event.mimeData().urls():
      event.acceptProposedAction()
        
    else:
      event.ignore()


  def dragMoveEvent(self, event):
    
    self.dragEnterEvent(event)


  def dragLeaveEvent(self, event):
  
    event.ignore()

    
  def dropEvent(self, event):
     
    if event.mimeData().urls():
      mods = event.keyboardModifiers()
      haveCtrl = mods & QtCore.Qt.ControlModifier
      haveShift = mods & QtCore.Qt.ShiftModifier
      
      if haveShift or haveCtrl:
        replace = True
      else:
        replace = False
      
      filePaths = [url.path() for url in event.mimeData().urls() if path.isfile(url.path())]
      
      if filePaths and self.openFunc:
        self.openFunc(filePaths, replace)
        event.acceptProposedAction()     
      else:
        event.ignore()
        
    else:
      event.ignore()
      
      
  def minimumSizeHint(self):
  
    return QtCore.QSize(50, 50)


  def sizeHint(self):
  
    return QtCore.QSize(400, 300) 


  def _updateView(self, x1, x2, y1, y2):
      
    x1 = max(0.0, x1)
    y1 = max(0.0, y1)
    
    x2 = min(1.0, x2)
    y2 = min(1.0, y2)
    
    if (x1, x2, y1, y2) == self.viewRegion:
      return
    
    if y2 == y1:
      raise Exception()
    self.viewRegion = (x1, x2, y1, y2)
    
    self.draw(False)
    
  
  def zoom(self, z):
    
    x1, x2, y1, y2 = self.viewRegion
    
    pos = self.mapFromGlobal(QtGui.QCursor.pos())
    fx = pos.x() / float(self.width())
    fy = pos.y() / float(self.height())
    
    gx = 1.0 - fx
    gy = 1.0 - fy
    
    dx = x2-x1
    dy = y2-y1
    
    cx = x1 + fx * dx
    cy = y1 + fy * dy
    
    dx *= z
    dy *= z    
    
    x1 = cx - (fx*dx)
    x2 = cx + (gx*dx)
    y1 = cy - (fy*dy)
    y2 = cy + (gy*dy)
    
    self._updateView(x1, x2, y1, y2)


  def resetZoom(self):
  
    self._updateView(0.0, 1.0, 0.0, 1.0)
    
    
  def resetView(self):
    # To centre
    
    x1, x2, y1, y2 = self.viewRegion
    
    dx = (x2-x1)/2.0
    dy = (y2-y1)/2.0
    
    self._updateView(0.5-dx, 0.5+dx,0.5-dy, 0.5+dy)

    
  def translate(self, tx, ty):
          
    if tx or ty:
      w, h = self.viewSize
      x1, x2, y1, y2 = self.viewRegion
      dx = x2-x1
      dy = y2-y1
      
      tx *= dx/w # pixels * value / pixel
      ty *= dy/h      
      
      x1 -= tx
      x2 -= tx
      y1 -= ty
      y2 -= ty
      
      if x1 < 0:
        x1 = 0
        x2 = dx
       
      elif x2 > 1.0:
        x2 = 1.0
        x1 = 1.0 - dx
      
      if y1 < 0:
        y1 = 0
        y2 = dy
      
      elif y2 > 1.0:
        y2 = 1.0
        y1 = 1.0 - dy  
      
      self._updateView(x1, x2, y1, y2)
    
  def moveLeft(self):

    self.translate(-10,0)
    
    
  def moveRight(self):
    
    self.translate(10,0)


  def moveUp(self):

    self.translate(0,-10)
    
     
  def moveDown(self):
  
    self.translate(0,10)
    
    
  def wheelEvent(self, event, callSister=True):
    
    if event.delta() < 0:
      z = 1.2
    else:
      z = 0.8333333333333333
    
    self.zoom(z)
    self._updatePosCoords( (event.x(), event.y()) )
    
    event.accept()
    
  
  def enterEvent(self, event):
    
    self.setFocus(Qt.OtherFocusReason)
    event.accept()


  def leaveEvent(self, event):
    
    self._updatePosCoords()
    event.accept()
  
  
  def _getClickChromoPos(self, x, y):  
    
    tickRegions, xRegions = self.genomeBrowserX._getChromoRegions()
    tickRegions, yRegions = self.genomeBrowserY._getChromoRegions()
 
    xSeq = ySeq = 0
    xChromo = yChromo = ''
 
    for chromo, pix1, nPix, seqStart, seqEnd, j in xRegions:
      pix2 = pix1 + nPix
 
      if pix1 <= x < pix2:
        f = (x-pix1) / float(pix2-pix1)
        w = seqEnd-seqStart
        xSeq = seqStart + int(f * w)
        xChromo = chromo
        break
 
    for chromo, pix1, nPix, seqStart, seqEnd, j in yRegions:
      pix2 = pix1 + nPix
 
      if pix1 <= y < pix2:
        f = (y-pix1) / float(pix2-pix1)
        w = seqEnd-seqStart
        ySeq = seqStart + int(f * w)
        yChromo = chromo
        break
    
    return xChromo, xSeq, yChromo, ySeq
    
    
  def mousePressEvent(self, event):
    
    QtGui.QWidget.mousePressEvent(self, event)
  
    pos = event.pos()
    buttons = event.buttons()
    mods = event.modifiers()
    have_ctrl = mods & Qt.CTRL
    self.prevPos = QPoint(event.pos())
  
    if buttons & LeftButton:
      x, y = event.x(), event.y()
      xChromo, xSeq, yChromo, ySeq = self._getClickChromoPos(x, y)
      print('Map X {}:{:,} Y {}:{:,}'.format(xChromo, xSeq, yChromo, ySeq))
  
    elif buttons & RightButton:

      x, y = event.x(), event.y()
      seq_pos = self._getClickChromoPos(x, y)
      self.select_start = seq_pos
  
  def keyPressEvent(self, event):
    
    QtGui.QWidget.keyPressEvent(self, event)
    nuc = self.mainApp.nuc
    
    x, y = self.prevPos.x(), self.prevPos.y()
    chromo, xSeq, yChromo, ySeq = self._getClickChromoPos(x, y)
    
    mods = event.modifiers()
    have_ctrl = mods & Qt.CTRL
   
    if have_ctrl and event.key() == Qt.Key_Z:
      code = 'Recording'
      typ = 'derived'
    
      if code in nuc.dataTracks[typ]:
        group = nuc.dataTracks[typ][code]
        chromo_group = nuc._getGroup(chromo, group)

        if 'regions' in chromo_group:
          regions = array(chromo_group['regions'], int32) # start, end
          values  = array(chromo_group['values']) # origValue, normValue
 
          regions = regions[:-1]
          values = values[:-1]
        
          if len(regions):
            nuc._setData('regions', chromo_group, uint32, regions)
            nuc._setData('values', chromo_group, float32, values)
          
          else:
            nuc._delete('regions', chromo_group)
            nuc._delete('values', chromo_group)
            nuc._delete(chromo, group)
            
          self.draw() 
          
   
  def mouseReleaseEvent(self, event):
   
    QtGui.QWidget.mouseMoveEvent(self, event)
    nuc = self.mainApp.nuc
  
    if self.select_start:
      x, y = event.x(), event.y()
      seq_pos = self._getClickChromoPos(x, y)
      
      xChromo1, xSeq1, yChromo1, ySeq1 = self.select_start
      xChromo2, xSeq2, yChromo2, ySeq2 = seq_pos
      
      print('{}:{:,}  {}:{:,} - {}:{:,}  {}:{:,}'.format(xChromo1, xSeq1, yChromo1, ySeq1, xChromo2, xSeq2, yChromo2, ySeq2))
      
      if len(set([xChromo1, yChromo1, xChromo2, yChromo2])) == 1:
        
        chromo = xChromo1
        code = 'Recording'
        typ = 'derived'
        
        if code in nuc.dataTracks[typ]:
          group = nuc.dataTracks[typ][code]
          chromo_group = nuc._getGroup(chromo, group)
 
          a = (xSeq1 + ySeq1)/2
          b = (xSeq2 + ySeq2)/2
 
          a = int(round(a, -4))
          b = int(round(b, -4))
          
          if 'regions' in chromo_group:
            regions = array(chromo_group['regions'], int32) # start, end
            values  = array(chromo_group['values']) # origValue, normValue
 
            regions = concatenate([regions, [[a,b],]])
            values  = concatenate([values, [[1.0,1.0],]])
 
          else:
            regions  = array([[a,b],], int32)
            values   = array([[1.0,1.0],])
           
          nuc._setData('regions', chromo_group, uint32, regions)
          nuc._setData('values', chromo_group, float32, values)
 
          self.draw()
   
    self.select_start = None
    
    
  def mouseMoveEvent(self, event):
    
    QtGui.QWidget.mouseMoveEvent(self, event)
    
    x = event.x()
    y = event.y()
    dx = x - self.prevPos.x()
    dy = y - self.prevPos.y()
    buttons = event.buttons()
    pos = event.pos()
    
    if (buttons & LeftButton) or (buttons & MiddleButton):
      mods = event.modifiers()
      haveCtrl = mods & Qt.CTRL
      haveShift = mods & Qt.SHIFT
      
      if self.flipY:
        dy *= -1
      
      if haveCtrl or haveShift:
        self.translate(-dx, -dy)
      else:
        self.translate(dx, dy)
      
    elif buttons & RightButton:
      pass
    
    self.prevPos = QPoint(pos)
    self._updatePosCoords( (x,y) )
  
  
  def _updatePosCoords(self, pos=None):
  
    return
  
    if pos is None:
      pass
  
    x, y = pos
  
    xRegions = self.genomeBrowserX.chromoRegions
    yRegions = self.genomeBrowserY.chromoRegions
    
    xSeq = ySeq = xWidth = yWidth = 0 
    xChromo = yChromo = ''
    
    for chromo, pix1, nPix, seqStart, seqEnd, j in xRegions:
      pix2 = pix1 + nPix
      
      if pix1 <= x < pix2:
        f = (x-pix1) / float(pix2-pix1)
        w = seqEnd-seqStart
        xSeq = seqStart + int(f * w)
        xChromo = chromo
        xWidth += w
        break
    
    for chromo, pix1, nPix, seqStart, seqEnd, j in yRegions:
      pix2 = pix1 + nPix
      
      if pix1 <= y < pix2:
        f = (y-pix1) / float(pix2-pix1)
        w = seqEnd-seqStart
        ySeq = seqStart + int(f * w)
        yChromo = chromo
        yWidth += w
        break
 
    #self.regionCallback(xChromo, xSeq, yChromo, ySeq)

  
  def updateMatrix(self, edgeColor=(48, 48, 64)):
    
    nuc = self.mainApp.nuc
    
    if nuc:
      contactGroups = nuc.getSelectedContactGroups()
      chromosomes = nuc.getDisplayedChromosomes()
      
      nC = len(contactGroups)
      if not nC:
        self.matrix = self.blank
        return
      
      gamma = 10 ** self.gammaWidget.get()
      
      t0 = time()
      tN = 0.0
 
      seqLen = 0
      chromoSizes = {}
      sortChromos = []
      chromoStarts = {}
      chromoSeqStarts = {}
      cisRegions = {}
      chromo_idx = {}
 
      for i, chromo in enumerate(chromosomes):
        start, end = nuc.getChromosomeLimits(chromo)
        
        size = (end-start)
        if not size:
          continue
 
        sortChromos.append(chromo)
        chromoSizes[chromo] = size
        chromoStarts[chromo] = seqLen
        chromoSeqStarts[chromo] = start
        seqLen += size
        chromo_idx[chromo] = i
        
      if not seqLen:      
        self.matrix = self.blank
        return
      
      x1, x2, y1, y2 = self.viewRegion
      dx = x2-x1
      dy = y2-y1
      
      nBinsY = self.height()
      binSize = (dy * seqLen) / nBinsY

      nBinsX = min(self.width(), int(seqLen/binSize))
      
      #x2 = min(1.0, x1 + (dx * nBinsX) / nBinsY )
      x2 = min(1.0, x1 + (nBinsX*binSize/float(seqLen)))      
      
      nBinsAllX = nBinsX / dx
      nBinsAllY = nBinsY / dy
 
      matrix    = [zeros((nBinsY, nBinsX), int32) for group in contactGroups]
      
      startBinX = x1 * nBinsAllX
      endBinX   = startBinX + nBinsX
      startBinY = y1 * nBinsAllY
      endBinY   = startBinY + nBinsY
      
      for chromo in sortChromos:
        chromoSizes[chromo] /= binSize
        chromoStarts[chromo] /= binSize
             
      cache = {}
      pX = 0.0
      pY = 0.0
      divsX = []
      divsY = []
      
      if self.flipY:
        chromosA = sortChromos[::-1]
      else:
        chromosA = sortChromos
      
      for chrA in chromosA:
        startA = chromoStarts[chrA]
        if startA > endBinY:
          continue
 
        sizeA = chromoSizes[chrA]
        endA = startA + sizeA
        if endA < startBinY:
          continue
 
        wholeA = True
 
        if startA < startBinY:
          pY1 = startBinY-startA
          wholeA = False
        else:
          pY1 = 0.0
 
        if endA > endBinY:
          pY2 = sizeA - (endA-endBinY)
          wholeA = False
        else:
          pY2 = sizeA
 
        seqStartA = chromoSeqStarts[chrA] + pY1 * binSize # chromoSeqStarts[chrA]
        deltaY = pY2-pY1
        divsX = []
        pX = 0.0
        for chrB in sortChromos:
          startB = chromoStarts[chrB]
          if startB > endBinX:
            continue
 
          sizeB = chromoSizes[chrB]
          endB = startB + sizeB
          if endB < startBinX:
            continue
 
          wholeB = True
 
          if startB < startBinX:
            wholeB = False
            pX1 = startBinX-startB
          else:
            pX1 = 0.0
 
          if endB > endBinX:
            wholeB = False
            pX2 = sizeB - (endB-endBinX)
          else:
            pX2 = sizeB
 
          seqStartB = chromoSeqStarts[chrB] + pX1 * binSize # chromoSeqStarts[chrB]
          deltaX = pX2-pX1
 
          if pX < nBinsX:
            divsX.append(pX)
              
          chromoPair = frozenset([chrA, chrB])
          pixX = min(nBinsX, int(deltaX))
          pixY = min(nBinsY, int(deltaY))
          
          if not (pixX * pixY):
            pX += deltaX
            continue  
          
          if wholeA and wholeB and (chrA != chrB) and (chromoPair in cache):
            subMat = [(s, b, a, x.T) for s, a, b, x in cache[chromoPair]]
            del cache[chromoPair]
          
          else:
            
            subMat = []
            for g, groupName in enumerate(contactGroups):
              if self.splitMap:
                split_side = nuc.getContactsDiagonalSide(groupName)
              else:
                split_side = 0
                
              sumData = zeros((pixY, pixX), int32)
              contDict = nuc.getCachedContacts(groupName)
 
              if chrA == chrB:
                if (chrA in contDict) and (chrB in contDict[chrA]):
                  data = array(contDict[chrA][chrB], int32)
 
                  if wholeA and wholeB:
                    binContacts(data, sumData, seqStartA, seqStartB, binSize, 1, 0, split_side)
 
                  else:
                    binContacts(data, sumData, seqStartA, seqStartB, binSize, 0, 0, split_side)
                    binContacts(data, sumData, seqStartB, seqStartA, binSize, 0, 1, split_side)
 
              else:
                if (chrA in contDict) and (chrB in contDict[chrA]):
                  data = array(contDict[chrA][chrB], int32)
                  binContacts(data, sumData, seqStartA, seqStartB, binSize, 0, 0)
 
                if (chrB in contDict) and (chrA in contDict[chrB]):
                  data = array(contDict[chrB][chrA], int32)
                  binContacts(data, sumData, seqStartB, seqStartA, binSize, 0, 1)
 
              subMat.append((split_side, chrA, chrB, sumData))
            
            if (chrA != chrB) and wholeA and wholeB:
              cache[chromoPair] = subMat
            
          if not len(subMat):
            pX += deltaX
            continue  

          sx1 = int(pX)
          sy1 = int(pY)
          pX1 = int(pX1)
          pY1 = int(pY1)
          
          if self.flipY:
            for i, (split_side, chrA, chrB, data) in enumerate(subMat):
              if (split_side < 0) and (chromo_idx[chrA] > chromo_idx[chrB]):
                continue
 
              elif (split_side > 0) and (chromo_idx[chrA] < chromo_idx[chrB]):
                continue
            
              matrix[i][sy1:sy1+pixY, sx1:sx1+pixX] = data[::-1,:]
          
          else: 
            for i, (split_side, chrA, chrB, data) in enumerate(subMat):
              if (split_side < 0) and (chromo_idx[chrA] > chromo_idx[chrB]):
                continue
 
              elif (split_side > 0) and (chromo_idx[chrA] < chromo_idx[chrB]):
                continue
              
              matrix[i][sy1:sy1+pixY, sx1:sx1+pixX] = data
           
          if chrA == chrB:
            # seqStarts, pixStarts, pixEnds, binSize
            cisRegions[chrA] = (seqStartA, seqStartB,
                                sy1, sx1, sy1+pixY, sx1+pixX,
                                binSize)
 
          pX += deltaX
 
        if pY < nBinsY:
          divsY.append(pY)
 
        pY += deltaY
      
      pixmap = zeros((int(pY)+1, int(pX)+1, 3), int32)
      
      # Main contact matrix
      
      for i, groupName in enumerate(contactGroups):
        m = matrix[i]
        r, g, b, a = nuc.getContactGroupColor(groupName).rgba32
        
        if nuc.areContactsSingleCell(groupName):
          addContactMapPoints(pixmap, m, int32(r), int32(g), int32(b))
        
        else:
          binSizeD = nuc.getContactsBinSize(groupName)
          nPix = int32(ceil(binSizeD/binSize))
          addContactMapBlocks(pixmap, m, int32(r), int32(g), int32(b), nPix, int32(m.max()), gamma)
      
      # Add any data track boxes
      
      if cisRegions:
        for typ, code in nuc.getDisplayedDataTracks():
          if nuc.getDataTrackPeakType(typ, code) != 2:
            continue # Displayed in genome browser instead
 
          pixRegions, pixValues = self._getDataTrackRegions(nuc, cisRegions, typ, code)
          r, g, b = nuc.getDataTrackColor(typ, code).rgb24
 
          if pixRegions is not None:
            #print  typ, code, cisRegions.keys(), pixRegions.shape, pixValues.shape, r, g, b
            addContactMapRegions(pixmap, pixRegions, pixValues, int32(r), int32(g), int32(b))
      
      # Chromosome edge grid lines
      
      r, g, b = [int32(x) for x in edgeColor]
      addContactMapMesh(pixmap, array(divsX, int32), array(divsY, int32), r, g, b)
      pixmap = array(pixmap, uint8)
      
      #print "Time taken: %.3f" % (time()-t0)
      
    else:
      self.matrix = self.blank
      h, w = self.blank.shape[:2]   
      self.genomeBrowserX.setRegion(0.0, 1.0, w)
      self.genomeBrowserY.setRegion(0.0, 1.0, h)
      return
    
    h, w = pixmap.shape[:2]
    self.genomeBrowserX.setRegion(x1, x2, w)
    self.genomeBrowserY.setRegion(y1, y2, h, self.flipY)
    
    self.viewRegion = (x1, x2, y1, y2)
    self.matrix = pixmap
    
    
  def draw(self, clearCache=True):
    
    self.makeCurrent()
    self.updateMatrix()
    sy, sx, sz = self.matrix.shape
       
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGB,
                    sx, sy, 0, GL.GL_BGR, GL.GL_UNSIGNED_BYTE,
		    self.matrix)

    GL.glNewList(self.glList, GL.GL_COMPILE)
    GL.glBegin(GL.GL_QUADS)
    
    GL.glTexCoord2f(0.0, 0.0)
    GL.glVertex2f(0.0, 0.0)
    
    GL.glTexCoord2f(0.0, 1.0) 
    GL.glVertex2f(0.0, sy)
    
    GL.glTexCoord2f(1.0, 1.0)
    GL.glVertex2f(sx, sy)
    
    GL.glTexCoord2f(1.0, 0)
    GL.glVertex2f(sx, 0)
    GL.glEnd()
    GL.glEndList()
   
    if clearCache:
      self.parent().pixmapCache = {}
  
    self.update()
  
  def _getDataTrackRegions(self, nuc, chromoRegions, typ, code):
    
    pixRegions = None
    pixValues  = None
    minValue, maxValue = nuc.getDataTrackThresholds(typ, code)
    
    for chromo in chromoRegions:
      regions = nuc.getDataTrackRegions(code, chromo, typ, minValue=minValue) # Only data regions above threshold
      
      if regions is None:
        continue # No track data for this chromosome
      
      values = nuc.getDataTrackValues(code, chromo, typ, minValue=minValue)
      seqStartA, seqStartB,  pixStartA, pixStartB, pixEndA, pixEndB, binSize = chromoRegions[chromo] # Display region
      seqEndA = seqStartA + (pixEndA-pixStartA) * binSize
      seqEndB = seqStartB + (pixEndB-pixStartB) * binSize
      
      
      regs, vals = getContactMapRegions(regions, array(values[:,1], float),
                                        seqStartA, seqStartB,
                                        pixStartA, pixStartB,
                                        pixEndA, pixEndB, binSize)
      
      if pixRegions is not None:
        pixRegions = concatenate([pixRegions, regs], axis=0)
        pixValues = concatenate([pixValues, vals], axis=0)
      
      else:
        pixRegions = regs
        pixValues = vals
      
    return pixRegions, pixValues
