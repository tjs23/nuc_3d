from colorsys import hsv_to_rgb
from numpy import abs, array, float32, empty, arange, log10, searchsorted, ones, full
from numpy import zeros, int32, vstack, eye, sqrt, clip, dot, linalg, concatenate, cos, sin
from OpenGL import GL, GLU
from OpenGL import GLUT # Temporary for postional labels
from cUtil import OpenGlUtil, drawing, contour
from cUtil import dataLayer, apiUtil   

from gui.Gl3dPanel import Gl3dPanel
from gui.qtgui.ToolBar import ToolBar
from gui.qtgui.Label import Label
from gui.qtgui.Layout import FlowLayout
from gui.qtgui.Button import Button
from gui.qtgui.Entry import IntRangesEntry, IntEntry, Entry
from gui.qtgui.Frame import FlowFrame
from gui.qtgui.Menu import Menu
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.Slider import FloatSlider, Slider
from gui.qtgui.SpinBox import FloatSpinBox, IntSpinBox
from gui.qtgui.Colors import ColorDialog, GradientEditor
from gui.qtgui.Base import Icon

from NucApi import DATA_TRACK_SYMBOL_PARAMS
from NucApi import DISPLAY_MODES, RESTRAINT_COLOR_MODES

from PySide import QtCore, QtGui
from time import time

# Add fog option, - maybe disable on mouse move
# Large chromosome letters at centres - float in upper layer

TAU = 2.0 * 3.14159265358979323846

# Temporary for postional labels
GLUT_BITMAP_HELVETICA_12 = GLUT.GLUT_BITMAP_HELVETICA_12
GLUT_BITMAP_HELVETICA_18 = GLUT.GLUT_BITMAP_HELVETICA_18
glColor4f = GL.glColor4f
glRasterPos3d = GL.glRasterPos3d
glutBitmapCharacter = GLUT.glutBitmapCharacter

"""
Struc panel:
  Sparse points icon
    Select which set
    
GL panel:
  Separately scale and rotate sphere
  Mouse Shift+Ctrl
    
"""

"""
# Example sphere instancing
# IMplement with shaders

sphVertexBuf = GL.glGenBuffers(1)
GL.glBindBuffer(GL.GL_ARRAY_BUFFER, sphVertexBuf)
GL.glBufferData(GL.GL_ARRAY_BUFFER, 4*ico.size, ico.ravel(), GL.GL_STATIC_DRAW)

# The VBO containing the centre positions of spheres
sphPosBuff = GL.glGenBuffers(1)
GL.glBindBuffer(GL.GL_ARRAY_BUFFER, sphPosBuff)
GL.glBufferData(GL.GL_ARRAY_BUFFER, 4*coords.size, coords.ravel(), GL.GL_DYNAMIC_DRAW)

# The VBO containing the colors of the spheres
sphColorBuf = GL.glGenBuffers(1)
GL.glBindBuffer(GL.GL_ARRAY_BUFFER, sphColorBuf)
GL.glBufferData(GL.GL_ARRAY_BUFFER, 4*colorArray.size, colorArray.ravel(), GL.GL_DYNAMIC_DRAW)

# Attribute buffer : vertices
GL.glEnableVertexAttribArray(0)
GL.glBindBuffer(GL.GL_ARRAY_BUFFER, sphVertexBuf)
GL.glVertexAttribPointer(0, 3, GL.GL_FLOAT,  GL.GL_FALSE, 0, 0) # Attrib, pSize, dType, normed, stride, offset

# Attribute buffer : particle centers
GL.glEnableVertexAttribArray(1)
GL.glBindBuffer(GL.GL_ARRAY_BUFFER, sphPosBuff)
GL.glVertexAttribPointer(1, 3, GL.GL_FLOAT, GL.GL_FALSE, 0, 0)

# Attribute buffer : particles colors
GL.glEnableVertexAttribArray(2)
GL.glBindBuffer(GL.GL_ARRAY_BUFFER, sphColorBuf)
GL.glVertexAttribPointer(2, 4, GL.GL_FLOAT, GL.GL_FALSE, 0, 0) # Could use normalised unsigned bytes instead

# Set attribut index advancement
GL.glVertexAttribDivisor(0, 0) # Mesh reused, no index advance
GL.glVertexAttribDivisor(1, 1) # Posotion (center) indices advance 1
GL.glVertexAttribDivisor(2, 1) # Color indices advance 1

# Draw mesh at instanced positions
GL.glDrawArraysInstanced(GL.GL_TRIANGLES, 0, len(ico), len(coords))
""" 

def _calcCoordMesh(coords, inColors, radius, nVoxels):

  return apiUtil.calcCoordMesh(coords, inColors, radius, nVoxels)


class StructureOuterPanel(QtGui.QWidget):
  
  def __init__(self, parent, mainApp, openFunc, **kw):
    
    QtGui.QWidget.__init__(self, parent)
  
    self.mainApp = mainApp
    self.openFunc = openFunc
    self.getIcon = mainApp.getIcon
    
    layout = QtGui.QGridLayout(self)
    layout.setSpacing(2)
    layout.setContentsMargins(2,2,2,2)
    layout.setRowStretch(2, 1)
    layout.setColumnStretch(0, 1)
    self.setLayout( layout )

    # Structure toolbar
    
    icons = ['zoom-in.png', 'zoom-out.png',
             'zoom-center.png', 'move-center.png', #'move-left.png', 'move-right.png', 'move-up.png', 'move-down.png',
             'display-line.png','display-tube.png',
             'display-ball.png','solid.png','display-restraints.png',
             'color-seq.png', 'color-chromo.png',
             'color-data.png', #'color-faint.png',
             'color-model.png',
             'rmsd.png',]
    
    icons = [self.getIcon(i) for i in icons]
    
    texts = ['Zoom in','Zoom out','Reset zoom','Center', #'Move left','Move right','Move up','Move down'
             'Line display', 'Tube display',
             'Ball and stick display', 'Surface display',
             'Toggle restraints display',
             'Sequence postion colors', 'Chromosome colours',
             'Data track colours', # 'Faint colouring',
             'Different model colours',
             'Structure RMSD colours']
    
    funcs = [self.zoomIn, self.zoomOut,
             self.zoomReset, self.moveCenter, #self.moveLeft, self.moveRight, self.moveUp, self.moveDown,
             self.displayLine, self.displayTube, 
             self.displayBallStick, self.displayMesh,
             self.toggleRestraints,
             self.colorSeq, self.colorChromo,
             self.colorGenData, # self.colorFaint,
             self.colorModel, self.colorRmsd]

    self.strucToolbar = ToolBar(self, 'Structure toolbar', funcs, icons, texts, 
                                objName='strucToolbar', areas='tblr', iconSize=32)
    
    # Model selection
    
    self.modelButton = Button(self.strucToolbar, ' ', iconSize=32,
                              icon=self.getIcon('model.png'),
                              tipText='Set coordinate model')
    
    self.modelMenu = Menu(self.modelButton, 'Configure rendering options', 
                          setupFunc=self._setupModelMenu)
                     
    self.modelButton.setMenu(self.modelMenu)
    self.strucToolbar.addWidget(self.modelButton)
    

    #clipping in space
    clipButton = Button(self.strucToolbar, '', icon=self.getIcon('clip.png'),
                        tipText='Clip front and back planes', iconSize=32)
    clipMenu = Menu(clipButton, 'Structure clipping')
    clipButton.setMenu(clipMenu)
    
    self.clipSliderFrontH = FloatSlider(self.strucToolbar, startVal=-1.0, endVal=1.0, value=1.0,
                              direction='h', step=0.01, bigStep=0.1, callback=self.setClipF,
                              tracking=True, showNumber=False, tickInterval=None,
                              tickPosition=None, listener=None, decimals=2)
                              
    self.clipSliderBackH = FloatSlider(self.strucToolbar, startVal=1.0, endVal=-1.0, value=1.0,
                              direction='h', step=0.01, bigStep=0.1, callback=self.setClipB,
                              tracking=True, showNumber=False, tickInterval=None,
                              tickPosition=None, listener=None, decimals=2)

    self.clipSliderSphere = FloatSlider(self.strucToolbar, startVal=0.0, endVal=1.0, value=1.0,
                              direction='h', step=0.005, bigStep=0.01, callback=self.setClipS,
                              tracking=True, showNumber=False, tickInterval=None,
                              tickPosition=None, listener=None, decimals=3)
                              
    self.clipSliderSphere.setMinimumWidth(200)

    self.clipSliderHBackAction = clipMenu.addItem('Back:', widget=self.clipSliderBackH)
    self.clipSliderHFrontAction = clipMenu.addItem('Front:', widget=self.clipSliderFrontH)
    self.clipSliderSphereAction = clipMenu.addItem('Sphere:', widget=self.clipSliderSphere)
    clipMenu.addItem('Reset clipping', callback=self.resetClip, icon=self.getIcon('edit-undo.png'))
    
    self.strucToolbar.addWidget(clipButton)

    
    #expanding in space
    expandButton = Button(self.strucToolbar, '', icon=self.getIcon('expand.png'),
                        tipText='Expand structure', iconSize=32)
    expandMenu = Menu(expandButton, 'Expand')
    expandButton.setMenu(expandMenu)
    
    self.expansionSliderH = FloatSlider(self.strucToolbar, startVal=1.0, endVal=5.0, value=1.0,
                              direction='h', step=0.1, bigStep=0.5, callback=self.setExpansion,
                              tracking=True, showNumber=True, tickInterval=None,
                              tickPosition=None, listener=None, decimals=1)

    self.expansionSliderHAction = expandMenu.addItem('Expansion:', widget=self.expansionSliderH)
    
    self.strucToolbar.addWidget(expandButton)
    
    self.restActionB = self.strucToolbar.getActions()[7]
    self.restActionB.setCheckable(True)
    self.restActionB.setChecked(True)     
    
    renderOptButton = Button(self.strucToolbar, '  ', iconSize=32,
                            icon=self.getIcon('configure.png'),
                            tipText='Configure rendering options')
    renderOptMenu = Menu(renderOptButton, 'Configure rendering options', 
                         setupFunc=self._setupDispDetailsMenu)
    renderOptButton.setMenu(renderOptMenu)
    self.strucToolbar.addWidget(renderOptButton)

    # Colour toolbar
    
    colorConfigButton = Button(self.strucToolbar, '  ', icon=self.getIcon('colors.png'),
                               tipText='Colour settings', iconSize=32)
    colorConfigMenu = Menu(colorConfigButton, 'Colours configure',
                           setupFunc=self._setupColorConfigMenu)
    colorConfigButton.setMenu(colorConfigMenu)
    self.strucToolbar.addWidget(colorConfigButton)
  
    layout.addWidget(self.strucToolbar, 0, 0, 1, 1)

    
    # Structure selection toolbar
    
    self._button_structures = {}
    self.strucSelectButtons = []
    bg = self.structureButtonGroup = QtGui.QButtonGroup(self)
    bg.setExclusive(True)
    bg.connect(bg, QtCore.SIGNAL('buttonClicked(int)'),  self.toggleStructure)
    
    #tb = ToolBar(self, 'Structure selection toolbar', [], None, None, 
    #             objName='strucSelectFrame', areas='tblr',
    #             isVertical=False, iconSize=(32,32))
 
    self.strucSelectFrame = FlowFrame(self)
    #tb.addWidget(self.strucSelectFrame)
 
    layout.addWidget(self.strucSelectFrame, 1, 0)
    
    
    
    # Main 3D widget
                           
    self.innerPanel = StructurePanel(self, mainApp, openFunc, **kw)
    
    layout.addWidget(self.innerPanel, 2, 0)
      
    
  def _updateModelButton(self):
  
    nuc = self.mainApp.nuc
    
    if nuc:
      models = nuc.getDisplayedModels()
      nModels = len(models)
      
      if nModels  == 1:
        text = ' %s' % (1+models[0],)
        
      elif nModels > 1:
        text = ' *'
        
      else:
        text = ' ' 
      
      
    else:
      text = ' '
  
    self.modelButton.setText(text)
  
  
  def _setDisplayModels(self, models):
    
    if self.mainApp.nuc.setDisplayedModels(models):
      self.updateContents()
    
    
  def _setDisplayModelsRange(self, models):
    
    models = [x-1 for x in models]
    
    
    if self.mainApp.nuc.setDisplayedModels(models):
      self.modelMenu.setActiveAction(None)
      if self.modelMenu.isVisible():
        self.modelMenu.hide()
      
      self.updateContents()

    
  def resetClip(self):
  
    self.clipSliderFrontH.set(1.0, doCallback=False)
    self.clipSliderBackH.set(1.0, doCallback=False)
    self.clipSliderSphere.set(1.0, doCallback=False)
    
    self.innerPanel.clipPlaneFront = 1.0
    self.innerPanel.clipPlaneBack = 1.0
    self.innerPanel.clip_sphere = 1.0
      
    self.mainApp.nuc.notifiers.add('/display/attrs/models')
    self.updateContents()
 
  
  def setClipS(self, radius):
  
    self.clipSliderSphere.set(radius, doCallback=False)
    self.innerPanel.clip_sphere = radius 
       
    self.mainApp.nuc.notifiers.add('/display/attrs/models')
    self.updateContents()
    
      
  def setClipF(self, planePos):
  
    self.setClip(planePos, True)
  
  
  def setClipB(self, planePos):
  
    self.setClip(planePos, False)
    

  def setClip(self, zPos, isFront):
    
    pfh = self.clipSliderFrontH.get()
    pbh = self.clipSliderBackH.get()
    
    if isFront:   
      if pbh < -zPos:
        zPos = -pbh
              
      self.clipSliderFrontH.set(zPos, doCallback=False)
      self.innerPanel.clipPlaneFront = zPos
      
    else:

      if pfh < -zPos:
        zPos = -pfh
         
      self.clipSliderBackH.set(zPos, doCallback=False)
      self.innerPanel.clipPlaneBack = zPos
      
    self.innerPanel.update()


  def setExpansion(self, expansion):
    
    #set values
    self.expansionSliderH.set(expansion, doCallback=False)
    self.innerPanel.expansion = expansion

    #add notifier and redraw
    self.mainApp.nuc._setAttr(self.mainApp.nuc.structure,'expansion',expansion)
    #self.mainApp.nuc.notifiers.add('/structures/attr/expansion')
    self.innerPanel.update()
    self.refreshStructure()
  
  
  def refreshStructure(self):

    self.innerPanel.construct(self.mainApp.nuc)
    self.innerPanel.paintEvent(None)
 
 
  def zoomIn(self):

    self.innerPanel.zoom(0.66666666667)
    
    
  def zoomOut(self):
  
    self.innerPanel.zoom(1.5)
    
    
  def zoomReset(self):
  
    self.innerPanel.resetZoom()
  
  
  def moveCenter(self):
  
    self.innerPanel.resetView()
    
    
  def moveLeft(self):
  
    self.innerPanel.translate(25, 0)
    
    
  def moveRight(self):
  
    self.innerPanel.translate(-25, 0) 
    
    
  def moveUp(self):
  
    self.innerPanel.translate(0, 25)
    
    
  def moveDown(self):
  
    self.innerPanel.translate(0, -25)


  def _setDisplayMode(self, value):
    
    refresh = False
    
    nuc = self.mainApp.nuc
    for chromo in nuc.getChromosomes():
      if nuc.setChromoDisplayParams(chromo, displayMode=value):
        refresh = True
    
    if refresh:
      self.updateContents()


  def _setColorMode(self, value):
    
    refresh = False
    
    nuc = self.mainApp.nuc
    for chromo in nuc.getChromosomes():
      if nuc.setChromoDisplayParams(chromo, colorMode=value):
        refresh = True
    
    if refresh:
      self.updateContents()
    
    return refresh
    
 
  def colorSeq(self):

    self._setColorMode(0)
 
 
  def colorChromo(self):

    self._setColorMode(1)
  
 
  def colorGenData(self):

    self._setColorMode(2)
 
 
  def colorModel(self):

    self._setColorMode(3)
 
 
  def colorRmsd(self):

    self._setColorMode(4)


  #def colorFaint(self):
  #
  #  self._setColorMode(6)
 

  def displayBallStick(self):

    self._setDisplayMode(0)
 
 
  def displayLine(self):

    self._setDisplayMode(1)


  def displayTube(self):

    self._setDisplayMode(2)
    

  def displayMesh(self):

    self._setDisplayMode(3)


  def toggleRestraints(self, *args):
    
    if self.mainApp.nuc:
      attrs = self.mainApp.nuc.display.attrs
      opts = list(attrs['options'])      
      showCis    = opts[2]
      showTrans  = opts[3]   
      
      if showCis or showTrans:
        self.mainApp.nuc.setRestraintsDisplayed(False, False)
        self.mainApp.restActionA.setChecked(False)
        self.restActionB.setChecked(False)
      
      else:
        self.mainApp.nuc.setRestraintsDisplayed(True, True)
        self.mainApp.restActionA.setChecked(True)
        self.restActionB.setChecked(True)
       
      self.updateContents()
  
  
  def updateContents(self):
  
    nuc = self.mainApp.nuc
    
    if len(nuc.structures.keys()) > 1:
      self.strucSelectFrame.show()
      self.updateStrucSelectToolbar()
    else:
      self.strucSelectFrame.hide()

    nuc.notifiers.add('/display/attrs/models')

    self._updateModelButton()
    
    attrs = self.mainApp.nuc.display.attrs
    
    showCis    = attrs['options'][2]
    showTrans  = attrs['options'][3]   
    
    showRest = bool(showCis or showTrans)
    self.restActionB.setChecked(showRest)     
    self.mainApp.restActionA.setChecked(showRest)     
    
    alpha = nuc.display.attrs['chromoAlpha']
    
    self.innerPanel.alpha = alpha
    self.innerPanel.construct(nuc)
    self.innerPanel.update()


  def _alphaUpdate(self, alpha):
    
    nuc = self.mainApp.nuc
    nuc.display.attrs['chromoAlpha'] = alpha
    nuc.notifiers.add('/display/attrs/chromoAlpha')
 
    self.innerPanel.alpha = alpha
    self.innerPanel.construct(nuc)
    self.innerPanel.update()
    
  
  def _setupColorConfigMenu(self, menu):
  
    menu.clear()
    nuc = self.mainApp.nuc
      
    if nuc:
      r, g, b, a = self.innerPanel.bgColor
      icon = Icon(None, '#%02X%02X%02X' % (255*r,255*g,255*b), 12)
      menu.addItem('Structure background', icon=icon, callback=self.setStrucBgColor)
    
      attrs = nuc.display.attrs
      restColorMode = attrs['options'][6]
      restColorMenu = Menu(menu, 'Restraint colouring')
      
      for i, name in enumerate(RESTRAINT_COLOR_MODES):
        restColorMenu.addItem(name, checked=bool(i == restColorMode), object=i,
                              callback=lambda j:self.setRestraintColoring(j))
      
      alphaMenu = Menu(menu, 'Chromosome opacity')
      self.alphaWidget =  FloatSlider(alphaMenu, startVal=0.1, endVal=1.0,
                                      value=self.innerPanel.alpha,
                                      direction='h', step=0.05, bigStep=None,
                                      callback=self._alphaUpdate)                      
      alphaMenu.addItem(u'\u03B1', widget=self.alphaWidget)
      
      menu.addSeparator()
                              
      schemeAttrs = nuc.display['colorSchemes'].attrs
      names = sorted([n for n in schemeAttrs])
      
      for name in names:
        colors = ['#%02X%02X%02X' % (255*r,255*g,255*b) for r,g,b in schemeAttrs[name]]
        subMenu = Menu(menu, name.title() + ' scheme', icon=None)
        
        callback=lambda c=colors, n=name:self.setGradientColors(c, n)
        widget = GradientEditor(None, colors=colors, callback=callback)
        subMenu.addItem('', widget=widget)
      
      menu.addSeparator()


  def setStrucBgColor(self):
  
    prev = self.innerPanel.bgColor
    colorObj = ColorDialog(self).getColor(prev)
            
    if colorObj:
      rgb =  list(colorObj.getRgbF())
      self.innerPanel.bgColor = rgb
      self.mainApp.interactomePanel.bgColor = rgb
      self.updateContents()
  
  
  def setRestraintColoring(self, mode):
    
    nuc = self.mainApp.nuc
    
    if nuc:
      changed = nuc.setRestraintColoring(mode)
      
      if changed:
        self.updateContents()
  
    
  def setGradientColors(self, colors, name):
    
    nuc = self.mainApp.nuc
    
    rgbs = [[int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)] for c in colors]
    rgbs = (array(rgbs) / 255.0).tolist()
    
    group = nuc.display['colorSchemes']
    schemeAttrs = nuc._setAttr(group, name, rgbs)
    
    self.updateContents()    


  def recalcDensity(self):
    
    nuc = self.mainApp.nuc
    
    if nuc:
      attrs = nuc.display.attrs
      densityRad =  attrs['sizes'][4]                       
      chromos = nuc.getDisplayedChromosomes()
      
      nuc.calcDensity(0, chromos, densityRad)
  
      for chromo in chromos:
        chrColorMode = nuc.getChromoDisplayParams(chromo)[2]

        if chrColorMode == 2:
          self.updateContents()
          break

  
  def toggleStructure(self, toggled):
    
    nuc = self.mainApp.nuc
    
    if nuc:
      buttons = self.structureButtonGroup.buttons()
      selection = [self._button_structures[b] for b in buttons if b.isChecked()]
      
      if selection:
        structCode = selection[0]
 
        if structCode != nuc.structure.name.split('/')[-1]:
          nuc.setCurrentStructure(structCode)
          self.updateContents()


  def updateStrucSelectToolbar(self):
        
    nuc = self.mainApp.nuc
    
    self._button_structures = {}
    
    if nuc:
      structGroup = nuc.structures
      codes = [x for x in structGroup.keys() if nuc.getNumModels(x)]
      
      c = len(codes)
      
      if c < 2:
        self.strucSelectFrame.hide()
        return
      elif self.innerPanel.isVisible():
        self.strucSelectFrame.show()
              
      buttons = self.structureButtonGroup.buttons()
      
      try:
        codes.sort(key=lambda a:'%5d' % int(a))
      except ValueError:
        codes.sort()
      
      for i, code in enumerate(codes):
        attrs = structGroup[code].attrs
        
        #isSelected = bool(attrs['isActive'])
        isSelected = structGroup[code] == nuc.structure
        
        text = attrs['name']
        tipText = 'Hide/show structure %s' % text
        
        if i < len(self.strucSelectButtons):
          button = buttons[i]
          button.setChecked(isSelected)
          button.setText(text)
          #button.setIcon(None)
          button.setToolTip(tipText)
          self.strucSelectButtons[i].setVisible(True)
        
        else:
          button = Button(self.strucSelectFrame, text, icon=None,
                          tipText=tipText, grid=None)
          button.setCheckable(True)
          button.setChecked(isSelected)
 
          self.structureButtonGroup.addButton(button)
          self.structureButtonGroup.setId(button, i)
          self.strucSelectButtons.append(button)
        
        self._button_structures[button] = code
      
      a = len(self.strucSelectButtons)
      
      if a > c:
        for i in range(c, a):
          self.strucSelectButtons[i].setVisible(False)

  def _setupModelMenu(self, menu):
  
    menu.clear()
    nuc = self.mainApp.nuc
    
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
    nuc = self.mainApp.nuc
    
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
                                  step=0.5, callback=self.mainApp.setSphereSize,
                                  multiplier=1.1)

      stickSizeEntry = FloatSpinBox(menu, stickSize, minValue=0.0, maxValue=1.0,
                                    step=0.25, callback=self.mainApp.setStickSize)

      sphDetailEntry = IntSpinBox(menu, sphDetail, minValue=0, maxValue=4, step=1,
                                  callback=self.mainApp.setSphereDetail)
 
      lineWidthEntry = FloatSpinBox(menu, lineWidth, minValue=0.5, maxValue=5.0,
                                    step=0.5, callback=self.mainApp.setLineWidth)

      lineSmoothEntry = IntSpinBox(menu, lineSmooth, minValue=0, maxValue=4, step=1,
                                   callback=self.mainApp.setLineSmooth)

      tubeSizeEntry = FloatSpinBox(menu, tubeSize, minValue=0.001, maxValue=10.0,
                                   step=0.5, callback=self.mainApp.setTubeSize, multiplier=1.1)

      tubeSmoothEntry = IntSpinBox(menu, tubeSmooth, minValue=0, maxValue=4, step=1,
                                   callback=self.mainApp.setTubeSmooth)

      tubeDetailEntry = IntSpinBox(menu, tubeDetail, minValue=0, maxValue=4, step=1,
                                   callback=self.mainApp.setTubeDetail)

      menu.addItem('Cis restraints', checked=bool(showCis),
                   callback=lambda a:self.mainApp.setRestraintLinesCis(a.isChecked()))
      menu.addItem('Trans restraints', checked=bool(showTrans), 
                   callback=lambda a:self.mainApp.setRestraintLinesTrans(a.isChecked()))
      
                   
      menu.addItem('Chromosome labels', checked=bool(showChromoLabels), 
                   callback=lambda a:self.mainApp.setChromoLabels(a.isChecked()))
      menu.addItem('Position labels', checked=bool(showLabels), 
                   callback=lambda a:self.mainApp.setTextLabels(a.isChecked()))
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



class StructurePanel(Gl3dPanel):

  def __init__(self, parent, mainApp, openFunc, **kw):
    
    Gl3dPanel.__init__(self, parent, mainApp, openFunc)
    self.rGlList = 1
    self.iGlList = 2
    self.cGlList = 5
    self.glWidget = self
    self.mainApp = mainApp
    self.alpha = 1.0
    
    self.chromoGlLists = {}
    self.chromoCenters = {}
    self.chromoGenomePos = {}
    self.genomeSize = 0
    self.restraintGlList = None
    self.interactGlList = None
    self.origin = zeros(3)
    self.setAutoFillBackground(False)
    self.setAttribute(QtCore.Qt.WA_NoSystemBackground)       
    self.fontS = QtGui.QFont('Helvetica', 10)
    self.fontM = QtGui.QFont('Helvetica', 14)
    self.fontL = QtGui.QFont('Helvetica', 20, 1)
    self.clip_sphere = 1.0
    
    self.image_scale = 1.0
    self.image_rotate = eye(3)
    self.image_translate = array([0.0, 0.0, 0.0])
      

  def mouseReleaseEvent(self, event):
    
    selectIndex = Gl3dPanel.mouseReleaseEvent(self, event)
    
    if selectIndex:
      p = round(self.genomeSize * (selectIndex/4228250625.0))
 
      for chromo in self.chromoGenomePos:
        start, span = self.chromoGenomePos[chromo]
 
        if start <= p < (start+span):
          bp = int(p)-start
          print('Pick point: Chr%s:%d' % (chromo, bp))
          break
    
    
  #draw the chromosome labels on each chromosome
  def updateQtLayer(self, painter):
    
    nuc = self.mainApp.nuc
    
    if not nuc:
      return
    
    chrs = nuc.getDisplayedChromosomes(nuc.structure.name)
    nChrs = len(chrs)
    #size of canvas
    w0, h0 = self.glGeometry[:2]
        
    if nChrs == 1:
      chromoText = 'Chromosome : %s' % (chrs[0])
      
    elif nChrs:
      chromoText = 'Chromosomes : %s' % (','.join(chrs))
    
    else:
      return
    
    dispAttrs = nuc.display.attrs
    showText = dispAttrs['options'][4]
    
    if not showText:
      return
    
    painter.setRenderHint(painter.Antialiasing)
    
    # TBD HUD ...
    #painter.setFont(self.fontM)    
    #painter.setPen(QtGui.QColor(255, 255, 128, 128))
    #painter.drawText(QtCore.QPoint(10, 20), chromoText)
    
    painter.setFont(self.fontL)
    fm = painter.fontMetrics()
    
    r0, g0, b0 = self.bgColor[:3]
    fg_color = QtGui.QColor( 255-int(255*r0), 255-int(255*g0), 255-int(255*b0), 255)
    bg_color = QtGui.QColor( int(255*r0), int(255*g0), int(255*b0), 192)
    
    p = 5.0

    #find the mean of the chromosome centres in the z direction
    #set the centre of each chromosome
    modelCentres = {}
    t = self.modelMatrix    
    zMean = 0.0
    for chrA in self.chromoCenters:
      x, y, z = self.chromoCenters[chrA]           
      cen = dot([x,y,z,1.0], t)
      zMean += cen[2]
      modelCentres[chrA] = cen

    zMean /= float(len(self.chromoCenters))
    
    #draw each chromosome label
    for chrA in self.chromoCenters:
      x, y, z, w = modelCentres[chrA]
      
      #if z < zMean:
      #  continue
      
      r, g, b, a = nuc.getChromoColor(chrA).rgba
      
      # prespective 
      x *= -h0/z
      y *= h0/z
      x += w0/2.0
      y += h0/2.0   
      
      bbox = fm.tightBoundingRect(chrA)
      bw = float(bbox.width())
      bh = float(bbox.height())
      
      color = QtGui.QColor(int(255*r), int(255*g), int(255*b), 255)
      painter.setPen(color)
      painter.setBrush(bg_color)
      painter.drawRoundedRect(x-p, y-bh/2.0-p, bw+p+p, bh+p+p, p, p)
      
      painter.setPen(fg_color)
      painter.drawText(x, y+bh/2.0, chrA)
  
  
  def construct(self, nuc=None, models=None):   
    """
    Setup OpenGL vertex, colour and normal data - using GlLists for now...
    """
    
    self.makeCurrent()
    GL.glEnable(GL.GL_DEPTH_TEST)
    self.lighting = set()
    self.pickBuffer = set()
    glLists = []
    glListsPerChromo = 5 # Data track, labels, chromo coordinates, ROI, position pick buffer
    
    if not nuc:
      return Gl3dPanel.construct(self)
    
    coordsGroup = nuc._getCoordsGroup()
    if not coordsGroup:
      return Gl3dPanel.construct(self)
       
    particGroup = nuc._getParticleGroup()
    chromosomes = nuc.getChromosomes(nuc.structure.name)
    
    zAxis = array([0.0, 0.0, 1.0])   
    nModels = nuc.getNumModels()
    models = range(nModels) 
    pickColors = {} # Chromosome particle colours for screen position pick buffer
    alpha = self.alpha
    
    # Get display parameters
    params = self._getDisplayParams(nuc)
    
    colorMode, dispMode, showCis, showTrans, showChromoText, showText, \
    restColorMode, radius, stickFrac, lineWidth, tubeWidth, showModels, \
    sphDetail, lineSmooth, tubeSmooth, tubeDetail, chromoParams = params
    
    models = [m for m in models if m in showModels]
    
    # Determine updates - what will be redrawn
    
    #some actions have to be done on all chromosomes, some on restraints too.
    doChromos = set()
    doRestraints = False
    doInteractions = False
    for gName in nuc.notifiers:
      if gName.startswith('/structures/') and ('/restraints' in gName):
        doRestraints = True
      
      elif gName.startswith('/chromosomes/'):
        if gName.endswith('/attrs/display'):
          doRestraints = True
          c = gName.split('/')[2]
          doChromos.add(c) 
          
        if gName.endswith('/attrs/region_of_interest'):
          doRestraints = True
          c = gName.split('/')[2]
          doChromos.add(c) 

      elif gName.startswith('/structures/') and ('/coords' in gName):
        doRestraints = True
        parts = gName.split('/')
        if len(parts) > 4:
          doChromos.add(parts[4]) 
                 
      elif gName.startswith('/dataTrack/'):
        if gName.startswith('/dataTrack/derived/density/'):
          for c in chromosomes:
            chrShow, chrText, chrColorMode, chrDispMode = nuc.getChromoDisplayParams(c)
            if chrShow and chrColorMode == 2:
              doChromos.add(c)
           
        else:
          for c in chromosomes:
            chrShow, chrText, chrColorMode, chrDispMode = nuc.getChromoDisplayParams(c)
            if chrShow and chrColorMode == 4:
              doChromos.add(c)
      
      elif gName.startswith('/interactions/'):
        doInteractions = True

      elif gName.startswith('/display/attrs/detailLevel'):
        doChromos.update(chromosomes)
        # could redraw only if the current mode changed

      elif gName.startswith('/display/attrs/sizes'):
        doChromos.update(chromosomes)

      elif gName.startswith('/display/attrs/chromoAlpha'):
        doRestraints = True
        doChromos.update(chromosomes)

      elif gName.startswith('/display/attrs/interactome'):
        doRestraints = True
        doChromos.update(chromosomes)

      elif gName.startswith('/display/colorSchemes'):
        doChromos.update(chromosomes)
        doRestraints = True
      
      elif gName.startswith('/display/attrs/options'):
        doRestraints = True
        doInteractions = True

      elif gName.startswith('/display/attrs/models'):
        doRestraints = True
        doChromos.update(chromosomes)

      elif gName.startswith('/structures/') and gName.endswith('/attrs/expansion'):
        doRestraints = True
        doChromos.update(chromosomes)
          
    # Get regions of interest      
    chromo_rois = {}
    roi = None
    
    for chromo in chromosomes:
      rois = nuc.getChromoRegionsOfInterest(chromo)
      if rois is not None:
        selected = rois[:,4].nonzero()[0]
        rois = rois[selected]
        
        if len(rois):
          alpha = 0.75
        else:
          rois = []
          
      else:
        rois = []
      
      chromo_rois[chromo] = rois    
      
      if len(rois) and not roi:
        roi = (chromo, rois[0][0], rois[0][1])
      
    # Centre and cache all model coords
    # for each chromosome, its centre is the mean of all model centres
    args = (nuc, coordsGroup, particGroup, models, self.clip_sphere, roi)
    n, coordCache, self.origin, self.clipLimit, self.chromoCenters, idx_segments, seq_segments, pos_cache = self._getChromoCoords(*args)
    
    # Check for blank
    
    if self.reference and not self.mainApp.coordImage:
      for i in self.reference:
        GL.glDeleteLists(i, 1)
      
      self.reference = set()
      
    if not n:
      if self.restraintGlList:
        GL.glDeleteLists(self.rGlList, 1)

      if self.interactGlList:
        GL.glDeleteLists(self.iGlList, 1)
        
      for chromo in self.chromoGlLists:
        glList = self.chromoGlLists[chromo]
        GL.glDeleteLists(glList, glListsPerChromo)
      
      return # nothing to draw - use a default?
    
    
    # Draw restraints
         
    if showCis or showTrans:
      if doRestraints or not self.restraintGlList:
        if self.restraintGlList:
          GL.glDeleteLists(self.rGlList, 1)
      
        self.restraintGlList = self.rGlList
        GL.glNewList(self.rGlList, GL.GL_COMPILE)
        self._drawRestraints(nuc, coordCache, models, showCis,
                             showTrans, restColorMode, alpha, idx_segments, lineWidth)
        GL.glEndList()
      
      if self.restraintGlList:
        glLists.append(self.rGlList)
    
    # Draw interactions
    
    iCodes = nuc.getDisplayedInteractions()
    if iCodes:
      if doInteractions or doRestraints or not self.interactGlList:
        if self.interactGlList:
          GL.glDeleteLists(self.iGlList, 1)
      
        self.interactGlList = self.iGlList
        GL.glNewList(self.iGlList, GL.GL_COMPILE)
        self._drawInteractions(nuc, particGroup, coordCache, iCodes, models, alpha, seq_segments)
        GL.glEndList()
      
      if self.interactGlList:
        glLists.append(self.iGlList)

    
    # Get genomic position of chromosomes  for pick colours
    genomeSize = 0
    genomePos = {}
    for chromo in coordCache:
      start, end = nuc.getChromosomeLimits(chromo)
      span = end - start     
      genomePos[chromo] = (genomeSize, span) 
      genomeSize += span
      
    self.genomeSize = genomeSize  
    self.chromoGenomePos = genomePos
        
    # Draw tracks and labels, setup colours
    
    drawChromos = []
    for chromo in coordCache:
      chrShow, chrText, chrColorMode, chrDispMode = chromoParams[chromo]
      
      if not chrShow:
        continue
           
      skip = False
      if chromo in self.chromoGlLists:
        glList = self.chromoGlLists[chromo]
        if chromo in doChromos:
          #delete the layers that will be redrawn
          GL.glDeleteLists(glList, glListsPerChromo)
        else:
          skip = True
        
      elif self.chromoGlLists:
        glList = max(self.chromoGlLists.values()) + glListsPerChromo
        self.chromoGlLists[chromo] = glList
      
      else:
        glList = self.cGlList
        self.chromoGlLists[chromo] = glList
          
      #get linear bead position data for drawing the chromosome
      chrColor = list(nuc.getChromoColor(chromo).rgba)
      pos = pos_cache[chromo]
      nCoords = len(pos)
        
      #
      # DATA TRACK DRAWING
      #
      
      if chrColorMode != 5: # Not faint
        dataSymbols = []
        self._addDataSymbols(nuc, chromo, dataSymbols)
 
        if dataSymbols:
          GL.glNewList(glList, GL.GL_COMPILE)
          self._drawDataTrack(nuc, chromo, pos, coordCache, models,
                              dataSymbols, lineWidth, lineSmooth, seq_segments[chromo])
          GL.glEndList()
          glLists.append(glList)
      
      #
      # POSITIONAL LABELS
      #
      
      glList += 1
      if showText:
        labelData = nuc.getChromoLabels(chromo)
        
        if labelData is not None:
          GL.glNewList(glList, GL.GL_COMPILE)
          self._drawPosLabels(nuc, chromo, labelData, pos, coordCache,
                              models, lineWidth, lineSmooth, seq_segments[chromo])
          GL.glEndList()
          glLists.append(glList)               
      
      #
     
      if skip:
        for i in range(glListsPerChromo):
          glLists.append(glList+i)
        
        if chrDispMode in (0,2):
          self.lighting.add(glList+2)
          
        continue
 
      #
      # CHROMOSOME COLOURS
      #

      colorArray = zeros((nCoords, 4), float)
      faint = array([0.4, 0.4, 0.4, 0.4 * alpha])
      
      if chrColorMode == 0: # Sequence
        scheme = nuc.getColorScheme('sequence')
        colorArray = nuc.calcInterpolatedColors(scheme, nCoords, alpha)
     
      elif colorMode == 1: # Chromo
        chrColor[3] = alpha
        colorArray[:] = chrColor
 
      elif chrColorMode == 2: # Data track
        dataColors = self._getDataColors(nuc, chromo, pos, alpha, None)
        colorArray = dataColors
        
      elif chrColorMode == 3: # Model
        scheme = nuc.getColorScheme('model')
        schemeColors = nuc.calcInterpolatedColors(scheme, nModels, alpha)

      elif chrColorMode == 4: # RMSD
        scheme = nuc.getColorScheme('rmsd')
        schemeColors = nuc.calcInterpolatedColors(scheme, 256, alpha)
        colorArray = self._getRmsdColors(nuc, chromo, schemeColors)
 
      elif chrColorMode == 5: # Faint
        colorArray[:] = faint
        
      else: # Default to Chromo
        chrColor[3] = alpha
        colorArray[:] = chrColor
        
      for start, end, color_opt, thickness, selected, color_r, color_g, color_b in chromo_rois[chromo]:
        if color_opt:
          p1, p2 = searchsorted(pos, [start, end])
 
          if color_opt == 1:
            roi_color = chrColor

          elif color_opt == 3:
            scheme = nuc.getColorScheme('density')
            roi_color = nuc.calcInterpolatedColors(scheme, p2-p1, alpha)
            
          else:
            roi_color = array([color_r/255.0, color_g/255.0, color_b/255.0, 1.0])
 
          colorArray[p1:p2] = roi_color
        
      drawChromos.append( (chromo, colorArray, glList) )
      glList += 2
    
    # Parallel pre-calculations
    
    #t0 = time()
    meshDict = {}
    meshCalc = []
    meshColors = []
    
    for chromo, colorArray, glList in drawChromos:
      chrShow, chrText, chrColorMode, chrDispMode = chromoParams[chromo]

      if chrShow and chrDispMode == 3: # Surface mesh
        meshCalc.append(chromo)
        meshColors.append(colorArray)
 
    if meshCalc:
      paraEngine = nuc._getParallelEngine()
      radius, nVoxels = 2.0, 75
      
      if paraEngine is None:
        for i, chromo in enumerate(meshCalc):
          coords = coordCache[chromo][0]
          meshDict[chromo] = apiUtil.calcCoordMesh(coords, meshColors[i], radius, nVoxels)
        
      else:
        coords = [coordCache[chromo][0] for chromo in meshCalc]
        sameArgs={'radius':radius, 'nVoxels':nVoxels, }
        diffArgs={'coords':coords, 'inColors':meshColors}
        job = paraEngine.run(_calcCoordMesh, sameArgs, diffArgs)
 
        for i, data in enumerate(job.getResult()):
          chromo = meshCalc[i]
          meshDict[chromo] = data
    
    #
    # Draw models
    #
    
    vertex_cache = {}
    
    for chromo, colorArray, glList in drawChromos:
              
      glList += 1
      GL.glNewList(glList, GL.GL_COMPILE)
      glLists.append(glList)   
      GL.glEnable(GL.GL_NORMALIZE)
      lighting = False
      vList = None
      
      rois = chromo_rois[chromo]
      chrShow, chrText, chrColorMode, chrDispMode = chromoParams[chromo]
      vertex_cache[chromo] = {}
      
      for model in coordCache[chromo]:
        vertex_cache[chromo][model] = []
        
        # Draw coordinates
        
        for a, b in idx_segments[chromo]:
          if b-a < 2:
            continue
        
          coords = coordCache[chromo][model][a:b+1]
          segment_colors = colorArray[a:b+1]
          
          if len(coords) < 2:
            continue
 
          if chrColorMode == 3:
            segment_colors[:,:] = schemeColors[model]
 
          #
          # CHROMOSOME DRAWING
          #
 
          if chrDispMode == 1: # Lines
 
            # Interpolate coords for smoother path
            if lineSmooth:
              for i in range(lineSmooth):
                coords = drawing.smoothPath(coords)
 
              colors = apiUtil.getInterpolatedColors(segment_colors, len(coords))

            else:
              colors = segment_colors
            
            if chrColorMode == 5:# Faint
              GL.glLineWidth(1.0)
            else:
              GL.glLineWidth(lineWidth)
 
            GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
            GL.glEnableClientState(GL.GL_COLOR_ARRAY)
 
            # GL vertices/colours must be flat array
            cList = colors.ravel()
            vList = coords.ravel()
 
            # Cache vertex coords for pick layer
            vertex_cache[chromo][model].append(vList)
 
            GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
            GL.glColorPointer(4, GL.GL_FLOAT, 0, cList)
 
            # draw here
            GL.glDrawArrays(GL.GL_LINE_STRIP, 0, len(vList)//3)
 
            GL.glDisableClientState(GL.GL_COLOR_ARRAY)
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
 
          else:  # Ball+stick, tubes or surface mesh
            lighting = True
 
            if chrDispMode == 0: # Ball+stick
 
              radius2 = radius * stickFrac
              vList, cList, nList = drawing.createBallAndStick(coords, segment_colors, radius, radius2, sphDetail)
              drawType = GL.GL_TRIANGLES

            elif chrDispMode == 2: # Tubes
 
              if tubeSmooth:
                for i in range(tubeSmooth):
                  coords = drawing.smoothPath(coords)
 
                colors = apiUtil.getInterpolatedColors(segment_colors, len(coords))
 
              else:
                colors = segment_colors
 
              vList, cList, nList  = drawing.createTube(coords, colors, tubeWidth, 4 + tubeDetail*4)
              drawType = GL.GL_TRIANGLE_STRIP

            elif chrDispMode == 3: # Surface mesh
 
              if meshDict:
                vList, cList, nList = meshDict[chromo]
              else:
                vList, cList, nList = apiUtil.calcCoordMesh(coords, array(segment_colors), radius=2.0, nVoxels=75)
 
              #print "Mesh calc time:%.3f" % (time()-t0)
 
              drawType = GL.GL_TRIANGLES
 
            # Cache vertex coords for pick layer
            vertex_cache[chromo][model].append(vList)
 
            GL.glLineWidth(1.0)
            GL.glDisable(GL.GL_LINE_SMOOTH)
 
            GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
            GL.glEnableClientState(GL.GL_COLOR_ARRAY)
            GL.glEnableClientState(GL.GL_NORMAL_ARRAY)
 
            GL.glNormalPointer(GL.GL_FLOAT, 0, nList)
            GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
            GL.glColorPointer(4, GL.GL_FLOAT, 0, cList)
 
            GL.glDrawArrays(drawType, 0, len(vList)/3)
 
            GL.glDisableClientState(GL.GL_NORMAL_ARRAY)
            GL.glDisableClientState(GL.GL_COLOR_ARRAY)
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
 
            GL.glEnable(GL.GL_LINE_SMOOTH)
        
        # draw all of ROI, irrespective of clipped segments etc
      
      GL.glEndList()
      
      if vList is None:
        continue

      if lighting:
        self.lighting.add(glList)
      
      glList += 1
      GL.glNewList(glList, GL.GL_COMPILE)
      glLists.append(glList)
      
        
      for model in coordCache[chromo]:
         
        # Draw ROIs
 
        for start, end, color_opt, thickness, selected, color_r, color_g, color_b in rois:
          
          #print start, end, color_opt, thickness, selected, color_r, color_g, color_b
            
          p1, p2 = searchsorted(pos_cache[chromo], [start, end])
          
          if p1 == p2:
            continue
          
          if p1 > p2:
            p1, p2, = p2, p1

          coords = coordCache[chromo][model][p1:p2]
          segment_colors = colorArray[p1:p2]
          segment_colors[:,3] = 1.0 # No transparency
          
          #print p1, p2, coords, chrDispMode
          
          if chrDispMode == 9:
            lighting = False
            
            # Interpolate coords for smoother path
            if lineSmooth:
              for i in range(lineSmooth):
                coords = drawing.smoothPath(coords)
 
              colors = apiUtil.getInterpolatedColors(segment_colors, len(coords))

            else:
              colors = segment_colors
 
            GL.glLineWidth(lineWidth * 3)
 
            GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
            GL.glEnableClientState(GL.GL_COLOR_ARRAY)
 
            # GL vertices/colours must be flat array
            cList = colors.ravel()
            vList = coords.ravel()
 
            GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
            GL.glColorPointer(4, GL.GL_FLOAT, 0, cList)
 
            GL.glDrawArrays(GL.GL_LINE_STRIP, 0, len(vList)//3)
 
            GL.glDisableClientState(GL.GL_COLOR_ARRAY)
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
          
          elif chrDispMode in (0,1,2):
            lighting = True
            
            if chrDispMode == 0:# Ball+stick
 
              radius2 = radius * stickFrac
              vList, cList, nList = drawing.createBallAndStick(coords, segment_colors, radius, radius2, sphDetail)
              drawType = GL.GL_TRIANGLES

            elif chrDispMode in (1,2): # Tubes
 
              if tubeSmooth:
                for i in range(tubeSmooth):
                  coords = drawing.smoothPath(coords)
 
                colors = apiUtil.getInterpolatedColors(segment_colors, len(coords))
 
              else:
                colors = segment_colors
 
              vList, cList, nList  = drawing.createTube(coords, colors, tubeWidth, 4 + tubeDetail*4)
              drawType = GL.GL_TRIANGLE_STRIP
 
            GL.glLineWidth(1.0)
            GL.glDisable(GL.GL_LINE_SMOOTH)
 
            GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
            GL.glEnableClientState(GL.GL_COLOR_ARRAY)
            GL.glEnableClientState(GL.GL_NORMAL_ARRAY)
 
            GL.glNormalPointer(GL.GL_FLOAT, 0, nList)
            GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
            GL.glColorPointer(4, GL.GL_FLOAT, 0, cList)
 
            GL.glDrawArrays(drawType, 0, len(vList)/3)
 
            GL.glDisableClientState(GL.GL_NORMAL_ARRAY)
            GL.glDisableClientState(GL.GL_COLOR_ARRAY)
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
 
            GL.glEnable(GL.GL_LINE_SMOOTH)

      if lighting:
        self.lighting.add(glList)
      
      GL.glDisable(GL.GL_NORMALIZE)
      GL.glEndList()
                         
      del coordCache[chromo] # free some memory
      
      # Simplified point picking with positional indexing RGBA colours
      # Drawn into into invisible picking buffer
      
      # Each chromo pos has a genome fraction
      # convert each fraction to a pick colour
      start, span = genomePos[chromo]
      fracs = (pos + start)/float(genomeSize)
      pickColors = drawing.getFractionPickColors(fracs)
      
      # Interpolate to current number of vertices, assume even expansion
      # Positional colours same for all models
      pickColors = apiUtil.getInterpolatedColors(pickColors, len(vList)//3, blend=False)
      pickColors = pickColors.ravel()
      
      glList += 1
      GL.glNewList(glList, GL.GL_COMPILE)
      glLists.append(glList)
      self.pickBuffer.add(glList)
      
      for model in vertex_cache[chromo]:
        
        for vList in vertex_cache[chromo][model]: # Already smoothed vertices etc for eacg chromo segment     
        
          if len(vList) < 2:
            continue

          if chrDispMode == 1: # Lines
            GL.glLineWidth(lineWidth)
            drawType = GL.GL_LINE_STRIP
          else: # Ball+stick or tubes
            GL.glLineWidth(1.0)
            drawType = GL.GL_TRIANGLE_STRIP
 
          GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
          GL.glEnableClientState(GL.GL_COLOR_ARRAY)
 
          GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
          GL.glColorPointer(4, GL.GL_FLOAT, 0, pickColors)
 
          GL.glDrawArrays(drawType, 0, len(vList)//3)
 
          GL.glDisableClientState(GL.GL_COLOR_ARRAY)
          GL.glDisableClientState(GL.GL_VERTEX_ARRAY)

      GL.glEndList()
    
    #
    #  Image point clouds
    #
    
    if self.mainApp.coordImage:
      coords = nuc.getImageCoords(self.mainApp.coordImage)
      
      #cen = coords.mean(axis=0)
      
      coords -= coords.mean(axis=0)
      #coords *= 10.0/coords.max()
      coords /= 1e3
      
      #coords -= cen
      coords *= self.image_scale
      coords = dot(coords, self.image_rotate)
      
      coords += self.image_translate
      
      radius = 0.25 * self.image_scale
      detail = 1
      alpha = 1.0
      scheme = array([[1.0, 0.8, 0.0],[1.0, 0.8, 0.0]])
      colorArray = nuc.calcInterpolatedColors(scheme, len(coords), alpha)
      
      vList, cList, nList = drawing.createBallAndStick(coords, colorArray, radius, 0, detail)
      #vList = coords.ravel()
      #cList = colorArray.ravel()
      
      drawType = GL.GL_TRIANGLES
      #drawType = GL.GL_POINTS
     
      glList += 1
      GL.glNewList(glList, GL.GL_COMPILE)
      glLists.append(glList)
      self.reference = set([glList])
      self.lighting.add(glList)
      
      GL.glPointSize(5.0)            
      GL.glEnable(GL.GL_NORMALIZE)
      GL.glDisable(GL.GL_LINE_SMOOTH)
      
      GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
      GL.glEnableClientState(GL.GL_COLOR_ARRAY)
      GL.glEnableClientState(GL.GL_NORMAL_ARRAY)
 
      GL.glNormalPointer(GL.GL_FLOAT, 0, nList)
      GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
      GL.glColorPointer(4, GL.GL_FLOAT, 0, cList)
        
      GL.glDrawArrays(drawType, 0, len(vList)/3)
      
      GL.glDisableClientState(GL.GL_NORMAL_ARRAY)
      GL.glDisableClientState(GL.GL_COLOR_ARRAY)
      GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
      
      GL.glEnable(GL.GL_LINE_SMOOTH)   
      GL.glDisable(GL.GL_NORMALIZE)
      
      GL.glEndList()
    
    else:
      self.image_coords = None
     
    self.glLists = glLists
    
    nuc.notifiers = set()
  
  
  def _getDisplayParams(self, nuc):
      
    dispAttrs = nuc.display.attrs
    
    colorMode, dispMode, showCis, showTrans, \
      showChromoText, showText, restColorMode = dispAttrs['options'][:7]
    
    radius, stickFrac, lineWidth, tubeWidth = dispAttrs['sizes'][:4]
    
    showModels = nuc.getDisplayedModels()
    
    sphDetail, lineSmooth, tubeSmooth, tubeDetail = dispAttrs['detailLevel'][:4]
  
    chromoParams = {}
    for chromo in nuc.chromosomes:
      chromoParams[chromo] = nuc.getChromoDisplayParams(chromo)
  
    params = [colorMode, dispMode, showCis, showTrans, showChromoText, showText, 
              restColorMode, radius, stickFrac, lineWidth, tubeWidth, showModels, 
              sphDetail, lineSmooth, tubeSmooth, tubeDetail, chromoParams]
    
    return params
    

  def _getChromoCoords(self, nuc, coordsGroup, particGroup, models,
                       clip_radius=None, region_of_interest=None):
      
    n = 0
    coordCache = {}
    pos_cache = {}
    maxCoords = []
    minCoords = []
    centre = array([0.0, 0.0, 0.0])
    chromoCenters = {}
    idx_segments = {}
    seq_segments = {}
    
    for chromo in coordsGroup:
 
      if nuc.getChromoDisplayParams(chromo)[0]: # visible
        chromoCenters[chromo] = array([0.0, 0.0, 0.0])
        availModels = len(coordsGroup[chromo])
        coordCache[chromo] = {}
        
        pos = array(particGroup[chromo]['positions'], int32)
        pos_cache[chromo] = pos
        nPos = len(pos)
 
        for model in models:
          if model < availModels:
            coords = array(coordsGroup[chromo][model], float)
            if len(coords) > nPos: # Bug in CSM import
              coords = coords[:nPos]
 
            c = coords.sum(axis=0)
            chromoCenters[chromo] += c/len(coords)
            centre += c
 
            maxCoords.append(coords.max(axis=0))
            minCoords.append(coords.min(axis=0))
            n += len(coords)
            
            coordCache[chromo][model] = coords
        chromoCenters[chromo] /= float(len(models))
        
        # default segments whole region
        idx_segments[chromo] = array([[0, nPos-1],], int32)
        seq_segments[chromo] = array([[pos[0], pos[-1]],], int32)
    
    if n > 0:
      maxCoords = array(maxCoords).max(axis=0)
      minCoords = array(minCoords).min(axis=0)
      clipLimit = (maxCoords - minCoords).max()
        
      if 0.0 < clip_radius < 1.0:
        clip_radius = clipLimit * clip_radius
        
        if region_of_interest:
          chromo, region_start, region_end = region_of_interest
          a, b = searchsorted(pos_cache[chromo], [region_start, region_end])
          origin = zeros(3, float)
          
          for model in coordCache[chromo]:
            region_coords = coordCache[chromo][model][a:b+1]
            origin += region_coords.mean(axis=0)
          
          origin /= float(len(models))
          
        else:
          origin = centre / float(n)
        
        clipLimit = clip_radius * 1.5
        
        # spherical clip
        
        for chromo in coordCache:
          coords = array([coordCache[chromo][m] for m in sorted(coordCache[chromo])])
          
          idx_regions = apiUtil.getRadialClipRegions(coords, origin, clip_radius)
          seq_regions = [[pos_cache[chromo][a], pos_cache[chromo][b]] for a, b in idx_regions]
           
          idx_segments[chromo] = idx_regions
          seq_segments[chromo] = array(seq_regions, int32)
            
      else:
        clipLimit *= 0.7071067811865476    
        origin = centre / float(n)
      
      #idx_regions = [[0,182], [302,1646]]
      #seq_regions = [[pos_cache['X'][a], pos_cache['X'][b]] for a, b in idx_regions]
      #idx_segments['X'] = array(idx_regions)
      #seq_segments['X'] = array(seq_regions, int32)
      
      if 'expansion' in nuc.structure.attrs:
        expansion = nuc.structure.attrs['expansion']
        
        if expansion > 1.0:
          for chromo in coordsGroup:
            
            if nuc.getChromoDisplayParams(chromo)[0]: # visible
              delta = (expansion - 1.0) * (chromoCenters[chromo] - origin)

              for model in coordCache[chromo]:
              
                coordCache[chromo][model][:] += delta

              chromoCenters[chromo] += delta
              
          maxCoords += (expansion - 1.0) * (maxCoords - origin)
          minCoords += (expansion - 1.0) * (minCoords - origin)
          clipLimit = (maxCoords - minCoords).max() * 0.7071067811865476

    else:
      origin = centre
      clipLimit = 512.0
 
    return n, coordCache, origin, clipLimit, chromoCenters, idx_segments, seq_segments, pos_cache


  def _getDepthColors(self, nuc, chromo, nCoords, colors):

    depths = nuc.getDataTrackValues('density', chromo, 'derived')                       
    
    if (depths is None) or (len(depths) != nCoords):
      attrs = nuc.display.attrs
      structure = nuc.structure.name
      densityRad =  attrs['sizes'][5]                   
      # TBD Parallelise  per chromo
      nuc.calcDensity(0, nuc.getDisplayedChromosomes(structure), densityRad, structure) 
      depths = nuc.getDataTrackValues('density', chromo, 'derived') 
            
    indices = array(255*depths[:,1], int) # Normalised vals -> 0..255 indices
    
    return colors[indices]
  
    
  def _getRmsdColors(self, nuc, chromo, colors):

    modelRmsds, atomRmsds = nuc.calcModelRmsds(None, chromosomes=[chromo,])           
    # TBD Parallelise  per chromo
    
    atomRmsds /= 4.0 * atomRmsds.mean() or 1.0
    atomRmsds = clip(atomRmsds, 0.0, 1.0)
    
    indices = array(255*atomRmsds, int) # Normalised vals -> 0..255 indices
    
    return colors[indices]
  
    
  def _getDataColors(self, nuc, chromo, pos, alpha=0.5, segments=None):
    
    bgColor = (array(self.bgColor) + array([0.5, 0.5, 0.5, 1.0])) / 2.0
    bgColor[3] = alpha
    
    nCoords = len(pos)
    dataColors = zeros((nCoords, 4), float)
    signals = zeros(nCoords, float)
    nAdded = 0.0
    
    if chromo[-1] in 'AB':
      chromoB = chromo[:-1]
    else:
      chromoB = None
    
    for typ in nuc.dataTracks:
      typGroup = nuc.dataTracks[typ]
      
      for code in typGroup:
        attrs = typGroup[code].attrs
        
        if 'display' not in attrs:
          continue
        
        if not attrs['options'][0]: # Not selected
          continue
        
        if attrs['options'][2]: # Uses symbol
          continue
        
        refGroup = nuc.getRefDataTrackGroup(typ, code)
        
        if chromo in refGroup:
          dataChromo = chromo
 
        elif chromoB and chromoB in refGroup:
          dataChromo = chromoB
 
        else:
          continue
        
        group   = refGroup[dataChromo]
        regions = group['regions']
        
        if not len(regions):
          continue
        
        color   = attrs['display'][:4]
        thrsh   = attrs['display'][5]
        
        regions = array(regions, int32)
        values  = array(group['values'], float)
        
        if segments is not None:
          idx = apiUtil.pairRegionsIntersection(regions, segments)
          regions = regions[idx]
          values = values[idx]
        
        dataLayer.addLayerColor(dataColors, signals, pos, regions,
                                values, color, thrsh)
        nAdded += 1.0
    
    if nAdded:
      signals = vstack([signals, signals, signals, signals,]).T
 
      dataColors = clip(dataColors, 0.0, 1.0) # Additive RGB
      
      dataColors *= signals
      dataColors += (1.0-signals) * bgColor
      dataColors[:,3] = alpha
      
    else:
      dataColors[:] = bgColor
    
    return  dataColors 
  
  
  def _addDataSymbols(self, nuc, chromo, dataSymbols):
            
    if chromo[-1] in 'AB':
      chromoB = chromo[:-1]
    else:
      chromoB = None
    
    for typ in nuc.dataTracks:
      typGroup = nuc.dataTracks[typ]
 
      for code in typGroup:
        attrs = typGroup[code].attrs
 
        if 'display' not in attrs:
          continue
 
        if not attrs['options'][0]: # Not selected
          continue
 
        if not attrs['options'][2]: # Uses color
          continue
        
        refGroup = nuc.getRefDataTrackGroup(typ, code)
               
        if chromo in refGroup:
          dataChromo = chromo
 
        elif chromoB and chromoB in refGroup:
          dataChromo = chromoB
 
        else:
          continue
 
        group = refGroup[dataChromo]
        dataSymbols.append((group, attrs))
  
  
  def _drawInteractions(self, nuc, particGroup, coordCache, iCodes, models, alpha, segments=None):
    
    bgColor = (array(self.bgColor) + array([0.5, 0.5, 0.5, 1.0])) / 2.0
    bgColor[3] = alpha
    
    chromos = set(coordCache.keys())
    iGroup = nuc.interactions
    m = models[0]
    
    for code in iCodes:
      group = iGroup[code]
      rgba = group.attrs['display'][:4]
      vmin = group.attrs['display'][5]
      vmax = group.attrs['display'][6] or 1.0
      lineWidth = group.attrs['display'][4]
      style = group.attrs['options'][2]

      refGroup = nuc.getRefInteractionsGroup(code)
      
      GL.glLineWidth(lineWidth)        
      
      for chrA in refGroup:
        if chrA not in chromos:
          continue
 
        groupA = refGroup[chrA]
        pos_a = array(particGroup[chrA]['positions'])
        
        for chrB in groupA:
          if chrB not in chromos:
            continue
          
          groupB = groupA[chrB]     
          pos_b = array(particGroup[chrB]['positions'])
          
          regions = array(groupB['regions'], int32)
          values  = array(groupB['values'])[:,1] # normed
          
          idx = (values >= vmin).nonzero() # clip
          
          regions = regions[idx]
          values = values[idx]
          
          values = clip(values, vmin, vmax) / vmax
          values = vstack([values, values, values, values,]).T
          
          if segments is not None: # Filter chromosome segments
            if not len(segments[chrA]):
              continue

            if not len(segments[chrB]):
              continue

            pos_pairs = regions[:,(0,2)]
            idx = apiUtil.pairsDualRegionIntersection(pos_pairs, segments[chrA], segments[chrB])
            posA = pos_pairs[idx,0]
            posB = pos_pairs[idx,1]
            values = values[idx]
         
          else:
            posA = regions[:,0]
            posB = regions[:,2]
           
          ni = len(posA)
          
          colors = values * rgba # color is scaled by value
          colors += (1.0-values) * bgColor
         
          iColors = empty((ni, 2, 4), float) 
          iColors[:,0] = colors
          iColors[:,1] = colors
          iColors = iColors.ravel()
     
          iCoords = empty((ni, 2, 3), float)
          
          for model in models:
            coordsA = nuc.getPositionCoords(model, posA, chrA) # defaults to displayed structure
            coordsB = nuc.getPositionCoords(model, posB, chrB)

            iCoords[:,0] = coordsA
            iCoords[:,1] = coordsB
 
            vertexList = iCoords.ravel()
 
            GL.glLineStipple(1, style)
            GL.glEnable(GL.GL_LINE_STIPPLE)
 
            GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
            GL.glEnableClientState(GL.GL_COLOR_ARRAY)
            GL.glVertexPointer(3, GL.GL_FLOAT, 0, vertexList)
            GL.glColorPointer(4, GL.GL_FLOAT, 0, iColors)
 
            #GL.glColor4f(1.0, 1.0, 1.0, 0.5)
            GL.glDrawArrays(GL.GL_LINES, 0, len(vertexList)/3)
 
            GL.glDisableClientState(GL.GL_COLOR_ARRAY)
            GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
            GL.glDisable(GL.GL_LINE_STIPPLE)
 
  
  def _drawRestraints(self, nuc, coordCache, models, cis, trans, colouring,
                      alpha, segments=None, lineWidth=2.0, nColors=10):
      
    chromosomes = list(coordCache.keys())
    restraintDict = nuc.getRestraints(chromosomes)
    
    if not restraintDict:
      restraintDict = {}
      contact_dict = nuc.getContacts(chromosomes=chromosomes)
      particle_size = nuc.getBackboneSpacing()
      
      for key in sorted(contact_dict):
        chr1, chr2 = key
        start1, end1 = nuc.getChromosomeLimits(chr1)
        start2, end2 = nuc.getChromosomeLimits(chr2)
        rois1 = nuc.getChromoRegionsOfInterest(chr1)
        rois2 = nuc.getChromoRegionsOfInterest(chr2)
        d, n = contact_dict[key].shape
        
        if rois1 is not None:
          rois1 = array([r[:2] for r in rois1 if r[4]], int32)
          
          if not len(rois1):
            rois1 = None

        
        if rois2 is not None:
          rois2 = array([r[:2] for r in rois2 if r[4]], int32)
          
          if not len(rois2):
            rois2 = None
          
        contacts = contact_dict[key]
        pos1 = contacts[0].astype(int32)
        pos2 = contacts[1].astype(int32)
        
        if (rois1 is not None) or (rois2 is not None):
          if rois1 is not None:
            selected1 = zeros(n, int)
            idx1 = apiUtil.pointRegionsIntersection(pos1, rois1)
            selected1[idx1] = 1
          else:
            selected1 = ones(n, int)
          
          if rois2 is not None:
            selected2 = zeros(n, int)
            idx2 = apiUtil.pointRegionsIntersection(pos2, rois2)
            selected2[idx2] = 1
          else:
            selected2 = ones(n, int)
 
          selected = (selected1 * selected2).nonzero()[0]
          
          pos1 = pos1[selected]
          pos2 = pos2[selected]
          n = len(pos1)
 
        pos1 = ((pos1 - start1)//particle_size).astype(int32) +1
        pos2 = ((pos2 - start2)//particle_size).astype(int32) +1
        uppers = full((n,1), 2.0)
        lowers = full((n,1), 1.0)
       
        restraintDict[key] = concatenate([pos1[:,None], pos2[:,None], lowers, lowers, lowers, uppers], axis=1)

      
    particGroup = nuc._getParticleGroup()
    
    scheme = nuc.getColorScheme('restraint')
    colorArray = nuc.calcInterpolatedColors(scheme, nColors, 0.8)
    
    GL.glLineWidth(lineWidth)            
    
    for chromoPair in restraintDict:
      chromoA, chromoB = chromoPair
      
      print chromoPair
      
      if not cis and (chromoA == chromoB):
        continue
      
      if not trans and (chromoA != chromoB):
        continue
       
      if not (nuc.getChromoDisplayed(chromoA) and nuc.getChromoDisplayed(chromoB)):
        continue
      
      restraints = restraintDict[chromoPair]
      nRestraints = len(restraints)
      
      if not nRestraints:
        continue
      
      rest_pairs = array(restraints[:,:2], int32)
      
      if segments is not None: # Filter chromosome segments
      
        if not len(segments[chromoA]):
          continue

        if not len(segments[chromoB]):
          continue
      
        idx = apiUtil.pairsDualRegionIntersection(rest_pairs, segments[chromoA], segments[chromoB])
        rest_pairs = rest_pairs[idx]
        uppers = array(restraints[idx,5], float)
     
      else:
        uppers = array(restraints[:,5], float)
      
      
      idxA = rest_pairs[:,0]
      idxB = rest_pairs[:,1]
      nRestraints = len(idxA)
      
      rCoords = empty((nRestraints, 2, 3), float)
      
       
      modelCoordsA = coordCache[chromoA]
      modelCoordsB = coordCache[chromoB]
      
      m = models[0]     
      
      # Setup colours
      
      if colouring == 2:
        chrPos = array(particGroup[chromoA]['positions'])
      
      rColors = zeros((nRestraints, 2, 4), float)
      
      if colouring == 0: # Cis/Trans colours
        if chromoA == chromoB:
          data = colorArray[0]
        else:
          data = colorArray[-1]
        
        data[3] = alpha
        rColors[:,0] = data
        rColors[:,1] = data
 
      elif colouring == 2: # Seq sep. colours
 
        if chromoA == chromoB:
          deltas = log10(clip( abs(chrPos[idxA]-chrPos[idxB]), 1e4, 1e8 )) - 4.0
          cIndices = array(deltas/4.0, int)
          data = colorArray[cIndices]
          data[:,3] = alpha
          
        else:
          data = colorArray[-1]
          data[3] = alpha

        rColors[:,0] = data
        rColors[:,1] = data
 
      elif colouring == 3: # Chromosome colours
        data1 = list(nuc.getChromoColor(chromoA).rgba)
        data2 = list(nuc.getChromoColor(chromoB).rgba)
        
        data1[3] = alpha
        data2[3] = alpha
        
        rColors[:,0] = data1
        rColors[:,1] = data2
      
      for model in models:
        coordsA = modelCoordsA[model]
        coordsB = modelCoordsB[model]
        
        rA = coordsA[idxA]
        rB = coordsB[idxB]
        
        if colouring == 1: # Dist colours
          deltas = rA-rB
          dists = sqrt((deltas*deltas).sum(axis=1)) - uppers
          dists = clip(dists/uppers, 0.0, 4.0)
          cIndices = array(dists*(nColors-1)/4.0, int)
          data = colorArray[cIndices]
          data[:,3] = alpha
          rColors[:,0] = data
          rColors[:,1] = data
        
        rCoords[:,0] = rA
        rCoords[:,1] = rB
        
        vertexList = rCoords.ravel()
        colourList = rColors.ravel()
        
        #GL.glLineStipple(1, 43690) # int('0b1010101010101010', 2)
        GL.glLineStipple(1, 52428) # int('0b1100110011001100', 2)
        GL.glEnable(GL.GL_LINE_STIPPLE)
        
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
        GL.glEnableClientState(GL.GL_COLOR_ARRAY)
        GL.glVertexPointer(3, GL.GL_FLOAT, 0, vertexList)
        GL.glColorPointer(4, GL.GL_FLOAT, 0, colourList)
        
        GL.glColor4f(1.0, 1.0, 1.0, 0.5)
        GL.glDrawArrays(GL.GL_LINES, 0, len(vertexList)/3)
 
        GL.glDisableClientState(GL.GL_COLOR_ARRAY)
        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
        GL.glDisable(GL.GL_LINE_STIPPLE)

  
  def _drawPosLabels(self, nuc, chromo, labelData, pos, coordCache, models, lineWidth, lineSmooth, segments=None):
  
    # text, pos, color[0], color[1], color[2], opt[0], opt[1], opt[2]
    
    if segments is not None:
      if not len(segments): # emty selection
        return
        
      pos_b = array([x[1]+1 for x in labelData], int32)
      idx = apiUtil.pointRegionsIntersection(pos_b, segments)
      labelData = [labelData[i] for i in idx]
      
    if not labelData:
      return
  
    GL.glLineWidth(lineWidth)
    
    values = array([[1.0,1.0],], float)  
    dAngle = 3.0 * TAU/7.0
    nPoints = int32(7)
    vList = []
    cList = []
    annos = []
    
    pMin = pos[0]
    pMax = pos[-1]
      
    for model in models:
      coords = coordCache[chromo][model]
      
      if len(coords) < 2:
        continue
      
      for i in range(lineSmooth):
        coords = drawing.smoothPath(coords)
      
      for text, seqPos, r, g, b, o1, o2, o3 in labelData:
        if not (pMin <= seqPos <= pMax):
          continue
        
        regions = array([[seqPos, seqPos+1],], int32)
             
        symbolCoords = dataLayer.getSymbolCoords(pos, coords,  regions, values, dAngle,
                                                 nPoints, int32(lineSmooth), 
                                                 0.4, 0.0, minRadius=0.1)
        vList.extend(symbolCoords)
        cList.extend([r, g, b, 1.0] * len(symbolCoords))
        
        textPos = symbolCoords.mean(axis=0) + [0.5, 0.5, 0.5]
        annos.append([textPos.tolist(), text, (r,g,b, 1.0)])      
                                                 
    GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
    GL.glEnableClientState(GL.GL_COLOR_ARRAY)
 
    vList = array(vList, float).ravel()
    cList = array(cList, float)
 
    GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
    GL.glColorPointer(4, GL.GL_FLOAT, 0, cList)
    GL.glDrawArrays(GL.GL_LINES, 0, len(vList)/3)
    
    GL.glDisableClientState(GL.GL_COLOR_ARRAY)
    GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
    
    
    for point, text, color in annos:
      glColor4f(*color)
      glRasterPos3d(*point)

      for char in text:
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, ord(char))


  def _drawDataTrack(self, nuc, chromo, pos, coordCache, models,
                      dataSymbols, lineWidth, lineSmooth, segments=None):
        
    GL.glLineWidth(lineWidth)

    for model in models:
      coords = coordCache[chromo][model]
      
      if len(coords) < 2:
        continue
      
      for i in range(lineSmooth):
        coords = drawing.smoothPath(coords)

      for group, attrs in dataSymbols:
        if not len(group['regions']):
          continue
 
        symbol = attrs['options'][2]
        
        if not symbol:
          continue
        
        color = attrs['display'][:4]
        scale = attrs['display'][4]
        thrsh = attrs['display'][5]
        
        regions = array(group['regions'], int32)
        values = array(group['values'], float)
        
        
        if segments is not None:
          if not len(segments):
            continue
          
          idx = apiUtil.pairRegionsIntersection(regions, segments)
          
          regions = regions[idx]
          values = values[idx]
        
        if not len(regions):
          continue

        dAngle, nPoints = DATA_TRACK_SYMBOL_PARAMS[symbol]
        
        if symbol == 1 and lineSmooth: # Special case for circle
          f = lineSmooth+1
          dAngle /= f
          nPoints *= f

        symbolCoords = dataLayer.getSymbolCoords(pos, coords,  regions,
                                                 values, dAngle,
                                                 int32(nPoints), int32(lineSmooth), 
                                                 0.2*scale, thrsh, minRadius=0.01,
                                                 maxRadius=2.0)
                                                 
        GL.glColor4f(*color)
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
 
        vertexList = symbolCoords.ravel()
 
        GL.glVertexPointer(3, GL.GL_FLOAT, 0, vertexList)
        GL.glDrawArrays(GL.GL_LINES, 0, len(vertexList)/3)
        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
    

  def _getSurfaceCoords(self, nuc, model, chromosomes, colorArray,
                        radius=1.0, nVoxels=50, structure=None):
    """Create a volumetric representation of the structural surface.
       A voxellated coordinate list in 3D contour order"""
    
    particGroup = nuc._getParticleGroup(structure)
    allCoords = []
    
    for chromo in chromosomes:
      if chromo not in particGroup:
        continue
    
      coords = nuc.getModelCoords(model, [chromo,], structure)
      allCoords.append(coords)
    
    allCoords = vstack(allCoords)
    vertices, colors, normals = apiUtil.calcCoordMesh(allCoords, colorArray, radius, nVoxels)
    
    return vertices, colors, normals

