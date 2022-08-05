from numpy import array, float32, empty, arange, zeros, append, argwhere
from numpy import cos, sin, ones, cross, dot, sqrt, outer, clip, log10, int32
from gui.StructurePanel import StructureOuterPanel, StructurePanel
from math import ceil
from PySide import QtCore, QtGui
from OpenGL import GL, GLU

from gui.Gl3dPanel import Gl3dPanel
from gui.qtgui.ToolBar import ToolBar
from gui.qtgui.Label import Label
from gui.qtgui.Button import Button
from gui.qtgui.Entry import IntRangesEntry, IntEntry, Entry
from gui.qtgui.Menu import Menu
from gui.qtgui.Slider import FloatSlider, Slider
from gui.qtgui.SpinBox import FloatSpinBox, IntSpinBox
from gui.qtgui.Colors import ColorDialog, GradientEditor
from gui.qtgui.Base import Icon

from util.Structure import getArcChromoCoords, getCircleChromoCoords, getLinearChromoCoords

PI = 3.14159265358979323846
TAU = 2.0 * PI

# # To-do # #
# Text as textured rectangles
# Restraint colour scheme to match chromo colour scheme?
# Should use raw contats not binned ones
# Minumum number of chromosomal particles
# Fix colour mode actions
# Check in display type icons
# Separate cis/trans restraint buttons
# No data tracks on tubes or chomo colors
# Setup multiple GL contexts
# Toggle 2D/3D projection (check mouse movement in 2D)
# On change tabs callbacks
# Arc options cis; up/down, trans up/down
# Arcs into cython - variable num segments?


class InteractomeOuterPanel(StructureOuterPanel):
  
  def __init__(self, parent, mainApp, openFunc, **kw):
    
    QtGui.QWidget.__init__(self, parent)
  
    self.mainApp = mainApp
    self.openFunc = openFunc
    self.getIcon = mainApp.getIcon
    
    layout = QtGui.QGridLayout(self)
    layout.setSpacing(2)
    layout.setContentsMargins(2,2,2,2)
    layout.setRowStretch(1, 2)
    layout.setColumnStretch(0, 2)
     
    # Interactome toolbar
    
    icons = ['zoom-in.png', 'zoom-out.png',
             'zoom-center.png', 'move-center.png',
             'display-line.png','display-tube.png',
             'display-restraints.png', 'linear.png',
             'arcs.png', 'circles.png',
             'color-seq.png', 'color-chromo.png',
             'color-data.png','color-faint.png',
             'color-density.png',
             'rmsd.png',]
    
    icons = [self.getIcon(i) for i in icons]
    
    texts = ['Zoom in','Zoom out','Reset zoom','Center',
             'Line display', 'Tube display',
             'Toggle restraints display',
             'Linear chromosome arrangement',
             'Arc chromosome arrangement',
             'Circle chromosome arrangement',
             'Sequence postion colors', 'Chromosome colours',
             'Data track colours', 'Faint colouring',
             'Density colours', 'Structure RMSD colours']
    
    funcs = [self.zoomIn, self.zoomOut, self.zoomReset, self.moveCenter,
             self.setInteractomeLine, self.setInteractomeTube, 
             self.toggleInteractomeRestraints, self.setInteractomeLinear,
             self.setInteractomeArcs, self.setInteractomeCircles,
             self.setInteractColorSeq, self.setInteractColorChromo,
             self.setInteractColorGenData, self.setInteractColorFaint, 
             self.setInteractColorDensity, self.setInteractColorRmsd]

    self.interactToolbar = ToolBar(self, 'Interactome toolbar', funcs, icons, texts, 
                                   objName='interactToolbar', areas='tblr', iconSize=32)
          
    colorConfigButton = Button(self.interactToolbar, ' ', icon=self.getIcon('colors.png'),
                               tipText='Colour settings', iconSize=32)
    colorConfigMenu = Menu(colorConfigButton, 'Colours configure',
                           setupFunc=self._setupColorConfigMenu)
                           
    colorConfigButton.setMenu(colorConfigMenu)
    self.interactToolbar.addWidget(colorConfigButton) 

    layout.addWidget(self.interactToolbar, 0, 0, 1, 2)
    
    # Main 3D widget
                           
    self.innerPanel = InteractomePanel(self, mainApp, openFunc, **kw)

    layout.addWidget(self.innerPanel, 1, 0)
    self.setLayout( layout )


  def updateContents(self):
  
    nuc = self.mainApp.nuc

    self.innerPanel.construct(nuc)
    self.innerPanel.update()
    
    
  def setInteractColorSeq(self):

    self.setInteractomeColorMode(0)


  def setInteractColorChromo(self):

    self.setInteractomeColorMode(1)


  def setInteractColorGenData(self):

    self.setInteractomeColorMode(2)


  def setInteractColorFaint(self):

    self.setInteractomeColorMode(4)


  def setInteractColorDensity(self):

    self.setInteractomeColorMode(5)


  def setInteractColorRmsd(self):

    self.setInteractomeColorMode(6)
  
  
  def toggleInteractomeRestraints(self):
    
    nuc = self.mainApp.nuc
    
    attrs = nuc.display.attrs
    opts = list(attrs['interactome'])      
    showCis    = opts[5]
    showTrans  = opts[6]   
    
    if showCis or showTrans:
      nuc.setInteractomeParams(showCis=False, showTrans=False)
    
    else:
      nuc.setInteractomeParams(showCis=True, showTrans=True)
     
    self.updateContents()

 
  def setInteractomeLine(self):
    
    self.setInteractomeRenderMode(1)
    
 
  def setInteractomeTube(self):
    
    self.setInteractomeRenderMode(2)


  def setInteractomeLinear(self):

    if self.mainApp.nuc.setInteractomeParams(displayMode=0):
      self.updateContents()


  def setInteractomeArcs(self):
  
    if self.mainApp.nuc.setInteractomeParams(displayMode=1):
      self.updateContents()


  def setInteractomeCircles(self):
  
    if self.mainApp.nuc.setInteractomeParams(displayMode=2):
      self.updateContents()


  def setInteractomeColorMode(self, val):
  
    if self.mainApp.nuc.setInteractomeParams(colorMode=val):
      self.updateContents()
  
  
  def setInteractomeRenderMode(self, val):
  
    if self.mainApp.nuc.setInteractomeParams(renderMode=val):
      self.updateContents()
  
  
  def setInteractomeThickness(self, val):
  
    if self.mainApp.nuc.setInteractomeParams(thickness=val): 
       self.updateContents()
 
  
  def setInteractomeIs3D(self, boolean):
  
    if self.mainApp.nuc.setInteractomeParams(is3d=boolean): 
       self.updateContents()
 
  
  def setInteractomeShowCis(self, boolean):
  
    if self.mainApp.nuc.setInteractomeParams(showCis=boolean): 
      self.updateContents()
  
  
  def setInteractomeShowTrans(self, boolean):
  
    if self.mainApp.nuc.setInteractomeParams(showTrans=boolean):  
      self.updateContents()
    
    
class InteractomePanel(StructurePanel):

  def __init__(self, parent, mainApp, openFunc, **kw):
    
    StructurePanel.__init__(self, parent, mainApp, openFunc)
    self.rGlList = 101
    self.iGlList = 102
    self.cGlList = 105
    self.scale = 10.0
    self.seqScale = 100e6  # i.e. self.scale units is 100 Mb
    self._arcCoords = None
  
  
  def getArcCoords(self, nSegments):
  
    if (self._arcCoords is None) or (len(self._arcCoords) != nSegments):
      
      dAngle = PI/nSegments
      angles = [(x, x+dAngle) for x in arange(0.0, PI-dAngle/2, dAngle)]
      angles = array(angles, float).ravel()
      
      coords = empty((2, 2*nSegments), float)      
      coords[0] = (1.0 + cos(angles)) / 2.0
      coords[1] = sin(angles) / 2.0 
      
      self._arcCoords = coords
  
    return self._arcCoords


  def _getDisplayParams(self, nuc):
    
    dispAttrs = nuc.display.attrs
    
    showText = 1
    showScalebar, restColorMode = dispAttrs['options'][5:7]
    
    interactMode, colorMode, dispMode, thickness, is3d, showCis, showTrans = dispAttrs['interactome'][:7]
    
    radius, stickFrac, lineWidth, tubeWidth = [1.0, 0.5, float(thickness), 1.0]
    
    showModels = set([0])
    
    sphDetail, lineSmooth, tubeSmooth, tubeDetail = [1, 0, 1, 1]
  
    chromoParams = {}
    for chromo in nuc.chromosomes:
      chrShow, chrText, chrColorMode, chrDispMode = nuc.getChromoDisplayParams(chromo)
      chromoParams[chromo] = [chrShow, chrText, colorMode, dispMode]
     
    params = [colorMode, dispMode, showCis, showTrans, showText, showScalebar, 
              restColorMode, radius, stickFrac, lineWidth, tubeWidth, showModels, 
              sphDetail, lineSmooth, tubeSmooth, tubeDetail, chromoParams]
    
    return params
   
   
  def _getChromoCoords(self, nuc, coordsGroup, particGroup, models,
                       clip_radius=None, region_of_interest=None):
     
    n = 0
    model = 0
    coordCache = {}
    chromoCenters = {}
    idx_segments = {}
    seq_segments = {}
    pos_cache = {}
    
    displayMode = nuc.display.attrs['interactome'][0]
    
    chromos = nuc.getDisplayedChromosomes(nuc.structure.name)
    
    offsets = {}
    totalSeqLen = 0
    maxLen = 0
    sortChromos = []
    
    for chromo in chromos:
      start, end = nuc.getChromosomeLimits(chromo)
      delta = end - start
      if not delta:
        continue
      
      delta  *= 1.05
      sortChromos.append(chromo)
      offsets[chromo] = totalSeqLen
      totalSeqLen += delta
      if delta > maxLen:
        maxLen = delta    
    
    nChromos = len(sortChromos) 
    for i, chromo in enumerate(sortChromos):
    
      seqPos = array(particGroup[chromo]['positions'], int32)
      pos_cache[chromo] = seqPos
      nPos = len(seqPos)
      
      coordCache[chromo] = {}
      if displayMode == 0:
        coords, centre = getLinearChromoCoords(self.seqScale, seqPos, i, nChromos, self.scale)
      
      elif displayMode == 1:
        coords, centre = getArcChromoCoords(seqPos, offsets[chromo], totalSeqLen, self.scale)

      else:
        coords, centre = getCircleChromoCoords(seqPos, i, nChromos, maxLen, self.scale)
      
      chromoCenters[chromo] = centre
      n += len(coords)
        # default segments whole region
        
      idx_segments[chromo] = array([[0, nPos-1],], int32)
      seq_segments[chromo] = array([[seqPos[0], seqPos[-1]],], int32)
      
      # restrict coords according to section selections
      #sectionCentre, sectionLength = nuc.getChromoSectionDisplayParams(chromo)
      
      #if sectionLength > 0:
      #  start = max(0, int(nPos * sectionCentre - sectionLength/2))
      #  end = min(nPos, int(nPos * sectionCentre + sectionLength/2))
      #  coords = coords[start:end]

      coordCache[chromo][model] = coords
    
    clipLimit = 2048
    origin = array([0.0, 0.0, 0.0])
       
    return n, coordCache, origin, clipLimit, chromoCenters, idx_segments, seq_segments, pos_cache
    
  def updateQtLayer(self, painter):
    
    nuc = self.mainApp.nuc
    
    if not nuc:
      return
    
    chrs = self.chromoCenters.keys()
    nChrs = len(chrs)
    #size of canvas
    w0, h0 = self.glGeometry[:2]
        
    if nChrs == 1:
      chromoText = 'Chromosome : %s' % (chrs[0])
      
    elif nChrs:
      chromoText = 'Chromosomes : %s' % (','.join(chrs))
    
    else:
      return    

    painter.setRenderHint(painter.Antialiasing)
    
    painter.setFont(self.fontL)
    white = QtGui.QColor(255, 255, 255, 255)
    
    fm = painter.fontMetrics()
    p = 5.0

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
      painter.setBrush(QtGui.QColor(0, 0, 0, 128))
      painter.drawRoundedRect(x-p, y-bh/2.0-p, bw+p+p, bh+p+p, p, p)
      
      painter.setPen(white)
      painter.drawText(x, y+bh/2.0, chrA)  
      
      
  def _drawRestraints(self, nuc, coordCache, models, cis, trans, colouring,
                      alpha, segments=None, lineWidth=2.0, nColors=10):
    
      
    chromosomes = list(coordCache.keys())
    
    # Will use currently active structure's restraints and particles
    restraintDict = nuc.getRestraints(chromosomes)
    particGroup = nuc._getParticleGroup()
    
    scheme = nuc.getColorScheme('restraint')
    colorArray = nuc.calcInterpolatedColors(scheme, nColors, 0.8)
    
    nSegments = 20
    arcCoords = self.getArcCoords(nSegments)
    arcX = arcCoords[0]
    arcY = arcCoords[1]
    
    #GL.glLineWidth(lineWidth)            
    GL.glLineWidth(1.0)            
    
    for chromoPair in restraintDict:
      chromoA, chromoB = sorted(chromoPair)
      
      if not cis and (chromoA == chromoB):
        continue
      
      if not trans and (chromoA != chromoB):
        continue
       
      chrShow, chrText, chrColorMode, chrDispMode = nuc.getChromoDisplayParams(chromoA)
      if not chrShow:
        continue

      chrShow, chrText, chrColorMode, chrDispMode = nuc.getChromoDisplayParams(chromoB)
      if not chrShow:
        continue
      
      restraints = restraintDict[chromoPair]
      nRestraints = len(restraints)
      
      if not nRestraints:
        continue
      
      rCoords = empty((nRestraints, 2*nSegments, 3), float)
      idxA = array(restraints[:,0], int)
      idxB = array(restraints[:,1], int)
      uppers = array(restraints[:,5], float)
      
      modelCoordsA = coordCache[chromoA]
      modelCoordsB = coordCache[chromoB]
      
      m = models[0]
      nA = len(modelCoordsA[m])
      nB = len(modelCoordsB[m])
      
      # TBD fix this
      idxA[(idxA >= nA).nonzero()] = nA -1
      idxB[(idxB >= nB).nonzero()] = nB -1


      zDir = -1.0 if chromoA == chromoB else 1.0
      nIdx = len(idxA)
      
      # Setup colours
      
      if colouring == 2:
        chrPos = array(particGroup[chromoA]['positions'])
      
      rColors = zeros((nRestraints, 2*nSegments, 4), float)
      if colouring == 0: # Cis/Trans colours
        if chromoA == chromoB:
          data = colorArray[0]
        else:
          data = colorArray[-1]
        
        rColors[:,:] = data
 
      elif colouring == 2: # Seq sep. colours
 
        if chromoA == chromoB:
          deltas = log10(clip( abs(chrPos[idxA]-chrPos[idxB]), 1e4, 1e8 )) - 4.0
          cIndices = array(deltas/4.0, int)
          data = colorArray[cIndices]
          rColors[:,:] = data.reshape(nIdx,1,4)
 
        else:
          rColors[:,:] = colorArray[-1]
 
      elif colouring == 3: # Chromosome colours
        scheme = [nuc.getChromoColor(chromoB).rgb,
                  nuc.getChromoColor(chromoA).rgb]
        colorArray = nuc.calcInterpolatedColors(scheme, 2*nSegments, 1.0)
        rColors[:] = colorArray
                     
      for model in models:
      
        coordsA = modelCoordsA[model]
        coordsB = modelCoordsB[model]
        rA = coordsA[idxA]
        rB = coordsB[idxB]
        deltas = rB-rA
        dists = sqrt((deltas*deltas).sum(axis=1))
        
        if colouring == 1: # Dist colours
          deltas = rA-rB
          dists1 = clip((dists-uppers)/uppers, 0.0, 4.0)
          cIndices = array(dists1*(nColors-1)/4.0, int)
          data = colorArray[cIndices]          
          rColors[:,:] = data.reshape(nIdx,1,4)

        startX = rA[:,0].reshape(nIdx, 1)
        startY = rA[:,1].reshape(nIdx, 1)        
        
        rCoords[:,:,0] = startX + outer(deltas[:,0], arcX)
        rCoords[:,:,1] = startY + outer(deltas[:,1], arcX)
        rCoords[:,:,2] = outer(zDir * sqrt(dists), arcY)
        
        vertexList = rCoords.ravel()
        colourList = rColors.ravel()
        
        GL.glLineStipple(1, 43690) # int('0b1010101010101010', 2)
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
    
      
  def old_2d_construct(self, nuc=None, binSize=int(1e6), cType='singleCell'):
    
    for i in self.glLists:
      GL.glDeleteLists(i, 1)
    
    self.glLists = []
    
    self.lighting = set()
    radius = 10.0
    
    glList = self.firstGlList
    if nuc and cType in nuc.contacts:      
      
      chromos = nuc.getChromosomes()
      
      tSize = 0
      limits = {}
      for chromo in chromos:
        start, end = nuc.getChromosomeLimits(chromo)
        limits[chromo] = start, end
        delta = end - start
        tSize += delta
        
      nBins = int(ceil(tSize/float(binSize))) + len(chromos)
      dAngle = TAU/nBins
      angle = 0.0
      lighting = False
      cAngles = {}
      GL.glNewList(glList, GL.GL_COMPILE)
      self.glLists.append(glList)
        
      for chromo in chromos:
        cAngles[chromo] = angle
        start, end = limits[chromo]
        delta = end - start
        cBins = int(delta/float(binSize)) + 1
        eAngle = angle + float(dAngle * cBins)
        #angles = arange(angle, eAngle, dAngle)

        #vList = zeros((len(angles), 3), float)
        #vList[:,0] = radius * cos(angles)
        #vList[:,1] = radius * sin(angles)
        #vList = vList.ravel()
        
        
        chrColor = nuc.getChromoColor(chromo).rgba
        GL.glLineWidth(4.0)
         
        OpenGlUtil.createArc(angle, dAngle, cBins, radius, chrColor)
        
        """
        GL.glColor4d(*chrColor)
 
        GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
 
        GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
 
        GL.glDrawArrays(GL.GL_LINE_STRIP, 0, len(vList)/3)
 
        GL.glDisableClientState(GL.GL_VERTEX_ARRAY)

        GL.glDisable(GL.GL_NORMALIZE)
        """
 
 
        angle = eAngle 
        
      if lighting:
        self.lighting.add(glList)
      
      GL.glEndList()
      glList += 1
      
      contacts = nuc.contacts
      subGroup = nuc.getCachedContacts(cType)
      chromoKeys = []
      
      
      for i, chromoA in enumerate(chromos):
        keyCis = ' '.join((chromoA, chromoA))
        if keyCis in subGroup:
          chromoKeys.append(keyCis)
        
        for chromoB in chromos[i+1:]:
          keyTrans = ' '.join(tuple(sorted([chromoA, chromoB])))
          
          if keyTrans in subGroup:
            chromoKeys.append(keyTrans)
       
      for chromoKey in chromoKeys:
        contacts = array(subGroup[chromoKey])
      
        posA = contacts[:,0] 
        posB = contacts[:,1]
        
        posA /= float(binSize)
        posB /= float(binSize)
        
        diffs = posA-posB
        idx = (diffs != 0).nonzero()
        
        posA = posA[idx]
        posB = posB[idx]
        
        chrA, chrB = chromoKey.split()
        
        startA, endA = limits[chrA]
        startB, endB = limits[chrB]
        
        angleA = cAngles[chrA]
        angleB = cAngles[chrB]
      
        
        anglesA = angleA + dAngle * posA
        anglesB = angleB + dAngle * posB
        
        nLines = len(anglesA)
        
        chrColorA = nuc.getChromoColor(chrA).rgba
        chrColorB = nuc.getChromoColor(chrB).rgba
        
        useCurves = False
        
        if useCurves:
          starts = zeros((nLines, 3), float)
          ends   = zeros((nLines, 3), float)
          starts[:,0] = radius * cos(anglesA)
          starts[:,1] = radius * sin(anglesA)
          ends[:,0]   = radius * cos(anglesB)
          ends[:,1]   = radius * sin(anglesB)
          
          mids = (starts + ends) / 2.0
          v1 = starts - mids
          sizes = sqrt((v1*v1).sum(axis=1))
          sizes[(sizes == 0).nonzero()] = 1.0
          
          v1[:,0] /= sizes
          v1[:,1] /= sizes
          v1[:,2] /= sizes
          
          zAxis = array([0.0, 0.0, -1.0])
          nAngles = 10
          dAngle2 = PI/(nAngles-1)
          
          vList = zeros((nAngles, 3), float)
          cList = zeros((nAngles, 4), float)
          
          glList += 1
          GL.glNewList(glList, GL.GL_COMPILE)
          self.glLists.append(glList)
          
          GL.glLineWidth(2.0)
          GL.glColor4f(1.0, 0.0, 0.0, 1.0)
          
          GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
          GL.glEnableClientState(GL.GL_COLOR_ARRAY)

          for i in range(nLines):
            v = v1[i]
            m = mids[i]
            s = sizes[i]

            axis = cross(v, zAxis)
            rMat = getRotationMatrix(axis, dAngle2)
          
            for j in range(nAngles):
              f = j/float(nAngles)
              
              vList[j] = m + v * s
              cList[j] = f * chrColorB + (1.0-f) * chrColorA
              
              v = dot(v, rMat)
            
            GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList.ravel())
            GL.glColorPointer(4, GL.GL_FLOAT, 0, cList.ravel())
 
            GL.glDrawArrays(GL.GL_LINE_STRIP, 0, len(vList))
          
          GL.glDisableClientState(GL.GL_COLOR_ARRAY)
          GL.glDisableClientState(GL.GL_VERTEX_ARRAY)
            
          GL.glDisable(GL.GL_NORMALIZE)
          GL.glEndList()
         
        else:
          vList = zeros((2*nLines, 3), float)
          vList[::2,0]  = radius * cos(anglesA)
          vList[::2,1]  = radius * sin(anglesA)
          vList[1::2,0] = radius * cos(anglesB)
          vList[1::2,1] = radius * sin(anglesB)
          
          cList = ones((len(vList), 4), float)
          cList[0::2] = chrColorA
          cList[1::2] = chrColorB
 
          vList = vList.ravel()
          cList = cList.ravel()
 
          glList += 1
          GL.glNewList(glList, GL.GL_COMPILE)
          self.glLists.append(glList)
 
          chrColor = nuc.getChromoColor(chromo).rgba
 
          GL.glLineWidth(2.0)
          GL.glEnableClientState(GL.GL_VERTEX_ARRAY)
          GL.glEnableClientState(GL.GL_COLOR_ARRAY)
 
          GL.glColorPointer(4, GL.GL_FLOAT, 0, cList)
          GL.glVertexPointer(3, GL.GL_FLOAT, 0, vList)
 
          GL.glDrawArrays(GL.GL_LINES, 0, len(vList)/3)
 
          GL.glDisableClientState(GL.GL_COLOR_ARRAY)
          GL.glDisableClientState(GL.GL_VERTEX_ARRAY)

 
          GL.glDisable(GL.GL_NORMALIZE)
          GL.glEndList()
        
       
    else:
      return Gl3dPanel.construct(self)

