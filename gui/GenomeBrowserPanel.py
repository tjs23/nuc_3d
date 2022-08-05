import sys, time
sys.path.append('../')

from math import log, ceil, floor
from gui.qtgui.Frame import Frame
from PySide import QtCore, QtGui

Qt = QtCore.Qt
QImage = QtGui.QImage
QPoint = QtCore.QPoint
LeftButton = QtCore.Qt.LeftButton
MiddleButton = QtCore.Qt.MiddleButton
RightButton = QtCore.Qt.RightButton

# Toolbar options:
#   Flip diagonal direction - icon
#   Regions coords could locate
#   Find - icon/menu
#     (centre) position
#     annotation (e.g. gene) - autocomplete?

# Annotation display
#   ? on top image at suitable scale, (45 degree?)

# If really zoomed in
#  Chr label takes MB, tick only show Kb

# Contact box display type
#   ContactMap.draw() - adds a GL_LINES vertex array for squares coords
#   Same GL List is OK for now
#   Boxes that are too small < 2 pix are ignored
#   Chromo dataTrack regions must map to [0,1]
#   Know already which chromo regions viewed
#    - sum total length
#    - calc starts and ends as fracs of total
#    - data region is a fraction of chomo sub-fraction with start offset
#    - top-left, bottom-right
#    - GL_TRIANGLE_STRIP
#    - fill in Cython
# 
#   Try semi-transparent

# Corner
#   "Chr:", Names of tracks in native colour, on black bg

class GenomeBrowserPanel(Frame):

  def __init__(self, parent, mainApp, regionCallback, isVertical=False, **kw):

    Frame.__init__(self, parent, **kw)
    
    self.mainApp = mainApp
    self.isVertical = isVertical
    self.reverse = False
    
    if isVertical:
      self.numPix = self.height()
    else:
      self.numPix = self.width()  
    
    self.viewRegion = (0.0, 1.0)
    self.trackWidth = 25
    self.prevPos = QPoint()
    self._trackRegions = []
    self.regionCallback = regionCallback

    fontMetric = QtGui.QFontMetricsF(self.font())
    self.textHeight = fontMetric.boundingRect('A').height()
    self.minSize = QtCore.QSize(2*self.textHeight, 2*self.textHeight)
       
    self.setMinimumSize(self.minSize)

    self.setStyleSheet("background-color: rgb(0, 0, 0);")
    self.setAutoFillBackground(True)
  

  def _updateRegionCoords(self, seqRegions):  
    
    if seqRegions:
      chromo1, pix1, nPix, seq1A, seq1B, j = seqRegions[0]
      chromo2, pix1, nPix, seq2A, seq2B, j = seqRegions[-1]
    
    else:
      chromo1, seq1A, chromo2, seq2B = ('', 0, '', 0)
    
    self.regionCallback(chromo1, seq1A, chromo2, seq2B)
    
  
  def _getScreenDataTrack(self, x, y):
     
    if self.isVertical:
      p = self.width() - x
    else:
      p = y
 
    for start, end, source, code in self._trackRegions:
      if start <= p < end:
        return source, code
        
  
  def wheelEvent(self, event):
    
    if event.delta() < 0:
      z = 1.2
    else:
      z = 0.8333333333333333
    
    mods = event.modifiers()
    haveCtrl = mods & Qt.CTRL
    haveShift = mods & Qt.SHIFT
    
    if haveCtrl or haveShift:
      nuc = self.mainApp.nuc      
      key = self._getScreenDataTrack(event.x(), event.y())
      
      if key is not None:
        source, code = key
      
        if haveCtrl:
          scale = nuc.getDataTrackScale(source, code)
          scale = min(10.0, max(0.01, scale*z))
          nuc.setDataTrackScale(source, code, scale)
        else:
          minVal, maxVal = nuc.getDataTrackThresholds(source, code)
          thrsh = min(1.0, max(0.0, minVal*z))
          nuc.setDataTrackThresholds(source, code, thrsh)
            
      self.parent().contactMap.draw()  
    
    else:  
      self.parent().zoom(z)
  
    event.accept()
    
    
  def mousePressEvent(self, event):
    
    QtGui.QWidget.mousePressEvent(self, event)
  
    pos = event.pos()
    self.prevPos = QPoint(event.pos())


  def mouseMoveEvent(self, event):
    
    QtGui.QWidget.mouseMoveEvent(self, event)
    
    dx = event.x() - self.prevPos.x()
    dy = event.y() - self.prevPos.y()
    buttons = event.buttons()
    pos = event.pos()
    
    if (buttons & LeftButton) or (buttons & MiddleButton):
      mods = event.modifiers()
      haveCtrl = mods & Qt.CTRL
      haveShift = mods & Qt.SHIFT
      
      if self.reverse:
        dy *= -1
      
      if haveCtrl or haveShift:
        self.parent().translate(-dx, -dy)
      else:
        self.parent().translate(dx, dy)
      
    elif buttons & RightButton:
      pass
      
    self.prevPos = QPoint(event.pos())
  
  
  def _getTickDelta(self, r0, r1, maxTicks):
    
    s0 = float(min(r0, r1))
    s1 = float(max(r0, r1))
    d = s1 - s0
    w = 0.001 * d
  
    n = int(floor(log(0.999*d)/log(10)))
    delta = pow(10.0, n)
    t = 0
    
    if 5*delta < d:
      delta = 5 * delta
      t = 2
      
    elif 2*delta < d:
      delta = 2 * delta
      t = 1
 
    nticks = int(floor((s1-s0-2*w) / delta))
    d = delta
    
    while True:
      if t == 2:
        d = (2 * d) / 5
      else:
        d = d / 2
        
      nticks = int(floor((s1-s0-2*w) / d))
      if nticks > maxTicks:
        break
        
      if t == 0:
        n = n - 1
        
      delta = d
      t = (t - 1) % 3
    
    return delta
    
    
  def _findTicks(self, r0, r1, width, delta=None, maxTicks=25, minDecimals=0):
 
    s0 = float(min(r0, r1))
    s1 = float(max(r0, r1))
    dd = float(max(abs(r0), abs(r1)))
 
    if not delta:
      delta = self._getTickDelta(r0, r1, maxTicks)

    n0 = int( ceil((s0 + 0.001 * (s1 - s0))/delta))
    
    if s0 == 0.0:
      n0 -= 1
    
    n1 = int(floor((s1 - 0.001 * (s1 - s0))/delta))
    n = int(floor(log(delta) / log(10)))
 
    tickPos = [ x * delta for x in xrange(n0, n1+1) ]
        
    if dd >= 1000:
      nn = int(floor(log(0.999*dd)/log(10)))
      nDecimals = max(0, nn-n)
      tickFormat = '%%.%de' % nDecimals
    
    else:
      nDecimals = max(minDecimals, -n)
      tickFormat = '%%.%df' % nDecimals
    
    coords = [width * (t-r0) / (r1-r0) for t in tickPos]
    last = len(coords)-1
    ticks = []
    
    for i, pos in enumerate(tickPos):
      coord = coords[i]
      if i == last:
        span = width-coord
        
      else:
        span = coords[i+1] - coord
      
      text = tickFormat % pos
      ticks.append([text, coord, int(ceil(span))])      
    
    if ticks and ticks[0][1] > 1:
      ticks.insert(0, ['', 0.0, ticks[0][1]])
    
    return ticks
    
    
  def _getDataTracks(self, chromoRegions):
  
    nuc = self.mainApp.nuc
    pixmapCache = self.parent().pixmapCache
    trackPixmaps = []
    
    for source, code in nuc.getDisplayedDataTracks():
      if nuc.getDataTrackPeakType(source, code) > 1:
        # Drawn on contact map not genome browser
        continue
      
      pixmaps = []
      for chromo, pos, span, start, end, b in chromoRegions:
        key = (source, code, chromo, span, start, end)
        
        if key in pixmapCache:
          pixmap = pixmapCache[key]
          
        else:
          pixmap = nuc.getDataTrackPixmap(source, code, chromo, start, end,
                                          span, self.trackWidth)
          pixmapCache[key] = pixmap
 
        pixmaps.append((pos, pixmap))
         
      trackPixmaps.append((pixmaps, source, code))  
         
    if len(pixmapCache) > 2000:
      self.parent().pixmapCache = {}
    
    return trackPixmaps
    
    
  def _getChromoRegions(self):
    
    nuc = self.mainApp.nuc
    p1, p2 = self.viewRegion
    chromoRegions = []
    tickRegions = []

    if nuc:
      chromosomes = nuc.getDisplayedChromosomes()
      
      if not chromosomes:
        return [], []
      
      nPix = self.numPix
      chromoSizes = {}
      sortChromos = []
      chromoStarts = {}
      chromoSeqStarts = {}
      
      seqLen = 0
      for chromo in chromosomes:
        start, end = nuc.getChromosomeLimits(chromo)
        size = (end-start)
        if not size:
          continue
        
        sortChromos.append(chromo)
        chromoSizes[chromo] = size
        chromoStarts[chromo] = seqLen
        chromoSeqStarts[chromo] = start
        seqLen += size
       
      genFrac = p2-p1
      binSize = (genFrac * seqLen) / nPix # bp per pixel
      pixStart = p1 * nPix / genFrac
      pixEnd = pixStart + nPix  

      for chromo in sortChromos:
        chromoSizes[chromo] /= binSize
        chromoStarts[chromo] /= binSize
      
      delta = self._getTickDelta(0.0, (genFrac*seqLen)/1e6, 25)
      
      pos = 0.0
      i = 0
      j = 0
      for chromo in sortChromos:
        start = chromoStarts[chromo]
        if start > pixEnd:
          continue  
        
        size = chromoSizes[chromo]
        end = start + size
        if end < pixStart:
          continue
        
        if start < pixStart:
          c1 = pixStart-start
        else:
          c1 = 0  
        
        if end > pixEnd:
          c2 = size - (end-pixEnd)
        else:
          c2 = size        
        
        span = int(c2-c1)
        if span < 1:
          continue
          
        seqStart = chromoSeqStarts[chromo]
        if delta < 100:
          if span >= 50:
            seqA = (seqStart + c1 * binSize) * 1e-6
            seqB = (seqStart + c2 * binSize) * 1e-6
          
            ticks = self._findTicks(seqA, seqB, span, delta)
            for text, coord, subSpan in ticks:
              tickRegions.append( (text, pos+coord, subSpan, i%2) )
              i += 1
 
          else:
            tickRegions.append( ('', pos, span, i%2) )
            i += 1
        
        chromoRegions.append( (chromo, int(pos), int(span), seqStart+int(c1*binSize), seqStart+int(c2*binSize), j%2))
        pos += span
        j += 1
       
    return tickRegions, chromoRegions
 
 
  def minimumSizeHint(self):
  
    return self.minSize


  def sizeHint(self):
    
    return self.minSize
  
  
  def setRegion(self, p1, p2, numPix, reverse=False):
  
    self.viewRegion = (p1, p2)
    self.numPix = numPix
    self.reverse = reverse
    
    if self.isVisible():
      self.update()
 
 
  def resizeEvent(self, event):
    
    return Frame.resizeEvent(self, event)


  def paintEvent(self, event):
    
    self.draw()
    
    return Frame.paintEvent(self, event)
     
     
  def draw(self):
  
    painter = QtGui.QPainter()
    painter.begin(self)
    fontMetric = QtGui.QFontMetricsF(painter.font())
    
    drawRect = painter.drawRect
    drawText = painter.drawText
    drawLine = painter.drawLine
    
    brushA = QtGui.QBrush(QtGui.QColor(48, 48, 64, 255))
    brushB = QtGui.QBrush(QtGui.QColor(64, 48, 48, 255))
    brushes = (brushA, brushB)
    penA = QtGui.QPen(QtGui.QColor(160, 160, 160, 255))
    penB = QtGui.QPen(QtGui.QColor(64, 64, 64, 255))
    penC = QtGui.QPen(QtGui.QColor(255, 255, 255, 255))

    tickRegions, chromoRegions = self._getChromoRegions()
    self._updateRegionCoords(chromoRegions)
    
    trackPixmaps = self._getDataTracks(chromoRegions)    
    
    bbox = fontMetric.boundingRect('A')
    thick = bbox.height()
    pos = 0
    span = self.width()
    
    totalSize = 2*thick + (self.trackWidth+2)*len(trackPixmaps)
   
    if self.isVertical:
      if self.reverse:
        painter.rotate(270)
        painter.translate(-self.height(), 0)
      
      else:
        painter.rotate(90)
        painter.translate(0, -totalSize)
      
      self.setFixedWidth(totalSize)
      
    else: 
      if self.reverse:
        painter.rotate(180)
        painter.translate(totalSize, 0)
      
      self.setFixedHeight(totalSize)
    
    if tickRegions:
      for text, pos, span, c in tickRegions:

        bbox = fontMetric.boundingRect(text)
        th = bbox.height()
        tw = bbox.width()

        brush = brushes[c]
        painter.setBrush(brush)
        painter.setPen(penB)

        drawRect(pos, 0, span, thick)
        painter.setPen(penA)
 
        if pos > 0:
           drawLine(pos, 0, pos, thick)
 
        if span > tw:
          drawText(pos+2, thick-3, text)
 
      for text, pos, span, a, b, c in chromoRegions:

        bbox = fontMetric.boundingRect(text)
        th = bbox.height()
        tw = bbox.width()

        brush = brushes[c]
        painter.setBrush(brush)
        painter.setPen(penB)
        drawRect(pos, thick, span, thick)
        
        if pos > 0:
          painter.setPen(penA)
          drawLine(pos, thick, pos, 2*thick)
        
        painter.setPen(penB)
        drawLine(pos, thick, pos+span, thick)
 
        if span > tw-2:
          painter.setPen(penC)
          drawText(pos+span/2-tw/2, 2*thick-2, text)
    
    else:
      for text, pos, span, a, b, c in chromoRegions:

        bbox = fontMetric.boundingRect(text)
        th = bbox.height()
        tw = bbox.width()

        brush = brushes[c]
        painter.setBrush(brush)
        painter.setPen(penB)
        drawRect(pos, 0, span, 2*thick)
        
        if pos > 0:
          painter.setPen(penA)
          drawLine(pos, 0, pos, 2*thick)
 
        if span > tw-2:
          painter.setPen(penC)
          drawText(pos+span/2-tw/2, thick+th/2, text)
    
    self._trackRegions = []
    totWidth = pos + span
    trackPixmaps = self._getDataTracks(chromoRegions)    
    painter.setPen(penB)
    
    y = 2*thick + 1
    for pixmaps, source, code in trackPixmaps:
      for x, pixmap in pixmaps:
        h, w = pixmap.shape[:2]
        qImage = QImage(pixmap.data, w, h, QImage.Format_ARGB32)
        painter.drawImage(x, y, qImage)
      
      self._trackRegions.append( (y, y+self.trackWidth, source, code))
      y += self.trackWidth + 2
      drawLine(0, y-2, totWidth, y-2)
      
    painter.end()

