from math import log
from PySide import QtGui, QtCore
Qt = QtCore.Qt
QPointF = QtCore.QPointF
QRectF = QtCore.QRectF
QColor = QtGui.QColor

from gui.qtgui.Base import Base, Icon
from gui.qtgui.Colors import ColorDialog
from gui.qtgui.FileSelect import selectSaveFile, FileType
from gui.qtgui.graphicsItems.LabelItems import LabelItem, MovableLabelItem
from gui.qtgui.InputDialog import askString
from gui.qtgui.Menu import Menu
from gui.qtgui.Print import PrintDialog

GREY_PEN = QtGui.QPen(QColor(128, 128, 128), 0.8)
TRANSPARENT = QColor(0, 0, 0, 0)
HIGHLIGHT_PEN = QtGui.QPen(QColor(255, 255, 255, 64), 2, Qt.DotLine)
HIGHLIGHT_BRUSH = QColor(255, 255, 255, 32)
NULL_RECT = QRectF()
NULL_POINT = QPointF()

class LegendItem(QtGui.QGraphicsItem):
  
  def __init__(self, parent):

    QtGui.QGraphicsItem.__init__(self)
    
    self.parent = parent
    self.setZValue(2)
    self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
    self.setAcceptedMouseButtons(Qt.LeftButton)
    self.fontMetric = QtGui.QFontMetricsF(self.parent.font())
    self.parent.scene().addItem(self)

  
  def setMatrix(self, matrix, boxSize, posColor, zeroColor, negColor):
    
    r1, g1, b1, a1 = posColor.getRgb()
    r2, g2, b2, a2 = zeroColor.getRgb()
    r3, g3, b3, a3 = negColor.getRgb()
    
    minVal = matrix.getMinValue()
    maxVal = matrix.getMaxValue() or 1.0
    
    # Rough step
    deltaVal0 = (maxVal-minVal)/20.0
    
    if not deltaVal0:
      deltaVal0 = 1.0
    
    # Neat, rounded step
    order = int(round(log(deltaVal0, 10)))
    deltaVal = round(deltaVal0,-order) + 10**order
    
    # Start from max, always go to at least zero
    start = round(maxVal/deltaVal) * deltaVal
    end   = min(0.0, (round(minVal/deltaVal)) * deltaVal)
    
    zDiff = max(abs(minVal),maxVal)
    bbox = self.fontMetric.tightBoundingRect
    
    t = bbox('x').height()
    nBlocks = int((start-end)/deltaVal)+1
    x = boxSize + t
    y = boxSize/2.0 + t/2.0 
    z = start
    
    if -5 < order < 6:
      textFormat  = '%%.%df' % max(1,-order)
    else:
      textFormat = '%.2e'
    
    rects = []
    texts = []
    widths = []
    colors = []
    for i in range(nBlocks):
      
      if z < 0:
        p = -z/zDiff
        p = min(1.0, max(0.0, p))
        q = 1.0 - p
        r = (p * r3) + (q * r2)
        g = (p * g3) + (q * g2)
        b = (p * b3) + (q * b2)
        text = textFormat % z
        
      else:
        p = z/zDiff
        p = min(1.0, max(0.0, p))
        q = 1.0 - p
        r = (p * r1) + (q * r2)
        g = (p * g1) + (q * g2)
        b = (p * b1) + (q * b2)
        text = ' ' + textFormat % z

      rects.append(QRectF(0, i * boxSize, boxSize, boxSize))
      texts.append((x, y, text))
      widths.append(bbox(text).width())
      colors.append(QColor(r, g, b))
      
      z -= deltaVal
      y += boxSize
    
    self.boundingRegion = QRectF(0, 0, boxSize+t+max(widths), y) 
    self.drawData = (rects, texts, colors) 
    self.update()
    
    
  def boundingRect(self):
    
    if self.boundingRegion:
      return self.boundingRegion # .adjust(-2,-2, 2, 2)
    
    else:
      return NULL_RECT
      
      
  def paint(self, painter, option, widget):
    
    drawRect = painter.drawRect
    drawText = painter.drawText
    setPen   = painter.setPen
    setBrush = painter.setBrush
    
    rects, texts, colors = self.drawData
    
    setPen(GREY_PEN)
    for i, rect in enumerate(rects):
      setBrush(colors[i])
      drawRect(rect)
      
    setPen(Qt.black)
    for t in texts:
      x, y, text = t
      drawText(x, y, text)
 
class MatrixItem(QtGui.QGraphicsItem):

  def __init__(self, parent, matrix=None,
               boxSize=20, valLimit=None, showGrid=True):

    QtGui.QGraphicsItem.__init__(self)

    self.parent = parent
    self.matrix = None
    self.valLimit = valLimit
    self.grid = showGrid
    self.boxSize = boxSize
    self.colors = (QColor('#0080FF'),
                   QColor('#000000'),
                   QColor('#FF8000'))
    
    self.gridPen =  QtGui.QPen(QColor(64, 64, 64), 0.5)               
    self.setZValue(1)
    self.fontMetric = QtGui.QFontMetricsF(self.parent.font())
    
    self.xAxisItems = [] # Set according to matrix
    self.yAxisItems = []
    
    self.boundingRegion = NULL_RECT
    self.plotRegion = NULL_RECT
    self.cursor = False
    self.cursorItem = MovableLabelItem(self.parent, pen=GREY_PEN)
    self.legendItem = LegendItem(self.parent)
    self.setAcceptHoverEvents(True)
    self.parent.scene().addItem(self)
      
    if matrix:
      self.setMatrix(matrix)
      
  
  def getPlotCoords(self):
  
    return self.plotRegion


  def getBoundingCoords(self):

    rect = self.boundingRegion
    
    x1 = rect.x()
    y1 = rect.y()
    
    w = rect.width()
    h = rect.height()

    x2 = x1+w
    y2 = y1+h
    
    cx = 0.5 * (x1 + x2)
    cy = 0.5 * (y1 + y2)
      
    return x1, y1, x2, y2, cx, cy
    
    
  def setColors(self, posColor, zeroColor, negColor):
  
    self.colors = (posColor, zeroColor, negColor)
    self.setMatrix()
  
  
  def setBoxSize(self, size):
    
    self.boxSize = size
    self.setMatrix(self.matrix) # have to refresh label rotation
  
  def setLegend(self, bool):
    
    if bool:
      self.legendItem.show()
    else:
      self.legendItem.hide()
  
  def setCoords(self, bool):
    
    if bool:
      self.cursorItem.show()
    else:
      self.cursorItem.hide()
    
  def setGrid(self, bool):
    
    self.grid = bool
    self.update()
  
  def setValueLimit(self, valLimit):
    
    self.valLimit = valLimit
    self.update()

  
  def setMatrix(self, matrix=None):
    
    s = self.boxSize
    posColor, zeroColor, negColor = self.colors
    
    r1, g1, b1, a1 = posColor.getRgb()
    r2, g2, b2, a2 = zeroColor.getRgb()
    r3, g3, b3, a3 = negColor.getRgb()
    
    bbox = self.fontMetric.tightBoundingRect
    
    if matrix:
      self.matrix = matrix
      parent = self.parent

      scene = parent.scene()
      
      data = matrix.matrix
      n = len(data)
      m = len(data[0])
      
      xWidths = [bbox(t).width() for t in matrix.xValues]
      
      if max(xWidths) > s-2:
        angle = 45
      else:
        angle = 0  
     
      if not self.xAxisItems:
        self.xAxisItems = [LabelItem(parent) for i in range(n)]
        self.yAxisItems = [LabelItem(parent) for j in range(m)]
      
      else:
        while len(self.xAxisItems) > n:
          item = self.xAxisItems.pop()
          scene.removeItem(item)
 
        while len(self.yAxisItems) > m:
          item = self.yAxisItems.pop()
          scene.removeItem(item)
        
        while len(self.xAxisItems) < n:
          self.xAxisItems.append(LabelItem(parent))
 
        while len(self.yAxisItems) < m:
          self.yAxisItems.append(LabelItem(parent))
      
      for i, xValue in enumerate(matrix.xValues):
        item = self.xAxisItems[i]
        item.setRotation(angle)
        item.setText(xValue)

      for j, yValue in enumerate(matrix.yValues):
        self.yAxisItems[j].setText(yValue)
        
    else:
      # Matrix is unchanged
      # Axes are also unchanged
      # Colors, valLimit may change
      matrix = self.matrix
      data = matrix.matrix
      n = len(data)
      m = len(data[0])
      xWidths = [bbox(t).width() for t in matrix.xValues]
      
      if max(xWidths) > s-2:
        angle = 45
      else:
        angle = 0  
      
      for i, xValue in enumerate(matrix.xValues):
        self.xAxisItems[i].setRotation(angle)
    
    
    if self.valLimit is None:
      zMax = float(matrix.getValueLimit())
    else:
      zMax = float(self.valLimit)
    
    if not zMax:
      zMax = 1.0
    
    textH = bbox('X').height()
    x0 = y0 = 2.0 * textH
    p = s/2.0
    yWidths = [bbox(t).width() for t in matrix.yValues]
    yWidth = max(yWidths)
    
    if angle:
      xHeight = max(xWidths)
      xHeight *= 0.7071 # sin(angle)
      xHeight += textH
    else:
      xHeight = textH
    
    x1 = x0 + yWidth + p
    x2 = x1 + n * s
    y1 = y0 + m * s
    y2 = y1 + xHeight + p
    
    y = y1 - s + p - textH/2.0
    
    for i, item in enumerate(self.yAxisItems):
      x = x0 + yWidths[i]/2.0
      item.setPos(x, y)
      y -= s
    
    x = x1 + p
    
    if angle:
      y = y1 + xHeight/2.0
    else: 
      y = y1 + p
      
    for item in self.xAxisItems:
      item.setPos(x, y)
      x += s
    
    rects = []
    rectsAppend = rects.append
    colors = []
    colorsAppend = colors.append
    
    if isinstance(data[0][0], (int, float)):
      x = x1
      for i in range(n):
        row = data[i]
        y = y1-s
 
        for j in range(m):
          z = row[j]
          if z < 0:
            p = max(z,-zMax)/-zMax
            q = 1.0 - p
            r = (p * r3) + (q * r2)
            g = (p * g3) + (q * g2)
            b = (p * b3) + (q * b2)
 
          else:
            p = min(z,zMax)/zMax
            q = 1.0 - p
            r = (p * r1) + (q * r2)
            g = (p * g1) + (q * g2)
            b = (p * b1) + (q * b2)
 
          rectsAppend( QRectF(x, y, s, s) )
          colorsAppend( QColor(r, g, b) )
 
          y -= s
 
        x += s
        
    else:
      x = x1
      l = len(data[0][0])
      lf = float(l)
      
      for i in range(n):
        row = data[i]
        y = y1-s
 
        for j in range(m):
          cell = row[j]
          
          r = 0.0
          g = 0.0
          b = 0.0
          
          for k in range(l):
            p1 = k/(lf-1.0)
            q1 = 1.0 - p1
            
            # Max val colour
            rc = (p1 * r1) + (q1 * r3)
            gc = (p1 * g1) + (q1 * g3)
            bc = (p1 * b1) + (q1 * b3)
            
            z = abs(cell[k])
            p = min(z,zMax)/zMax
            q = 1.0 - p
            
            # Encode value
            r += (p * rc) + (q * r2)
            g += (p * gc) + (q * g2)
            b += (p * bc) + (q * b2)
          
          # Blend  
          r /= lf
          g /= lf
          b /= lf
           
          rectsAppend( QRectF(x, y, s, s) )
          colorsAppend( QColor(r, g, b) )
 
          y -= s
 
        x += s
    
      
    self.legendItem.setMatrix(matrix, s, posColor,
                              zeroColor, negColor)
    self.legendItem.setPos(x2+s, y0)
    self.cursorItem.setPos(0.5*(x0+x1), y2+s)
    self.plotRegion = (x1, y0, x2, y1, 0.5*(x1+x2), 0.5*(y0+y1))
    self.boundingRegion = QRectF(x0, y0, x2-x0, y2-y0) 
    self.drawData = (rects, colors, QRectF(x1, y0, s * n, s * n)) 
    self.update()

  def _highlightOff(self):
  
    if self.cursor:
      i, j, rect = self.cursor
      self.xAxisItems[i].highlight(False)
      self.yAxisItems[j].highlight(False)
      self.cursor = None
      
    
  def hoverMoveEvent(self, event):
    
    x1, y1, x2, y2, cx, cy = self.plotRegion
    pos = event.pos()
    x = pos.x()
    y = pos.y()
    
    if (x1 <= x <= x2) and (y1 <= y <= y2):
      self._highlightOff()
      s = self.boxSize
      dx = (x-x1)
      dy = (y2-y)
      
      n = len(self.xAxisItems)
      m = len(self.yAxisItems)
      
      i = min(int( dx // s ), n-1)
      j = min(int( dy // s ), m-1)
      self.xAxisItems[i].highlight(True)
      self.yAxisItems[j].highlight(True)
      
      xb = x1 + (i * s)
      yb = y2 - (j * s) - s
      self.cursor = (i, j, QRectF(xb, yb, s, s))
      
      matrix = self.matrix
      z = matrix.matrix[i][j]
      xVal = matrix.xValues[i]
      yVal = matrix.yValues[j]
      
      if isinstance(z, (float, int)):
        zStr = '%.2f' % z
      else:
        zStr = ','.join(['%.2f' % v for v in z])
      
      text = '%s, %s: %s' % (xVal, yVal, zStr)
      self.cursorItem.setText(text)
      
      self.update()
    
    else: 
      self._highlightOff()
    
    return QtGui.QGraphicsItem.hoverMoveEvent(self, event)

  def hoverLeaveEvent(self, event):
    
    self._highlightOff()  
    self.cursor = False
    self.update()
    
    
  def boundingRect(self):
    
    if self.boundingRegion:
      return self.boundingRegion # .adjust(-2,-2, 2, 2)
    
    else:
      return NULL_RECT
      
      
  def paint(self, painter, option, widget):
    
    drawRect = painter.drawRect
    setPen   = painter.setPen
    setBrush = painter.setBrush
    
    rects, colors, outline = self.drawData
    
    if self.grid:
      setPen(self.gridPen)
      
      for k, rect in enumerate(rects):
        setBrush(colors[k])
        drawRect(rect)
    
    else:
      for k, rect in enumerate(rects):
        setPen(colors[k])
        setBrush(colors[k])
        drawRect(rect)
      
      setPen(GREY_PEN)
      setBrush(TRANSPARENT)
      drawRect(outline)
   
    if self.cursor:
      i, j, rect = self.cursor
      # ? inverse color ?
      setPen(HIGHLIGHT_PEN)
      setBrush(HIGHLIGHT_BRUSH)
      drawRect(rect)
      
    setPen(Qt.black)
    
    # Render as bitmap? - Far fewer objs.


INT = type(1)
FLOAT = type(1.0)

def binScatterDensity(dataSet, xBins=10, yBins=10):

  pairs = [point[:2] for point in dataSet] 
  matrix = [[0.0] * yBins for i in range(xBins)]

  xVals, yVals = zip(*pairs)
  
  xMin = min(xVals)
  yMin = min(yVals)
  
  xMax = max(xVals)
  yMax = max(yVals)
  
  xStep = float(xMax - xMin) / xBins
  yStep = float(yMax - yMin) / yBins
  
  for x,y in pairs:
  
    i = (x-xMin) // xStep
    j = (y-yMin) // yStep
    
    matrix[int(i)][int(j)] += 1.0

  return matrix


class MatrixDataSet:

  def __init__(self, matrix, xValues=None, yValues=None):

    if type(matrix) is type(dict):
      matrix, xValues, yValues = self._extractDict(matrix)
    
    # TBC: NumPy types
  
    self.xValues = self._checkValues(xValues, len(matrix))
    self.yValues = self._checkValues(yValues, len(matrix[0]))
    
    
    # TBC Qt 2D matrix
    self.matrix = matrix
  
  def getMaxValue(self):
    
    if isinstance(self.matrix[0][0], (int, float)):
      return max([max(row) for row in self.matrix])
    else:
      return max([max(min(cell for cell in row)) for row in self.matrix])
    
  def getMinValue(self):
    
    if isinstance(self.matrix[0][0], (int, float)):
      return min([min(row) for row in self.matrix])
    else:
      return min([min(min(cell for cell in row)) for row in self.matrix])
      
  def getValueLimit(self):
  
    matrix = self.matrix
    
    if isinstance(matrix[0][0], (int, float)):
      minZ = min([min(row) for row in matrix])
      maxZ = max([max(row) for row in matrix])
    
    else:
      minZ = min([min(min(cell for cell in row)) for row in matrix])
      maxZ = max([max(min(cell for cell in row)) for row in matrix])
       
    return max(maxZ, abs(minZ))
    
    
  def _extractDict(self, dataDict):
    
    keys = dataDict.keys()
    
    xValues, yValues = zip(*keys)
    
    xValues = list(set(xValues))
    yValues = list(set(yValues))
    
    n = len(xValues)
    m = len(yValues)
    
    xValues.sort()
    yValues.sort()
    
    matrix = [[0.0] * m for i in range(n)]
    
    for i, x in enumerate(xValues):
      for j, y in enumerate(yValues):
        z = dataDict.get((x,y), 0.0)
        matrix[i][j] = z
    
    return matrix, xValues, yValues
      
      
  def _checkValues(self, values, n):
   
    if values is None:
      values = list(range(1, n+1))
    
    for i, value in enumerate(values):
      if type(value) is INT:
        values[i] = u'%d' % (value,)
      
      elif type(value) is FLOAT:
        if value == 0.0:
          values[i] = u'0'
        
        elif abs(value) > 999999 or abs(value) < 0.01:
          values[i] = u'%5.2e' % (value,)
        
        else:
          values[i] = u'%s' % (value,)
          
      # Otherwise string or unicode is fine
       
    while len(values) < n:
      values.append('')
    
    return values    


class DensityPlot(QtGui.QGraphicsView, Base):

  def __init__(self, parent, matrix, title='Density Plot', valLimit=None,
               xAxisName='X Axis', yAxisName='Y Axis', boxSize=25,
               showGrid=True, showLegend=True, **kw):
      
    QtGui.QGraphicsView.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.parent = parent
    self.title = title
    self.xAxisName = xAxisName
    self.yAxisName = yAxisName
    self.showGrid = showGrid
    self.movePos = None

    self.setRenderHint(QtGui.QPainter.Antialiasing)
    #self.setResizeAnchor(QtGui.QGraphicsView.NoAnchor)
    #self.setResizeAnchor(QtGui.QGraphicsView.AnchorViewCenter)
    self.setResizeAnchor(QtGui.QGraphicsView.AnchorUnderMouse)
    self.setTransformationAnchor(QtGui.QGraphicsView.NoAnchor)
    self.setViewportUpdateMode(QtGui.QGraphicsView.FullViewportUpdate)

    self.gScene = QtGui.QGraphicsScene(self)
    self.setScene(self.gScene)   
    
    self.titleItem  = MovableLabelItem(self)
    self.matrixItem = MatrixItem(self, None, boxSize, valLimit, showGrid)
    self.coordsItem = MovableLabelItem(self)
    self.xLabelItem = MovableLabelItem(self)
    self.yLabelItem = MovableLabelItem(self, angle=-90)
    
    self.setMatrix(matrix)
    x1, y1, x2, y2, cx, cy = self.matrixItem.getPlotCoords()
    self.titleItem.set(cx, y1-boxSize, self.title)
    
    self.contextMenu = self.configMenu()
    self.show()
  
  def setBoxSize(self, size):
    
    self.matrixItem.setBoxSize(size)
    
    x1, y1, x2, y2, cx, cy = self.matrixItem.getPlotCoords()
    x3, y3, x4, y4, dx, dy = self.matrixItem.getBoundingCoords()
    
    # Top center
    self.xLabelItem.set(cx, y4+size, self.xAxisName)
    self.yLabelItem.set(x3-size, cy, self.yAxisName)
    self.titleItem.set(cx, y1-size, self.title)
  
  def setGrid(self, action):
  
    self.matrixItem.setGrid(action.isChecked())

  def setLegend(self, action):
  
    self.matrixItem.setLegend(action.isChecked())

  def setCoords(self, action):
  
    self.matrixItem.setCoords(action.isChecked())

  def setValueLimit(self, val):
  
    self.matrixItem.setValueLimit(val)
  
  def setMatrix(self, matrix, xAxisName=None, yAxisName=None):
    
    if xAxisName is not None:
      self.xAxisName = xAxisName
      
    if yAxisName is not None:
      self.yAxisName = yAxisName
    
    self.matrixItem.setMatrix(matrix)
    
    x1, y1, x2, y2, cx, cy = self.matrixItem.getPlotCoords()
    x3, y3, x4, y4, dx, dy = self.matrixItem.getBoundingCoords()
    
    # Top center
    s = self.matrixItem.boxSize
    self.xLabelItem.set(cx, y4+s, self.xAxisName)
    self.yLabelItem.set(x3-s, cy, self.yAxisName)
    
    
  def setColors(self, posColor, negColor, zeroColor):
    
    if not isinstance(posColor, QColor):
      posColor = QColor(posColor)
    
    if not isinstance(negColor, QColor):
      negColor = QColor(negColor)
    
    if not isinstance(zeroColor, QColor):
      zeroColor = QColor(zeroColor)
    
    self.matrixItem.setColors(posColor, negColor, zeroColor)
   

  def setTitle(self, title=''):
    
    self.title = title
    self.titleItem.setText(title)
    
    self.update()

  
  def askGraphTitle(self):
  
    text = askString('Text Entry','Enter plot title', self.title, self)
    self.setTitle(text)
  
  def setPosColor(self):
  
    p, z, n = self.matrixItem.colors
    
    dialog = ColorDialog()
    p = dialog.getColor() # opens the dialog
    
    if p:
      self.matrixItem.setColors(p, z, n)
      self.posColAction.setIcon(Icon(color=p))
  
  def setZeroColor(self):
  
    p, z, n = self.matrixItem.colors
    
    dialog = ColorDialog()
    z = dialog.getColor() # opens the dialog
    
    if z:
      self.matrixItem.setColors(p, z, n)
      self.zeroColAction.setIcon(Icon(color=z))
  
  
  def setNegColor(self):
  
    p, z, n = self.matrixItem.colors
    
    dialog = ColorDialog()
    n = dialog.getColor() # opens the dialog
    
    if n:
      self.matrixItem.setColors(p, z, n)
      self.negColAction.setIcon(Icon(color=n))
  
  def wheelEvent(self, event):
    
    if event.delta() < 0:
      fac = 0.8333
    else:
      fac = 1.2
    
    self.scale(fac, fac)
  
    event.accept()
    
    #QtGui.QGraphicsView.wheelEvent(self, event) 
  
  def mouseReleaseEvent(self, event):
    
    self.movePos = None
    
    return QtGui.QGraphicsView.mouseReleaseEvent(self, event)
        
  def mousePressEvent(self, event):
   
    button = event.button()
    
    # deal with inconsistency in Qt versions for button naming
    try:
      MiddleButton = Qt.MiddleButton
    except AttributeError:
      MiddleButton = Qt.MidButton

    if button == Qt.RightButton:
      self.popupContextMenu(event.globalPos())
    
    elif button == MiddleButton:
      pos = event.pos()
      h = self.horizontalScrollBar().sliderPosition()
      v = self.verticalScrollBar().sliderPosition()
      self.movePos = pos.x()+h, pos.y()+v
    
    return QtGui.QGraphicsView.mousePressEvent(self, event)
    
  def mouseMoveEvent(self, event):
      
    button = event.button()
    
    if self.movePos:
      x0, y0 = self.movePos
      pos = event.pos()
      self.horizontalScrollBar().setSliderPosition(x0-pos.x())
      self.verticalScrollBar().setSliderPosition(y0-pos.y())
    
    return QtGui.QGraphicsView.mouseMoveEvent(self, event)
   
  def popupContextMenu(self, pos):
    
    self.contextMenu.popup(pos)   
    
  def resetView(self):
  
    self.fitInView(self.scene().itemsBoundingRect(), Qt.KeepAspectRatio)
  
  def configMenu(self):

    # Context menu
    # - As squares of changing size?
    
    p, z, n = self.matrixItem.colors
    
    mouseMenu = Menu(None, 'Plot Options')
    mouseMenu.addItem('Show Grid?', self.setGrid, checked=True)
    mouseMenu.addItem('Show Legend?', self.setLegend, checked=True)
    mouseMenu.addItem('Show Cursor Coords?', self.setCoords, checked=True)
    sizeMenu = Menu(mouseMenu, 'Box Size')
    for sz in range(5,55,5):
      sizeMenu.addItem('%d' % sz,lambda x=sz: self.setBoxSize(x))

    colorMenu = Menu(mouseMenu, 'Colours')
   
    self.posColAction = colorMenu.addItem('Set positive color',
                                          self.setPosColor,
                                          icon=Icon(color=p))
    self.zeroColAtion = colorMenu.addItem('Set zero color',     
                                          self.setZeroColor,
                                          icon=Icon(color=z))
    self.negColAction = colorMenu.addItem('Set negative color',
                                          self.setNegColor,
                                          icon=Icon(color=n))

    imageMenu = Menu(mouseMenu, 'Export Image')
    
    imageMenu.addItem('JPEG', self.exportJpeg)
    imageMenu.addItem('PNG', self.exportPng)
    imageMenu.addItem('PDF', self.exportPdf)
    imageMenu.addItem('SVG', self.exportSvg)
        
    mouseMenu.addItem('Print', self.printWidget)
 
    mouseMenu.addItem('Set Title', self.askGraphTitle)
   
    mouseMenu.addItem('Reset View', self.resetView)
    
    
    return mouseMenu
  
  def printWidget(self):
  
    dialog = PrintDialog(self) 
    dialog.printWidget(self)
    

if __name__ == '__main__':

  from math import exp

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup

  app = Application()
  
  popup = BasePopup(title='Density Plot Example')
  popup.setSize(500, 400)
  
  xx = [x/10.0 for x in range(1,11)]
  yy = [x*x for x in xx]
  #yy = [x for x in xx]

  data = []
  for x in xx:
    data.append([])
    for y in yy:
      data[-1].append(x*y - 0.4)
      #data[-1].append([y, x])
  
  #xx = [x * 1e6 for x in xx]
  
  matrix = MatrixDataSet(data, xx, yy)

  plot = DensityPlot(popup, matrix)

  app.start()

