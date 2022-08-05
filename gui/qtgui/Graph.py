import sys

from math import sin, cos, ceil

from PySide import QtGui, QtCore
Qt = QtCore.Qt
QPointF = QtCore.QPointF
QRectF = QtCore.QRectF
QColor = QtGui.QColor

from gui.qtgui.Base import Base, Icon
from gui.qtgui.Menu import Menu

from gui.qtgui.InputDialog import askString
from gui.qtgui.graphicsItems.LabelItems import MovableLabelItem, HIGHLIGHT
from gui.qtgui.FileSelect import selectSaveFile, FileType
from gui.qtgui.Print import PrintDialog
from gui.qtgui.Colors import ColorDialog
from gui.qtgui.Entry import IntEntry

LINE = 'line'
SCATTER = 'scatter'
HISTOGRAM = 'histogram'

CIRCLE = 'circle'
CROSS = 'cross'
SQUARE = 'square'
TRIANGLE = 'triangle'
STAR = 'star'
HEXAGON = 'hexagon'

NULL_RECT = QRectF()
NULL_POINT = QPointF()

QApp = QtGui.QApplication.instance

PI = 3.1415926535898
INF = float('inf')

CROSSHAIR_COLOR = QtGui.QColor( 32, 32, 255, 192)

ItemIsMovable = QtGui.QGraphicsItem.ItemIsMovable
ItemPositionChange = QtGui.QGraphicsItem.ItemPositionChange
#ItemIsSelectable = QtGui.QGraphicsItem.ItemIsSelectable
#ItemSendsGeometryChanges = QtGui.QGraphicsItem.ItemSendsGeometryChanges

# # # # #  T O  D O  # # # # #
#
# Menu set symbol, type, symSize & color globally
# Error scalable circles


class DataSetItem(QtGui.QGraphicsItem):
        
  def __init__(self, parent, dataSet=None):
    
    QtGui.QGraphicsItem.__init__(self)
    self.setZValue(1)
    self.parent = parent
    self.drawData = None
    self.region = NULL_RECT
    self.dataSet = dataSet
    
    self.parent.scene().addItem(self)
  
  def refresh(self):
    
    self.setDataSet(self.dataSet)
      
  def setDataSet(self, dataSet):
    
    self.dataSet = dataSet
    parent = self.parent
  
    pr0, pr1, pr2, pr3 = parent.getPlotRegion()
    
    xAxis = parent.xAxis
    
    if dataSet.secondAxis:
      yAxis = parent.yAxis2
    else:
      yAxis = parent.yAxis
       
    xMin, xMax = xAxis.valRange
    yMin, yMax = yAxis.valRange
    
    if xAxis.isReverse:
      xMin, xMax = xMax, xMin
    
    if yAxis.isReverse:
      yMin, yMax = yMax, yMin
    
    xSize = xMax-xMin
    ySize = yMax-yMin
    
    plotWidth = pr2 - pr0
    plotHeight = pr3 - pr1

    ppvX = plotWidth/float(xSize or 1.0)
    ppvY = plotHeight/float(ySize or 1.0)
 
    if yMin < 0.0 < yMax:
      f = yMax / float(yMax-yMin)
      yBase = pr1 + (f * plotHeight)
 
    elif yMax < 0.0:
      yBase = pr1
 
    else: # yMin < 0.0
      yBase = pr3
    
    dataPoints = dataSet.dataPoints
    
    
    xVals = [dp[0] for dp in dataPoints]
    xVals.sort()
    xDeltas = [xVals[i+1]-x for i,x in enumerate(xVals[:-1])]
    
    if xDeltas:
      binWidth = ppvX * 0.4 * min(xDeltas)
    else:
      binWidth = 1.0
    
    coords = []
    coordsAppend = coords.append
    
    for x, y, e in dataPoints:
      # TBC: used to be a check for x, y Nones in prev version
      
      x0 = pr0 + (x - xMin) * ppvX
      y0 = pr3 - (y - yMin) * ppvY
      
       
      if int(y0+1) < pr1:
        continue
        
      if int(y0) > pr3:
        continue
      
      if e is None:
        pointU = None
        pointL = None
        
      else:
        e0 = e * ppvY
        e1 = y0 + e0
        e2 = y0 - e0
        pointU = QPointF(x0,e1)
        pointL = QPointF(x0,e2)
      
      
      point = QPointF(x0, y0)
      yZero = QPointF(x0, yBase)
      coordsAppend( (point, pointU, pointL, yZero) )

    self.dataSet = dataSet
    self.drawData = coords, binWidth
    self.region = QRectF(QPointF(pr0, pr1), QPointF(pr2, pr3))
    
    self.show()
    self.update()

  def boundingRect(self):
    
    if self.region:
      pad = 2
      return self.region.adjusted(-pad, -pad, pad, pad)
    
    else:
      return NULL_RECT
      
  def paint(self, painter, option, widget):
    
    if not self.drawData:
      return
    
    coords, binWidth = self.drawData
    dataSet = self.dataSet
    
    plotType = dataSet.plotType
    symbol = dataSet.symbol
    symSize = dataSet.symbolSize
    
    qColor = QtGui.QColor(dataSet.color)
    qPolyF = QtGui.QPolygonF
    pen = QtGui.QPen(qColor, 0.8, Qt.SolidLine)

    drawEllipse = painter.drawEllipse
    drawLine = painter.drawLine
    drawPoly = painter.drawPolygon
    drawRect = painter.drawRect
    
    n = len(coords)
    r = symSize/2.0
    sx = QPointF(r, 0)
    s2 = QPointF(r, r)
    ignoreZero = dataSet.ignoreZero
        
    # Error margin
    
    if dataSet.errorColors:
      qColor1, qColor2 = dataSet.errorColors
      for i in range(n-1):
        point1, pointU1, pointL1, yZero1 = coords[i]
        point2, pointU2, pointL2, yZero2 = coords[i+1]
        
        if ignoreZero and (point1 == yZero1) and (point2 == yZero2):
          continue
        
        painter.setPen(qColor1)
        painter.setBrush(qColor1)
        if pointU1 and pointU2:
          poly = qPolyF()
          poly.append(point1)
          poly.append(pointU1)
          poly.append(pointU2)
          poly.append(point2)
          drawPoly(poly)
          
        painter.setPen(qColor2)
        painter.setBrush(qColor2)
        if pointL1 and pointL2:
          poly = qPolyF()
          poly.append(point1)
          poly.append(pointL1)
          poly.append(pointL2)
          poly.append(point2)
          drawPoly(poly)
          
    painter.setPen(pen)
    painter.setBrush(qColor)
    
    # Main plot
     
    if plotType == LINE:
      for i in range(n-1):
        point1, pointU1, pointL1, yZero1 = coords[i]
        point2, pointU2, pointL2, yZero2 = coords[i+1]
        if ignoreZero and (point1 == yZero1) and (point2 == yZero2):
          continue
          
        drawLine(point1, point2)
         
    if plotType in (LINE, SCATTER):
      if symbol == TRIANGLE:
        angles = [i*PI/1.5 for i in range(3)]
        trig = [QPointF(r*sin(a),-r*cos(a)) for a in angles]
        
        for point, u, l, z in coords:        
          if ignoreZero and point == z:
            continue
          
          poly = qPolyF()
          
          for delta in trig:
            poly.append(point+delta)

          drawPoly(poly)
 
      elif symbol == HEXAGON:
        angles = [i*PI/3.0 for i in range(6)]
        trig = [QPointF(r*sin(a),-r*cos(a)) for a in angles]
        
        for point, u, l, z in coords:
          if ignoreZero and point == z:
            continue
            
          poly = qPolyF()
          for delta in trig:
            poly.append(point+delta)
          
          drawPoly(poly)
     
      elif symbol == STAR:
        # Pentagram order
        angles = [0.0, 0.8*PI, 1.6*PI, 0.4*PI, 1.2*PI] 
        trig = [QPointF(r*sin(a),-r*cos(a)) for a in angles]
        
        for point, u, l, z in coords:
          if ignoreZero and point == z:
            continue
            
          poly = qPolyF()
          for delta in trig:
            poly.append(point+delta)
          
          drawPoly(poly)
 
      elif symbol == SQUARE:
        for point, pointU, pointL, yZero in coords:
          if ignoreZero and point == yZero:
            continue
            
          drawRect(QRectF(point-s2, point+s2))

      elif symbol == CROSS:
        s3 = QPointF(r, -r)
        for point, pointU, pointL, yZero in coords:
          if ignoreZero and point == yZero:
            continue
            
          drawLine(point-s2, point+s2)
          drawLine(point-s3, point+s3)
 
      else: # CIRCLE
        for point, u, l, z in coords:
          if ignoreZero and point == z:
            continue
            
          drawEllipse(point, r, r)
               
    else: # Histogram
      pen = QtGui.QPen(qColor.darker(), 0.01, Qt.SolidLine)
      painter.setPen(pen)
      pr0, pr1, pr2, pr3 = self.parent.getPlotRegion()

      if dataSet.gradient:
        dataSet.gradient.setStart(0, pr1)
        dataSet.gradient.setFinalStop(0, pr3)
        painter.setBrush(QtGui.QBrush(dataSet.gradient))
      else:
        painter.setBrush(qColor)

      sb = QPointF(binWidth, 0)
      
      for point, pointU, pointL, yZero in coords:
        drawRect(QRectF(point-sb, yZero+sb))
     
    # Error Bars
    
    if not dataSet.errorColors:
      painter.setPen(pen)
      for point, pointU, pointL, z in coords:
        if ignoreZero and point == z:
          continue
          
        if pointU and pointL:
          drawLine(pointU, pointL)
          drawLine(pointU-sx, pointU+sx)
          drawLine(pointL-sx, pointL+sx)

class AxisItem(QtGui.QGraphicsItem):
        
  def __init__(self, parent, axis=None, isVertical=False, secondAxis=False):
    
    QtGui.QGraphicsItem.__init__(self)
    self.setZValue(1)
    self.parent = parent
    self.drawData = None
    self.gridData = None
    self.region = NULL_RECT
    self.isVertical = isVertical
    self.secondAxis = secondAxis
    self.axis = axis
    self.bgColor = parent.bgColor
    self.mgColor = parent.mgColor
    self.fgColor = parent.fgColor
    
    if isVertical:
      if secondAxis:
        angle = 90
      else:
        angle = -90
    else:
      angle = None
      
    self.parent.scene().addItem(self)
    self.label = MovableLabelItem(parent, angle, parent.fgColor)

  def refresh(self):
    
    self.setAxis(self.axis)      
      
  def setAxis(self, axis):
    
    if axis is not self.axis:
      self.axis = axis
    
    if axis.isReverse:
      dataMax, dataMin = axis.valRange
    else:
      dataMin, dataMax = axis.valRange
    
    dataWidth = (dataMax - dataMin) or 1.0
    
    box = self.parent.getPlotRegion()
    
    if self.isVertical:
      plotMin = box[1]
      plotMax = box[3]
      
      if self.secondAxis:
        rx0, ry0, rx1, ry1 = box[2], plotMax, box[2], plotMin
      else:
        rx0, ry0, rx1, ry1 = box[0], plotMax, box[0], plotMin
      
    else:  
      plotMin = box[0]
      plotMax = box[2]
      rx0, ry0, rx1, ry1 = plotMin, box[3], plotMax, box[3]
    
    plotWidth = plotMax - plotMin
  
    tickLen = 8.0
    labels = axis.labels
    ticks = []
    grid = []
     
    if axis.ticks:
      ppv = plotWidth/float(dataWidth) or 1.0
 
      # Closest two labels can get
      fontMetric = QtGui.QFontMetricsF(self.parent.font())
      
      if self.isVertical:
        minPlotSep = 3.0 * fontMetric.boundingRect('X').height()
      else:
        sepA = fontMetric.boundingRect(self.formatText(dataMax)).width()
        sepB = fontMetric.boundingRect(self.formatText(dataMin)).width()
        minPlotSep = 2.0 * max(sepA, sepB)
       
      minDataSep = minPlotSep / ppv
     
      sci = '%e' % abs(minDataSep)
      deci = int(sci[-3:]) # num decimal places
      sigD = int(sci[0]) # significant digit

      num = 10.0
      sX = abs(sigD-num)
      for n in (1.0, 2.0, 5.0):

        s = abs(sigD-n)
        if s < sX:
          sX = s
          num = n 
 
      # Check y not in opp dir (-num)
      
      if labels and len(labels) > 1:
        dataStep = dataWidth/(len(labels)-1)
      else:
        dataStep = (abs(dataWidth)/dataWidth) *  num * 10**(deci)
      
      dataVal = dataMin - (dataMin % dataStep)
      
      getBbox = fontMetric.tightBoundingRect
      
      for i in range(int(ceil(dataWidth/dataStep))+2):
        
        roundedVal = round(dataVal,-deci) # avoid FP errors
        
        if self.isVertical:
          plotVal = plotMax - (dataVal - dataMin) * ppv
        else:
          plotVal = plotMin + (dataVal - dataMin) * ppv
        
        dataVal += dataStep
        
        if not (plotMin <= plotVal <= plotMax):
          continue

        if labels:
          if i >= len(labels):
            continue
        
          text  = labels[i]
 
        else:
          text = self.formatText(roundedVal)
        
        bbox = getBbox(text)
        w = bbox.width()
        h = bbox.height()
 
        if self.isVertical:
          
          if self.secondAxis:
            x0 = box[2] + tickLen
            x1 = box[2]
            x2 = box[0]
            left = x0 + tickLen/2.0
            right = x0 + w
            
            if right > rx1:
              rx1 = right
            
          else:
            x0 = box[0] - tickLen
            x1 = box[0]
            x2 = box[2]
            left = x0 - w - tickLen/2.0
 
            if left < rx0:
              rx0 = left
          
          y0 = y1 = y2 = plotVal
          bottom = y0 + h/2.0
          top = y0+h/2
          
          if bottom > ry1:
            ry1 = bottom
          
          if top < ry0:
            ry0 = top
          
          textPoint = QPointF(left, bottom) # Check Y direction
 
        else: # Check y dir
          y0 = box[3] + tickLen
          y1 = box[3]
          y2 = box[1]
          x0 = x1 = x2 = plotVal
          left = x0 - w/2.0
          right = left + w
          bottom = y0 + h + tickLen/2.0
          
          if left < rx0:
            rx0 = left
          
          if bottom > ry1:
            ry1 = bottom
          
          if right > rx1:
            rx1 = right 
          
          textPoint = QPointF(left, bottom) # Check Y direction

        start = QPointF(x0, y0)
        end = QPointF(x1, y1)
        
        if i > 0:
          gridStart = QPointF(x1, y1)
          gridEnd = QPointF(x2, y2)
          grid.append( (gridStart, gridEnd) )
        
        ticks.append( (start, end, text, textPoint) )

    
    if self.isVertical:
      if self.secondAxis:
        rx1 += tickLen
        region = QRectF(QPointF(rx0, ry0), QPointF(rx1, ry1))
        self.label.set(rx1+tickLen, region.center().y(), axis.name)
        startPoint = QPointF(box[2], plotMin)
        endPoint = QPointF(box[2], plotMax)
      
      else:
        rx0 -= tickLen
        region = QRectF(QPointF(rx0, ry0), QPointF(box[2], ry1))
        self.label.set(rx0-tickLen, region.center().y(), axis.name)
        startPoint = QPointF(box[0], plotMin)
        endPoint = QPointF(box[0], plotMax)
 
    
    else:
      ry1 += tickLen
      region = QRectF(QPointF(rx0, box[1]), QPointF(rx1, ry1))
      startPoint = QPointF(box[2], box[3])
      endPoint = QPointF(box[0], box[3])
      self.label.set(region.center().x(), ry1+tickLen, axis.name)
      
    self.drawData = startPoint, endPoint, ticks
    self.gridData = grid
    self.region = region
    self.update()
  
  def formatText(self,text):
  
    if isinstance(text, float):
      if text == 0:
        text = '0'
      elif abs(text) > 999999 or abs(text) < 0.01:
        text = '%3.1e' % text
      else:
        text = str(text)
    
    elif isinstance(text, int):
      text = str(text)

    elif isinstance(text, bool):
      text = str(text and 1 or 0)
      
    if text and text[0:1] == '@':
      text = ''

    return text
    
  def boundingRect(self):
    
    # Need to account for grid?
    
    if self.region:
      pad = 2
      return self.region.adjusted(-pad, -pad, pad, pad)
    
    else:
      return NULL_RECT
        
  def paint(self, painter, option, widget):
    
    if self.drawData:
      startPoint, endPoint, ticks = self.drawData
      axis = self.axis
 
      painter.setPen(self.mgColor)
      if axis.grid:
        for start, end in self.gridData:
          painter.drawLine(start, end)

      painter.setPen(self.fgColor) # Could set width
      painter.drawLine(startPoint, endPoint)
 
      painter.setPen(self.fgColor)
      for start, end, text, point in ticks:
        painter.drawLine(start, end)
        painter.drawText(point, text)

      
      
class CrosshairItem(QtGui.QGraphicsItem):
        
  def __init__(self, parent):
    
    QtGui.QGraphicsItem.__init__(self)
    self.setZValue(2)
    self.parent = parent
    self.drawData = (NULL_POINT, NULL_POINT)
    self.region = NULL_RECT
    self.parent.scene().addItem(self)
    
  def set(self, x0, y0, x1, y1):
    
    begin = QPointF(x0, y0)
    end = QPointF(x1, y1)
    self.drawData = (begin, end)
    self.region = QRectF(begin, end) 
    
    self.show()
    self.update()
      
  def boundingRect(self):
    
    if self.region:
      return self.region
    
    else:
      return NULL_RECT
      
  def paint(self, painter, option, widget):
   
    startPoint, endPoint = self.drawData
    
    #pen = QtGui.QPen(CROSSHAIR_COLOR, 2, Qt.DotLine)
    painter.setPen(CROSSHAIR_COLOR)

    painter.drawLine(startPoint, endPoint)
    
class CoordsItem(QtGui.QGraphicsItem):
        
  def __init__(self, parent):
    
    QtGui.QGraphicsItem.__init__(self)
    self.setZValue(2)
    self.parent = parent
    self.textBox = NULL_RECT
    self.region = NULL_RECT
    self.text = ''
    self.fgColor = parent.fgColor
    self.bgColor = parent.bgColor
    self.mgColor = parent.mgColor
   
    pr0, pr1, pr2, pr3 = parent.getPlotRegion()
    self.setPos(pr2-50, pr1-10)
    
    #self.setFlag(ItemIsSelectable)
    #self.selected = False
    #self.setFlag(ItemSendsGeometryChanges)
          
    self.hover = False
    self.setFlag(ItemIsMovable)
    self.setAcceptHoverEvents(True)
    self.setAcceptedMouseButtons(Qt.LeftButton)
    self.parent.scene().addItem(self)

  def hoverEnterEvent(self, event):
  
    self.hover = True
    self.update()

  def hoverLeaveEvent(self, event):
  
    self.hover = False
    self.update()
  
  def itemChange(self, change, value):
    
    #if change == ItemPositionChange:
    #  x = value.x()
    #  y = value.y()
    #  # self.update()
      
    return QtGui.QGraphicsItem.itemChange(self, change, value)
    
  def set(self, vals=None):
    
    if vals:
      self.text = text = ', '.join(['%.3f' % x for x in vals])
    else:
      self.text = text = ''
    
    fontMetric = QtGui.QFontMetricsF(self.parent.font())
    
    pad = 2
    self.textBox = box = fontMetric.tightBoundingRect(text)
    self.region = box.adjusted(-pad, -pad, pad, pad)
    
    self.update()
      
  def boundingRect(self):
    
    if self.region:
      pad = 2
      return self.region.adjusted(-pad, -pad, pad, pad)
    
    else:
      return NULL_RECT
      
  def paint(self, painter, option, widget):
    
    if not self.text:
      return
    
    rect = self.region
    
    if self.hover:
      painter.setPen(HIGHLIGHT)
      painter.setBrush(self.mgColor)
    
    else:
      painter.setPen(self.mgColor)
      painter.setBrush(self.bgColor)
      
    painter.drawRect( self.region )
    
    painter.setPen(self.fgColor)
    painter.drawText(self.textBox.bottomLeft(), self.text)
    
class LegendItem(QtGui.QGraphicsItem):
        
  def __init__(self, parent, dataSets=None):
    
    QtGui.QGraphicsItem.__init__(self)
    self.setZValue(2)
    self.parent = parent
    self.region = NULL_RECT
    self.drawData = None
    self.fgColor = parent.fgColor
    self.bgColor = parent.bgColor
    self.mgColor = parent.mgColor
    
    pr0, pr1, pr2, pr3 = parent.getPlotRegion()
    self.setPos(pr2-50, pr1+5)
   
    if dataSets:
      self.set(dataSets)
    
    self.hover = False
    self.setFlag(ItemIsMovable)
    self.setAcceptHoverEvents(True)
    self.setAcceptedMouseButtons(Qt.LeftButton)
    self.parent.scene().addItem(self)

  def hoverEnterEvent(self, event):
  
    self.hover = True
    self.update()

  def hoverLeaveEvent(self, event):
  
    self.hover = False
    self.update()
  
  def itemChange(self, change, value):
    
    #if change == ItemPositionChange:
    #  x = value.x()
    #  y = value.y()
    #  # self.update()
      
    return QtGui.QGraphicsItem.itemChange(self, change, value)
    
  def set(self, dataSets):
    
    # Given this is movable all coords are relative to zero
    
    if dataSets:
      fontMetric = QtGui.QFontMetricsF(self.parent.font())
      
      data = [(QColor(ds.color), ds.name, ds.symbol, ds.gradient) for ds in dataSets]
      colors, texts, symbols, grads = zip(*data)
      
      symbols = list(symbols)
      tBoxes = [fontMetric.tightBoundingRect(t) for t in texts]
      widths = [rect.width() for rect in tBoxes]
      heights = [rect.height() for rect in tBoxes]
      
      size = max(heights)
      width = max(widths)
      
      pos = self.pos()
      r = size/2.0
      x0 = 0
      y0 = 0
       
      sPoints = []
      points = []
      
      pad = 4
      y = y0 + pad
      for i in range(len(texts)):
        
        if dataSets[i].plotType == HISTOGRAM:
          symbols[i] = SQUARE
          
        x = x0 + pad
        
        sPoints.append( QPointF(x+r, y+r) )
        
        x += 1.5 * size
        points.append( QPointF(x, y+heights[i]/2.0+r) )
        
        y += size + pad
      
      brushes = []
      for grad in grads:
        if grad:
          grad.setStart(0, 0.0 * size)
          grad.setFinalStop(0, 1.0 * size)
          brushes.append(QtGui.QBrush(grad))
        else:
          brushes.append(None)
      
      self.drawData = r, symbols, sPoints, colors, points, texts, brushes
      self.region = QRectF(x0, y0, width+2*size, y-y0)
 
      self.show()
      self.update()
    
    else:
      self.drawData = None
      self.hide()
      
  def boundingRect(self):
    
    if self.region:
      pad = 2
      return self.region.adjusted(-pad, -pad, pad, pad)
    
    else:
      return NULL_RECT
      
  def paint(self, painter, option, widget):
     
    if self.drawData:
      rect = self.region
      
      if self.hover:
        painter.setPen(HIGHLIGHT)
        painter.setBrush(self.mgColor)
      
      else:
        painter.setPen(self.mgColor)
        painter.setBrush(self.bgColor)
      
      drawPoly = painter.drawPolygon
      drawRect = painter.drawRect
      drawEllipse = painter.drawEllipse
      drawLine = painter.drawLine
      
      drawRect(rect)
      qPolyF = QtGui.QPolygonF
      
      r, symbols, sPoints, colors, points, texts, brushes = self.drawData
      
      for i, point in enumerate(points):
        symbol = symbols[i]
        sPoint = sPoints[i]
        color = colors[i]
        
        painter.setPen(self.fgColor)
        painter.drawText(point, texts[i])
        painter.setBrush(brushes[i] or color)
        painter.setPen(color.darker())
        
        if symbol == TRIANGLE:
          angles = [i*PI/1.5 for i in range(3)]
          trig = [QPointF(r*sin(a),-r*cos(a)) for a in angles]
          poly = qPolyF()
 
          for delta in trig:
            poly.append(sPoint+delta)

          drawPoly(poly)
 
        elif symbol == HEXAGON:
          angles = [i*PI/3.0 for i in range(6)]
          trig = [QPointF(r*sin(a),-r*cos(a)) for a in angles]
          poly = qPolyF()
          for delta in trig:
            poly.append(sPoint+delta)
 
          drawPoly(poly)
 
        elif symbol == STAR:
          # Pentagram order
          angles = [0.0, 0.8*PI, 1.6*PI, 0.4*PI, 1.2*PI]
          trig = [QPointF(r*sin(a),-r*cos(a)) for a in angles]
          poly = qPolyF()
          for delta in trig:
            poly.append(sPoint+delta)
 
          drawPoly(poly)
 
        elif symbol == SQUARE:
          s2 = QPointF(r, r)
          drawRect(QRectF(sPoint-s2, sPoint+s2))

        elif symbol == CROSS:
          s2 = QPointF(r, r)
          s3 = QPointF(r, -r)
          drawLine(sPoint-s2, sPoint+s2)
          drawLine(sPoint-s3, sPoint+s3)
 
        else: # CIRCLE
          drawEllipse(sPoint, r, r)
        
class GraphDataSet:

  def __init__(self, dataPoints=None, name=None, color=None,
               plotType=LINE, lineWidth=1.0, symbol=None,
               symbolSize=5.0, secondAxis=False, gradient=None,
               errorColors=None, ignoreZero=False): 
    
    # List of (x,y) or (x,y,dy)
    if dataPoints and (len(dataPoints[0]) == 2):
      yErrors = [None] * len(dataPoints)
      xVals, yVals = list(zip(*dataPoints))
      dataPoints = list(zip(xVals, yVals, yErrors))
      
    self.dataPoints = dataPoints or []
    
    self.name = name
    self.color = QtGui.QColor(color)
    self.plotType = plotType
    self.lineWidth = lineWidth
    self.symbol = symbol
    self.symbolSize = symbolSize
    self.xDataRange = (0.0, 1.0)
    self.yDataRange = (0.0, 1.0)
    self.secondAxis = secondAxis
    self.gradient = gradient
    self.errorColors = errorColors
    self.ignoreZero = ignoreZero

    self.setDataRanges()
  
  
  def setDataRanges(self):
    
    dataPoints = self.dataPoints
    dataPoints.sort()

    xMin = dataPoints[0][0]
    xMax = dataPoints[-1][0]
    
    yUpper = [y + (dy or 0.0) for x, y, dy in dataPoints]
    yLower = [y - (dy or 0.0) for x, y, dy in dataPoints]
    
    yMin = min(yLower)
    yMax = max(yUpper)
    
    self.xDataRange = (xMin, xMax)
    self.yDataRange = (yMin, yMax)    
    
    
  def setDataLists(self, xList, yList, yErrors=None):
  
    n = len(xList)
    assert(len(yList) == n)
    
    if not yErrors:
      yErrors = [None] * n
    else:
      assert(len(yErrors) == n)
    
    self.dataPoints = list(zip(xList, yList, yErrors))
    self.setDataRanges()
        

class GraphAxis:

  def __init__(self, name=None, labels=None, ticks=True, valRange=None,
               grid=True, isLog=False, isReverse=False):
    
    # ? Font ? 
    self.name = name
    self.labels = labels
    self.ticks = ticks
    self.valRange = valRange or (0.0, 1e-10)
    self.grid = grid
    self.isLog = isLog
    self.isReverse = isReverse
    
DEFAULT_COLORS = ['#800000','#000080',
                  '#008000','#808000',
                  '#800080','#008080',
                  '#808080','#000000',
                  '#804000','#004080']
                 
class Graph(QtGui.QGraphicsView, Base):
  
  linePlot = LINE
  scatterPlot = SCATTER
  histogramPlot = HISTOGRAM
  
  # TBC CanvasOrigin equiv; for making multi-plots

  def __init__(self, parent, axes, dataSets=None, title=None, size=(400, 300), zoom=1.0,
               showCoords=True, callback=None, motionCallback=None, showLegend=True,
               bgColor=QtGui.QColor(255, 255, 255, 255),
               fgColor=QtGui.QColor(0, 0, 0, 255),
               mgColor=QtGui.QColor(200, 200, 200, 255), **kw):
    
    QtGui.QGraphicsView.__init__(self, parent)
    Base.__init__(self, parent, **kw)
 
    self.axes = axes
    if len(axes) == 2:
      self.xAxis, self.yAxis = axes
      self.yAxis2 = None
    else:
      self.xAxis, self.yAxis, self.yAxis2 = axes
    
    self.bgColor = bgColor
    self.fgColor = fgColor
    self.mgColor = mgColor
    self.title = title or ''
    self.size = size
    self.zoom = zoom or 1.0
    self.showCoords = showCoords
    self.callback = callback
    self.motionCallback = motionCallback
    self.movePos = None
    self.showLegend = showLegend
 
    if not self.xAxis.name:
      self.xAxis.name = 'X axis'
    if not self.yAxis.name:
      self.yAxis.name = 'Y axis'
    if self.yAxis2 and not self.yAxis2.name:
      self.yAxis2.name = 'Y axis 2'
    
    self.setBackgroundBrush(QtGui.QBrush(bgColor)) 
    self.setAutoFillBackground(True)
      
    self.setMinimumSize(*size)
    self.setRenderHint(QtGui.QPainter.Antialiasing)
    #self.setCacheMode(QtGui.QGraphicsView.CacheBackground)
    
    self.setResizeAnchor(self.AnchorViewCenter)
    self.setViewportUpdateMode(QtGui.QGraphicsView.FullViewportUpdate)
    #self.setDragMode(QtGui.QGraphicsView.ScrollHandDrag)
    #self.setInteractive(True)

    self.graphicsScene = QtGui.QGraphicsScene(self)
    self.setScene(self.graphicsScene)   
    
    self.setTransformationAnchor(self.AnchorViewCenter)
    self.setSceneRect(self.viewport().rect())
    
    self.contextMenu = self.configMenu()
    
    self.titleItem = MovableLabelItem(self, pen=fgColor)
    self.legendItem = LegendItem(self)
    self.coordsItem = CoordsItem(self)
    self.xCrosshair = CrosshairItem(self)
    self.yCrosshair = CrosshairItem(self)
    self.xAxisItem = AxisItem(self, self.xAxis)
    self.yAxisItem = AxisItem(self, self.yAxis, isVertical=True)
    if self.yAxis2:
      self.yAxisItem2 = AxisItem(self, self.yAxis2, isVertical=True, secondAxis=True)
    
    self.indicatorLinesH = []
    self.indicatorLinesV = []
    self.dataSetItems = []

    dataSets = dataSets or []
    self.dataSets = []
    for dataSet in dataSets:
      self.addDataSet(dataSet, False)
      
    x0, y0, x1, y1 = self.getPlotRegion()
    self.titleItem.set(0.5 * (x0+x1), y0-(self.zoom*25), self.title)
    self.refresh()
    
    self.show()
    
    if not showLegend:
      self.legendItem.hide()
  
  def sizeHint(self):
  
    w, h = self.size
    return QtCore.QSize(w+200, h+100)
  
    
  def setGrid(self, action):
    
    bool = action.isChecked()
    self.xAxis.grid = bool
    self.xAxisItem.update()
    self.yAxis.grid = bool
    self.yAxisItem.update()

  def setLegend(self, action):

    if action.isChecked():
      self.legendItem.show()
    else:
      self.legendItem.hide()

  def setCoords(self, action):
  
    if action.isChecked():
      self.coordsItem.show()
    else:
      self.coordsItem.hide()
  
  def setSize(self, width, height):
    
    pos = self.legendItem.pos()
    w, h = self.size
    p1 = pos.x() * width/w
    p2 = pos.y() * height/h
    self.legendItem.setPos(p1, p2)
    
    self.setMinimumSize(width, height)
    self.size = (width, height)
    self.refresh()
  
  def refresh(self):
    
    for dataSetItem in self.dataSetItems:
      dataSetItem.refresh()
    
    x0, y0, x1, y1 = self.getPlotRegion()
    self.titleItem.set(0.5 * (x0+x1), y0-(self.zoom*25), self.title)
    self.legendItem.set(self.dataSets)
    self.xAxisItem.refresh()
    self.yAxisItem.refresh()
    
    if self.yAxis2:
      self.yAxisItem2.refresh()
  
  def addDataSet(self, dataSet, refresh=True):
    
    if not dataSet.color:
      dataSet.color = DEFAULT_COLORS[len(self.dataSets) % len(DEFAULT_COLORS)]

    dataSetItem = DataSetItem(self, dataSet)
    
    self.dataSetItems.append(dataSetItem)
    self.dataSets.append(dataSet)
    
    xMinA, xMaxA = self.xAxis.valRange
    xMinB, xMaxB = dataSet.xDataRange
    
    if dataSet.plotType == HISTOGRAM:
      binAdjust = 0.5 * float(xMaxB-xMinB) / len(dataSet.dataPoints)
      xMaxB += binAdjust
      xMinB -= binAdjust
    
    self.xAxis.valRange = (min(xMinA, xMinB), max(xMaxA, xMaxB))
    
    if dataSet.secondAxis and self.yAxis2:
      yMinA, yMaxA = self.yAxis2.valRange
      yMinB, yMaxB = dataSet.yDataRange
      self.yAxis2.valRange = (min(yMinA, yMinB), max(yMaxA, yMaxB))
   
    else:
      yMinA, yMaxA = self.yAxis.valRange
      yMinB, yMaxB = dataSet.yDataRange
      self.yAxis.valRange = (min(yMinA, yMinB), max(yMaxA, yMaxB))
    
    if refresh:
      self.refresh()
      
    self.update() # Redraw
  
  def updateData(self, dataSets, title=None, axes=None):
    # Existing axes can be updated internally before this call
    # or axes can be replaced

    if title:
      self.setTitle(title)
  
    if axes:
      self.xAxis, self.yAxis = axes[:2]
      self.xAxisItem.setAxis(self.xAxis)
      self.yAxisItem.setAxis(self.yAxis)
      
      if len(axes) == 3:
        self.yAxis2 = axes[2]
        self.yAxisItem2.setAxis(self.yAxis2)
      
    dataSetItems = []
    
    xMins = []
    xMaxs = []
    yMins = []
    yMaxs = []
    y2Mins = []
    y2Maxs = []
      
    for i, dataSet in enumerate(dataSets):
      if not dataSet.color:
        dataSet.color = DEFAULT_COLORS[len(self.dataSets) % len(DEFAULT_COLORS)]
      
      if i < len(self.dataSetItems):
        dataSetItem = self.dataSetItems[i]
        dataSetItem.setDataSet(dataSet)
      else:
        dataSetItem = DataSetItem(self, dataSet)
        self.dataSetItems.append(dataSetItem)

      a, b = dataSet.xDataRange
 
      if dataSet.plotType == HISTOGRAM:
        binAdjust = 0.5 * float(b-a) / len(dataSet.dataPoints)
        a += binAdjust
        b -= binAdjust
      
      xMins.append(a)
      xMaxs.append(b)
      
      if dataSet.secondAxis and self.yAxis2:
        a, b = dataSet.yDataRange
        y2Mins.append(a)
        y2Maxs.append(b)
 
      a, b = dataSet.yDataRange
      yMins.append(a)
      yMaxs.append(b)
       
    self.xAxis.valRange = (min(xMins or [0]), max(xMaxs or [1]))
    self.yAxis.valRange = (min(yMins or [0]), max(yMaxs or [1]))
    
    if y2Mins:
      self.yAxis2.valRange = (min(y2Mins), max(y2Maxs))
    
    while len(self.dataSetItems) > len(dataSets):
      dataSetItem = self.dataSetItems.pop()
      self.graphicsScene.removeItem(dataSetItem)
    
    self.dataSets = dataSets
    self.refresh()
    self.update() # Redraw
  
  def update(self):

    self.refresh()
  
    return QtGui.QGraphicsView.update(self)
   
  def popupContextMenu(self, pos):
    
    self.contextMenu.popup(pos)   
    
  def resetView(self):
  
    self.fitInView(self.sceneRect(), Qt.KeepAspectRatio)
  
  def _chooseColor(self, dataSet):
  
    dialog = ColorDialog(self)
    dialog.setColor(dataSet.color)
    
    color = dialog.getHexColor()
    
    if color:
      dataSet.color = color[:7]
      dataSet.gradient = None
      self.update()   
  
  def _setupColormenu(self, menu):
    
    menu.clear()
    for dataSet in self.dataSets:
      menu.addItem(dataSet.name, callback=self._chooseColor,
                   object=dataSet, icon=Icon(color=dataSet.color))

  def _setSymbol(self, obj):
    
    dataSet, symbol = obj
    dataSet.symbol = symbol
    
    self.update()
                   
  def _setupSymbolMenu(self, menu):
    
    menu.clear()
    for dataSet in self.dataSets:
      if dataSet.plotType == HISTOGRAM:
        continue
    
      subMenu = Menu(menu, dataSet.name, icon=Icon(color=dataSet.color))
    
      for opt in (CIRCLE, SQUARE, TRIANGLE, STAR, CROSS, HEXAGON):
        subMenu.addItem(opt.title(), group=1, checked=dataSet.symbol==opt,
                        callback=self._setSymbol, object=(dataSet, opt))

  def _setPlotType(self, obj):
    
    dataSet, plotType = obj
    dataSet.plotType = plotType
    
    self.update()
                   
  def _setupPlotTypeMenu(self, menu):
    
    menu.clear()
    for dataSet in self.dataSets:
      subMenu = Menu(menu, dataSet.name, icon=Icon(color=dataSet.color))
    
      for opt in (LINE, SCATTER, HISTOGRAM):
        subMenu.addItem(opt.title(), group=1, checked=dataSet.plotType==opt,
                        callback=self._setPlotType, object=(dataSet, opt))
  
  def _setSymbolSize(self, obj):
  
    dataSet, value = obj
    
    dataSet.symbolSize = value
    self.update()
  
  def _setupSymbolSizeMenu(self, menu):
    
    menu.clear()
    for dataSet in self.dataSets:
      if dataSet.plotType == HISTOGRAM:
        continue
      
      subMenu = Menu(menu, dataSet.name, icon=Icon(color=dataSet.color))
    
      for opt in range(1,21):
        subMenu.addItem(str(opt), group=1, checked=dataSet.symbolSize==opt,
                        callback=self._setSymbolSize, object=(dataSet, opt))


  def configMenu(self):
    
    mouseMenu = Menu(None, 'Graph Options')
    
    mouseMenu.addItem('Show Grid?', self.setGrid, checked=True)
    mouseMenu.addItem('Show Legend?', self.setLegend, checked=self.showLegend)
    mouseMenu.addItem('Show Cursor Coords?', self.setCoords, checked=True)
    
    colorMenu = Menu(mouseMenu, 'Colours', setupFunc=self._setupColormenu)
    
    typeMenu = Menu(mouseMenu, 'Graph Type', setupFunc=self._setupPlotTypeMenu)
      
    symTypeMenu = Menu(mouseMenu, 'Symbol Shape', setupFunc=self._setupSymbolMenu)
    
    symSizeMenu = Menu(mouseMenu, 'Symbol Size', setupFunc=self._setupSymbolSizeMenu)
 
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
    dialog.printWidgetPixmap(self)


  def askGraphTitle(self):
  
    text = askString('Text Entry','Enter graph title', self.title, self)
    self.setTitle(text)
  
 
  def drawBackground(self, painter, viewRect):
    
    QtGui.QGraphicsView.drawBackground(self, painter, viewRect) 
          
    pr0, pr1, pr2, pr3 = self.getPlotRegion()
       
    painter.setPen(self.mgColor)
    painter.setBrush(QtCore.Qt.NoBrush)
    painter.drawRect(pr0, pr1, pr2-pr0, pr3-pr1) # plot edge

  
  def wheelEvent(self, event):
    
    if event.delta() < 0:
      fac = 0.8333
    else:
      fac = 1.2
    
    self.scale(fac, fac)
  
    event.accept()

  def mouseDoubleClickEvent(self, event):
    
    QtGui.QGraphicsView.mouseDoubleClickEvent(self, event)
    
    if self.callback:
      pos = self.mapToScene(event.pos())
      x = pos.x()
      y = pos.y()
 
      plotRegion = self.getPlotRegion()
      x0, y0, x1, y1 = plotRegion
  
      if (x0 < x < x1) and (y0 < y < y1):
        px, py = self.getDataCoords(x, y)
        self.callback(self, px, py)
        
  def mousePressEvent(self, event):
   
    QtGui.QGraphicsView.mousePressEvent(self, event)
    
    # deal with inconsistency in Qt versions for button naming
    try:
      MiddleButton = Qt.MiddleButton
    except AttributeError:
      MiddleButton = Qt.MidButton
      
    button = event.button()

    if button == Qt.RightButton:
      self.popupContextMenu(event.globalPos())
    
    elif button == MiddleButton:
      pos = event.pos()
      h = self.horizontalScrollBar().sliderPosition()
      v = self.verticalScrollBar().sliderPosition()
      self.movePos = pos.x()+h, pos.y()+v

  def mouseMoveEvent(self, event):
    
    QtGui.QGraphicsView.mouseMoveEvent(self, event)
    
    if self.movePos:
      x0, y0 = self.movePos
      pos = event.pos()
      self.horizontalScrollBar().setSliderPosition(x0-pos.x())
      self.verticalScrollBar().setSliderPosition(y0-pos.y())

    else:  
      pos = self.mapToScene(event.pos())
      x = pos.x()
      y = pos.y()
 
      plotRegion = self.getPlotRegion()
      (x0,y0,x1,y1) = plotRegion
 
      if (x0 < x < x1) and (y0 < y < y1):
        if self.showCoords or self.motionCallback:
          px, py = self.getDataCoords(x, y)
 
          if self.showCoords:
            self.coordsItem.set((px, py))
 
          if self.motionCallback:
            self.motionCallback(px, py)

        self.xCrosshair.set(x, y0, x, y1)
        self.yCrosshair.set(x0, y, x1, y)
    
  def mouseReleaseEvent(self, event):
    
    self.movePos = None
    
    return QtGui.QGraphicsView.mouseReleaseEvent(self, event)
  
  def coordsOff(self):
    
    if self.coordsItem:
      self.coordsItem.set(None)
    
    self.xCrosshair.hide()
    self.yCrosshair.hide()
    
    for line in self.indicatorLinesH:
      line.hide()
    for line in self.indicatorLinesV:
      line.hide()
    
  def leaveEvent(self, event):
    
    QtGui.QGraphicsView.leaveEvent(self, event)
    self.coordsOff()
    
    # Popdown menu?

  def enterEvent(self, event):
  
    QtGui.QGraphicsView.enterEvent(self, event)    
    

  def drawVerticalLines(self, positions):
    """ draw vertical lines at specified data values
    """

    minX, maxX = self.xAxis.valRange
    minY, maxY = self.yAxis.valRange
    
    plotRegion = self.getPlotRegion()
    x0, y0, x1 ,y1 = plotRegion
    
    deltaXplot = x1 - x0
    
    data = []
    for position in positions:
      x, y = self.getPlotCoords(position, minY)
      if x0 <= x <= x1:
        data.append(x)

    if not data:
      self.coordsOff()
      
    else:
      if self.showCoords:
        self.coordsItem.set(positions)

      for n, x in enumerate(data):
        if n >= len(self.indicatorLinesV):
          line = CrosshairItem(self)
          self.indicatorLinesV.append(line)
        
        self.indicatorLinesV[n].set(x, y0, x, y1)

  def drawHorizontalLines(self, positions):
    """ draw vertical lines at specified data values
    """

    minX, maxX = self.xAxis.valRange
    minY, maxY = self.yAxis.valRange
    
    plotRegion = self.getPlotRegion()
    x0, y0, x1, y1 = plotRegion
    
    deltaYplot = y1 - y0
    
    data = []
    for position in positions:
      x, y = self.getPlotCoords(minX, position)
      if y0 <= y <= y1:
        data.append(y)

    if not data:
      self.coordsOff()
      
    else:
      if self.showCoords:
        self.coordsItem.set(positions)

      for n, y in enumerate(data):
        if n >= len(self.indicatorLinesH):
          line = CrosshairItem(self)
          self.indicatorLinesH.append(line)
        
        self.indicatorLinesH[n].set(x0, y, x1, y)
  
  def getPlotCoords(self, dataX, dataY):
    
    pr0, pr1, pr2, pr3 = self.getPlotRegion()
 
    xMin, xMax = self.xAxis.valRange
    yMin, yMax = self.yAxis.valRange
    
    if self.xAxis.isReverse:
      xMin, xMax = xMax, xMin
    
    if self.yAxis.isReverse:
      yMin, yMax = yMax, yMin
    
    xSize = xMax-xMin
    ySize = yMax-yMin
    
    plotWidth = pr2 - pr0
    plotHeight = pr3 - pr1

    ppvX = plotWidth/float(xSize or 1.0)
    ppvY = plotHeight/float(ySize or 1.0)
  
    x0 = pr0 + (dataX - xMin) * ppvX
    y0 = pr3 - (dataY - yMin) * ppvY
    
    return x0, y0
  
  def getDataCoords(self, x, y):

    xMin, xMax = self.xAxis.valRange
    yMin, yMax = self.yAxis.valRange
    plotRegion = self.getPlotRegion()

    x -= plotRegion[0]
    y -= plotRegion[1]
  
    deltaXplot = plotRegion[2] - plotRegion[0]
    deltaYplot = plotRegion[3] - plotRegion[1] 
    deltaXdata = xMax - xMin
    deltaYdata = yMax - yMin
  
    vppX = deltaXdata/float(deltaXplot)
    vppY = deltaYdata/float(deltaYplot)
  
    dataX = xMin + x * vppX
    dataY = yMax - y * vppY
    
    return dataX, dataY

  def getPlotRegion(self):
  
    x0, y0 = 0, 0    
    w, h = self.size
    w *= self.zoom/2.0
    h *= self.zoom/2.0
    xpad = 25
    ypad = 25
  
    return (x0+xpad-w, y0+ypad-h, x0+xpad+w, y0+h)     
  
  def setZoom(self, factor=1.0):
  
    if factor:
      self.scale(factor, factor)

  def setTitle(self, title):
    
    if title:
      self.title = str(title)
    else:
      self.title = ''
     
    self.titleItem.setText( self.title )
    
    self.update()

  def setSymbols(self, symbol):
  
    for dataSet in self.dataSets:
      dataSet.symbol = symbol
  
    self.update()

  def setSymbolSize(self, size):
  
    for dataSet in self.dataSets:
      dataSet.symbolSize = size

    self.update()

  def setDataNames(self, names):
    
    for i, name in enumerate(names):
      self.dataSets[i].name = name
    
    self.update()

  def setPlotType(self, plotType):

    for dataSet in self.dataSets:
      dataSet.plotType = plotType

    self.update()
    



if __name__ == '__main__':

  from math import exp

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup

  app = Application()
  
  popup = BasePopup(title='Scrolled Graph Example')
  popup.setSize(600, 400)
  
  def click(*args):
    print('Clicked', args)
  
  gradient = QtGui.QLinearGradient()
  c1 = QtGui.QColor(255, 0, 0)
  c2 = QtGui.QColor(0, 0, 255)
  gradient.setColorAt(0.0, c1)
  gradient.setColorAt(1.0, c2)  
    
  numbers = [float(x/10.0) for x in range(1,21)] 
  squared = GraphDataSet([(x,x*x) for x in numbers], u'x\xb2',
                         '#FFD0B0', plotType='histogram', gradient=gradient)
  exponent = GraphDataSet([(x,exp(x)) for x in numbers],
                          u'e\u02e3', '#008000', symbol='star', symbolSize=10)
   
  xAxis = GraphAxis('X Value', labels=None, ticks=True)
  yAxis = GraphAxis('Y Value A', labels=None, ticks=True, valRange=(0,16))
  yAxis2 = GraphAxis('Y Value B', labels=None, ticks=True, grid=False)
  
  axes = (xAxis, yAxis, yAxis2)
  dataSets = [squared, exponent]
  
  graph = Graph(popup, axes, dataSets, title='Example Graph',
                callback=click, grid=(0,0))  

  points = [(x,1.0/x,0.5*1.0/x) for x in numbers] #+ [(10, 0, 0),(12, 0, 0),]
  graph.addDataSet( GraphDataSet(points, u'x\u207b\xb9', '#0000C0',
                    symbolSize=7, symbol='triangle', ignoreZero=False, secondAxis=True,
                    errorColors=(QtGui.QColor(255, 200, 200),
                                 QtGui.QColor(200, 200, 255))) )

  app.start()



