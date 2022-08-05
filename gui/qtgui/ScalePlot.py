from math import log
from PySide import QtGui, QtCore
Qt = QtCore.Qt
QPointF = QtCore.QPointF
QRectF = QtCore.QRectF
QColor = QtGui.QColor

import sys, os
from os.path import dirname

guiDir = os.path.abspath(dirname(dirname(dirname(__file__))))


if guiDir not in sys.path:
  sys.path.append(guiDir)

from gui.qtgui.Base import Base, Icon
from gui.qtgui.Colors import ColorDialog
from gui.qtgui.FileSelect import selectSaveFile, FileType
from gui.qtgui.graphicsItems.LabelItems import LabelItem, MovableLabelItem
from gui.qtgui.InputDialog import askString
from gui.qtgui.Menu import Menu
from gui.qtgui.Print import PrintDialog

GREY_PEN = QtGui.QPen(QColor(0, 0, 0), 0.5)
TRANSPARENT = QColor(0, 0, 0, 0)
HIGHLIGHT_PEN = QtGui.QPen(QColor(255, 255, 255, 64), 2, Qt.DotLine)
HIGHLIGHT_BRUSH = QColor(255, 255, 255, 32)
NULL_RECT = QRectF()
NULL_POINT = QPointF()

class ScaleItem(QtGui.QGraphicsItem):
  
  def __init__(self, parent):

    QtGui.QGraphicsItem.__init__(self, scene=parent.scene())
    
    self.parent = parent
    self.setZValue(2)
    self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
    self.setAcceptedMouseButtons(Qt.LeftButton)
    self.fontMetric = QtGui.QFontMetricsF(self.parent.font())
    
    self.sync()
  
  def sync(self):
    
    bbox = self.fontMetric.tightBoundingRect
    getColor = self.parent.getColor
    scalePoints = self.parent.scalePoints
    minValue = min(scalePoints)
    maxValue = max(scalePoints)
    length = self.parent.length
    width = self.parent.width
    isVertical = self.parent.isVertical
    clipped = self.parent.clipped

    deltaVal0 = (maxValue-minValue)/20.0
    if not deltaVal0:
      deltaVal0 = 1.0
    order = int(round(log(deltaVal0, 10)))
    
    if -5 < order < 6:
      textFormat  = '%%.%df' % max(1,-order)
    else:
      textFormat = '%.2e'
        
    textHeight = bbox('x').height()
    pad = textHeight
    ppv = float(length)/(maxValue-minValue)
    
    rects = []
    texts = []
    colors = []
    j = len(scalePoints)-1
    
    for i, scalePoint in enumerate(scalePoints):
      colors.append(getColor(i, scalePoint))
    
    if isVertical:
      for i, scalePoint in enumerate(scalePoints):
        
        if i == j:
          scalePoint2 = scalePoints[i-1]
          delta = scalePoint-scalePoint2
        
        else:
          scalePoint2 = scalePoints[i+1]
          delta = scalePoint2-scalePoint
        
        y = (scalePoint-minValue)*ppv
        text = textFormat % scalePoint
        
        if clipped and i in (0,j):
          if scalePoint < 0:
            text = u'\u2264' + text
          else:
            text = u'\u2265' + text
                               
        rects.append(QRectF(0, y, width, delta*ppv))
        texts.append((width + pad , y+textHeight/2.0+ width/2, text) )
      
      textWidths = [bbox(t[2]).width() for t in texts]
      dy = length
      dx = width + pad + max(textWidths)
      outline = QRectF(0, 0, width, length)
    
    else:
      for i, scalePoint in enumerate(scalePoints):
        if i == j:
          scalePoint2 = scalePoints[i-1]
          delta = scalePoint-scalePoint2
        
        else:
          scalePoint2 = scalePoints[i+1]
          delta = scalePoint2-scalePoint
        
        x = (scalePoint-minValue)*ppv
        dx = delta*ppv
        text = textFormat % scalePoint
        
        if clipped and i in (0,j):
          if scalePoint < 0:
            text = u'\u2264' + text
          else:
            text = u'\u2265' + text

        rects.append(QRectF(x, 0, dx, width))
        tw = bbox(text).width()
        texts.append( (x-tw/2+dx/2, width + pad + pad, text) )
        
      #text = textFormat % scalePoint2
      #texts.append( (length - (bbox(text).width())/2.0, width + pad + pad, text) )
      dy = width + pad + textHeight
      dx = length
      outline = QRectF(0, 0, length, width)
    
    self.boundingRegion = QRectF(0, 0, dx, dy) 
    self.drawData = (outline, rects, texts, colors) 
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
    
    outline, rects, texts, colors = self.drawData
    
    setPen(GREY_PEN)
    drawRect(outline)
    
    for i, rect in enumerate(rects):
      setBrush(colors[i])
      drawRect(rect)
      
    setPen(Qt.black)
    for t in texts:
      x, y, text = t
      drawText(x, y, text)

class ScalePlot(QtGui.QGraphicsView, Base):

  def __init__(self, parent, scalePoints, getColor=None,
               length=300, width=25, title='Scale Plot',
               axisName='Scale', isVertical=False, showGrid=True, clipped=False, **kw):
      
    QtGui.QGraphicsView.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.getColor = getColor or self._defaultColor
    self.parent = parent
    self.title = title
    self.axisName = axisName
    self.showGrid = showGrid
    self.movePos = None
    self.isVertical = isVertical
    self.length = length
    self.width = width
    self.clipped = clipped
    
    scalePoints = list(scalePoints)
    scalePoints.sort()
    
    if len(scalePoints) == 2:
      minValue, maxValue = scalePoints
      scalePoints = []
      n = length/width
      d = (maxValue-minValue)/float(n)
      v = minValue
      while v < maxValue:
        scalePoints.append(v)
        v += d
    
    self.scalePoints = scalePoints

    self.setRenderHint(QtGui.QPainter.Antialiasing)
    self.setResizeAnchor(QtGui.QGraphicsView.AnchorUnderMouse)
    self.setTransformationAnchor(QtGui.QGraphicsView.NoAnchor)
    self.setViewportUpdateMode(QtGui.QGraphicsView.FullViewportUpdate)

    self.gScene = QtGui.QGraphicsScene(self)
    self.setScene(self.gScene)   
    
    self.titleItem = MovableLabelItem(self)
    self.scaleItem = ScaleItem(self)
    
    if isVertical:
      self.labelItem = MovableLabelItem(self, angle=-90)
      self.titleItem.set(0, -width/2, self.title)
      self.labelItem.set(-width/2, length/2, self.axisName)
    else:
      self.labelItem = MovableLabelItem(self)
      self.titleItem.set(0, -width, self.title)
      self.labelItem.set(length/2, -width/2, self.axisName)
    
    self.contextMenu = self.configMenu()
    self.show()
  
  def _defaultColor(self, i, region):
  
    if i % 2 == 0:
      return QColor(192,192,192,255)
    else:
      return QColor(64,64,64,255)
    
  
  def askTitle(self):
  
    text = askString('Text Entry','Enter plot title', self.title, self)
    self.setTitle(text)
  
  def setTitle(self, title=''):
    
    self.title = title
    self.titleItem.setText(title)
    
    self.update()

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
        
    mouseMenu = Menu(None, 'Plot Options')

    imageMenu = Menu(mouseMenu, 'Export Image')
    
    imageMenu.addItem('JPEG', self.exportJpeg)
    imageMenu.addItem('PNG', self.exportPng)
        
    mouseMenu.addItem('Print', self.printWidget)
 
    mouseMenu.addItem('Set Title', self.askTitle)
   
    mouseMenu.addItem('Reset View', self.resetView)
    
    
    return mouseMenu
  
  def printWidget(self):
  
    dialog = PrintDialog(self) 
    dialog.printWidgetPixmap(self)
  
  def exportJpeg(self):
  
  
    fileTypes = [FileType('JPEG', ['*.jpg','.jpeg','.jpr']),]
    filePath = selectSaveFile(self, caption='Select image file name',
                              directory=None, fileTypes=fileTypes)
                   
    if filePath:
      widget = self
      pixmap = QtGui.QPixmap.grabWidget(widget, widget.rect())
      pixmap.save(filePath, 'JPEG')
    
  def exportPng(self):
  
    fileTypes = [FileType('PNG', ['*.png']),]
    filePath = selectSaveFile(self, caption='Select image file name',
                              directory=None, fileTypes=fileTypes)
                   
    if filePath:
      widget = self
      pixmap = QtGui.QPixmap.grabWidget(widget, widget.rect())
      pixmap.save(filePath, 'PNG')

if __name__ == '__main__':
  
  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from numpy import array
  
  def rainbow(i, value):
    color = QColor() 
    color.setHsvF(value/1e7, 1.0, 1.0)
    
    return color
  
  app = Application()
  
  popup = BasePopup(title='Scale Plot Example')
  popup.setSize(800, 800)
  popup.setStyleSheet("background-color:white;")

  #ScalePlot(popup, (-10,-5,0,5,10), isVertical=False)
  
  #ScalePlot(popup, (0,1e7), rainbow, 600, 90, isVertical=False) 

  def orangeBlue(i, value):
    value /= 4.0
    if value < -1.0:
      value = -1.0
      
    if value > 1.0:
      value = 1.0  
    
    if value < 0.0:
      b = -value
      g = 0.5 * b
      r = 0.0
    else:
      r = value
      g = 0.5 * r
      b = 0.0
    
    color = QColor()
    color.setRgbF(r,g,b)
    
    return color
  
  nums =  array(range(-8,9))/2.0 

  ScalePlot(popup, nums, orangeBlue, 600, 60, isVertical=False, clipped=True) 

  app.start()
