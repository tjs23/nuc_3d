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
from gui.qtgui.FileSelect import selectSaveFile, FileType

GREY_PEN = QtGui.QPen(QColor(128, 128, 128), 0.8)
BLACK_PEN = QtGui.QPen(QColor(0, 0, 0), 0.8)
TRANSPARENT = QColor(0, 0, 0, 0)
HIGHLIGHT_PEN = QtGui.QPen(QColor(255, 255, 255, 64), 2, Qt.DotLine)
HIGHLIGHT_BRUSH = QColor(255, 255, 255, 32)
NULL_RECT = QRectF()
NULL_POINT = QPointF()

    
    
from numpy import array, empty, ones, uint8

class MarkersItem(QtGui.QGraphicsItem):
  
  def __init__(self, parent, x0, y0, length, points):

    QtGui.QGraphicsItem.__init__(self, scene=parent.scene())
    
    self.parent = parent
    self.setZValue(2)
    #self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
    #self.setAcceptedMouseButtons(Qt.LeftButton)
    self.fontMetric = QtGui.QFontMetricsF(self.parent.font())
    self.length = length
    self.points = points
        
    bbox = self.fontMetric.tightBoundingRect
    th = bbox('x').height()

    
    lines = []
    texts = []
    triangles = []
    
    l = float(length)
    
    for val, text in points:
      x = x0 + (val * l)
      tw = bbox(text).width()/2.0
      xL = x-th/2.0
      xR = xL+th
      
      triangles.append( QtGui.QPolygonF( (QPointF(xL,y0-th), QPointF(x,y0), QPointF(xR,y0-th)), ) ) 
      texts.append( (QPointF(x-tw,y0-th*2.6), text) )
      lines.append( (QPointF(x,y0-th), QPointF(x,y0-2.0*th)) )
    
    self.boundingRegion = QRectF(0, -4*th, length, 4*th) 
    self.drawData = (lines, texts, triangles) 
    self.update()
    
    
  def boundingRect(self):
    
    if self.boundingRegion:
      return self.boundingRegion # .adjust(-2,-2, 2, 2)
    
    else:
      return NULL_RECT      
      
  def paint(self, painter, option, widget):
    
    drawLine = painter.drawLine
    drawText = painter.drawText
    setPen   = painter.setPen
    setBrush = painter.setBrush
    drawPolygon = painter.drawPolygon
    
    lines, texts, triangles = self.drawData
    
    setPen(BLACK_PEN)
    
    for p1, p2 in lines:
      drawLine(p1, p2)
    
    for t in texts:
      p1, text = t
      drawText(p1, text)
      
    for p in triangles:
      drawPolygon(p)

class ScaleItem(QtGui.QGraphicsItem):
  
  def __init__(self, parent, x0, y0, length, values):

    QtGui.QGraphicsItem.__init__(self, scene=parent.scene())
    
    self.parent = parent
    self.setZValue(2)
    #self.setFlag(QtGui.QGraphicsItem.ItemIsMovable)
    #self.setAcceptedMouseButtons(Qt.LeftButton)
    self.fontMetric = QtGui.QFontMetricsF(self.parent.font())
    self.length = length
    self.values = values
        
    bbox = self.fontMetric.tightBoundingRect
    textHeight = bbox('x').height()
    w = bbox('x').width()/2
    
    lines = [(QPointF(x0,y0), QPointF(x0+length,y0))]
    texts = []
    
    dx = length / float(len(values)-1)
    x = x0
    for value in values:
      line = (QPointF(x,y0), QPointF(x,y0+textHeight*0.75))
      text = str(value)
      tw = bbox(text).width()/2.0
      texts.append((QPointF(x-tw,y0+textHeight*2.5), text))
      lines.append(line)
      x += dx
    
    self.boundingRegion = QRectF(-w, 0, length+w+w, 2.5*textHeight) 
    self.drawData = (lines, texts) 
    self.update()
    
    
  def boundingRect(self):
    
    if self.boundingRegion:
      return self.boundingRegion # .adjust(-2,-2, 2, 2)
    
    else:
      return NULL_RECT      
      
  def paint(self, painter, option, widget):
    
    drawLine = painter.drawLine
    drawText = painter.drawText
    setPen   = painter.setPen
    setBrush = painter.setBrush
    
    lines, texts = self.drawData
    
    setPen(BLACK_PEN)
    
    for p1, p2 in lines:
      drawLine(p1, p2)
    
    for t in texts:
      p1, text = t
      drawText(p1, text)

class PixmapPlot(QtGui.QGraphicsView, Base):

  def __init__(self, parent, array, scale=None, markers=None, **kw):
      
    QtGui.QGraphicsView.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    if not scale:
      scale = [x*10 for x in range(11)]
    
    if not markers:
      markers = []
      
    self.parent = parent

    self.setRenderHint(QtGui.QPainter.Antialiasing)
    self.setResizeAnchor(QtGui.QGraphicsView.AnchorUnderMouse)
    self.setTransformationAnchor(QtGui.QGraphicsView.NoAnchor)
    self.setViewportUpdateMode(QtGui.QGraphicsView.FullViewportUpdate)

    self.gScene = QtGui.QGraphicsScene(self)
    self.setScene(self.gScene)   

    n, m, d = array.shape
    qImage = QtGui.QImage(array.data, m, n, QtGui.QImage.Format_RGB32)

    pixmap = QtGui.QPixmap.fromImage(qImage) #.scaled(m,m)
    self.pixmapItem = QtGui.QGraphicsPixmapItem(pixmap)
    self.gScene.addItem(self.pixmapItem)
    
    bbox = self.pixmapItem.boundingRect() 
    x = self.pixmapItem.x()
    y = self.pixmapItem.y()
    length = bbox.width()
    
    self.markersItrem = MarkersItem(self,x,y,length,markers)
   
    y += bbox.height()
    self.scaleItem = ScaleItem(self,x,y,length,scale)
 
    self.contextMenu = self.configMenu()

  def exportJpeg(self, filePath=None):  
    
    if not filePath:
      fileTypes = [FileType('JPEG', ['*.jpg','.jpeg','.jpr']),]
      filePath = selectSaveFile(self, caption='Select image file name',
                                directory=None, fileTypes=fileTypes)
 
                   
    if filePath:
      widget = self
      pixmap = QtGui.QPixmap.grabWidget(widget, widget.rect())
      pixmap.save(filePath, 'JPEG')
    
  def exportPng(self, filePath=None):
  
    if not filePath:
      fileTypes = [FileType('PNG', ['*.png']),]
      filePath = selectSaveFile(self, caption='Select image file name',
                                directory=None, fileTypes=fileTypes)
                   
    if filePath:
      widget = self
      pixmap = QtGui.QPixmap.grabWidget(widget, widget.rect())
      pixmap.save(filePath, 'PNG')
    
  def wheelEvent(self, event):
    
    if event.delta() < 0:
      fac = 0.8333
    else:
      fac = 1.2
    
    self.scale(fac, fac)
  
    event.accept()
    
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
   
  def popupContextMenu(self, pos):
    
    self.contextMenu.popup(pos)   

  def configMenu(self):
    
    mouseMenu = Menu(None, 'Pixmap Options')
  
    imageMenu = Menu(mouseMenu, 'Export Image')
    
    imageMenu.addItem('JPEG', self.exportJpeg)
    imageMenu.addItem('PNG', self.exportPng)
    
    return mouseMenu

if __name__ == '__main__':

  from math import exp
  
  
  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup

  app = Application()
  
  popup = BasePopup(title='Pixmap Plot Example')
  popup.setSize(800, 600)
  
  data = empty((100, 1000, 4), float, order='C')
  
  for i in range(1000):
    r = int(255 * i/1000.0)
    g = 128
    b = 255 - r
    a = 255
    data[:,i] = (r,g,b,a)
  
  data -= data.min()
  data /= data.max()
  
  pixmapArray = array(255 * data, uint8, order='C')
  
  markers = [(0.25, 'Quarter'), (0.5, 'Half'), (0.75, 'Three-quarters')]
  
  plot = PixmapPlot(popup, pixmapArray, range(11), markers)

  app.start()

