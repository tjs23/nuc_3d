from PySide import QtGui, QtCore
Qt = QtCore.Qt
QPointF = QtCore.QPointF
QRectF = QtCore.QRectF
QColor = QtGui.QColor

NULL_RECT = QRectF()
NULL_POINT = QPointF()
HIGHLIGHT = QColor(200, 0, 0)

ItemIsMovable = QtGui.QGraphicsItem.ItemIsMovable
ItemPositionChange = QtGui.QGraphicsItem.ItemPositionChange

class LabelItem(QtGui.QGraphicsItem):
        
  def __init__(self, parent, angle=None, pen=Qt.black):

    QtGui.QGraphicsItem.__init__(self)
    
    self.setZValue(2)
    self.parent = parent
    self.drawPoint = NULL_POINT
    self.text = None
    self.region = NULL_RECT
    self.pen = pen
    
    if angle:
      self.setRotation(angle)
          
    self.hover = False
    self.setAcceptHoverEvents(True)
    
    self.parent.scene().addItem(self)
  
  def highlight(self, bool):
  
    if bool:
      self.pen = HIGHLIGHT
    else:
      self.pen = Qt.black  
    
    self.update()
  
  def hoverEnterEvent(self, event):
  
    self.hover = True
    self.update()

  def hoverLeaveEvent(self, event):
  
    self.hover = False
    self.update()
  
  def setText(self, text):
    
    fontMetric = QtGui.QFontMetricsF(self.parent.font())
    bbox = fontMetric.boundingRect(text)
    
    w = bbox.width()
    h = bbox.height()
    
    l = - w/2.0
    b = h/2.0
    t = -b

    self.drawPoint = QPointF(l, b)
    self.text = text
    self.region = QRectF(l, t, w, h).normalized()

    self.update()
    
  def set(self, x, y, text):
    
    self.setPos(x, y)
    self.setText(text)
      
  def boundingRect(self):
    
    if self.region:
      return self.region # .adjust(-2,-2, 2, 2)
    
    else:
      return NULL_RECT
      
  def paint(self, painter, option, widget):
  
    painter.setPen(self.pen)
    painter.drawText(self.drawPoint, self.text)
       

class MovableLabelItem(LabelItem):
        
  def __init__(self, parent, angle=None, pen=Qt.black, **kw):
    
    LabelItem.__init__(self, parent, angle, pen, **kw)

    self.setFlag(ItemIsMovable)
    self.setAcceptedMouseButtons(Qt.LeftButton)
  
  def itemChange(self, change, value):
    
    return QtGui.QGraphicsItem.itemChange(self, change, value)
      
  def paint(self, painter, option, widget):
  
    if self.hover:
      painter.setPen(HIGHLIGHT)
    else:
      painter.setPen(self.pen)
      
    painter.drawText(self.drawPoint, self.text)
