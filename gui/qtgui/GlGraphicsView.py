from PySide import QtGui, QtCore
from OpenGL import GL, GLU, GLUT

from gui.qtgui.Base import Base

# This is to some extent inspired by Rivo Laks's WidgetProxy ( http://websvn.kde.org/*checkout*/trunk/playground/libs/kgllib/extras/kgllib/widgetproxy.cpp )

class GlScene(QtGui.QGraphicsScene):
  
  def __init__(self, parent, glWidget):
    
    QtGui.QGraphicsScene.__init__(self, 0, 0, glWidget.width(), glWidget.height(), parent)
    self.widget = glWidget
    
  def drawBackground(self, painter, rect):
    
    self.widget.paintEvent(None, painterArg=painter)
    #self.widget.update()
    
    QtGui.QGraphicsScene.drawBackground(self, painter, rect)
        
  def resize(self):
    
    self.setSceneRect(0, 0, self.widget.width(),  self.widget.height())
    return self.sceneRect()

class GlGraphicsView(QtGui.QGraphicsView, Base):
        
  def __init__(self, parent, mainApp, glWidget, **kw):
    
    self.eventForwarding = False
    QtGui.QGraphicsView.__init__(self, parent)
    Base.__init__(self, parent,  **kw)
    self.mainApp = mainApp
    self.widget = glWidget
    self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
    self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
    self.setOptimizationFlags(QtGui.QGraphicsView.DontSavePainterState)
    self.setViewportUpdateMode(QtGui.QGraphicsView.SmartViewportUpdate)
    self.setRenderHints(QtGui.QPainter.Antialiasing)
    self.setStyleSheet("border: 0px solid")
    self.setViewport(glWidget)
    self.setScene(GlScene(self, glWidget))
    self.setViewportMargins(4, 4, 4, 4) # Otherwise there is a black margin from the glWidget that is not properly repainted.

#    self.setSizePolicy(glWidget.sizePolicy())
    
#    self.setMouseTracking(glWidget.hasMouseTracking())
    
  def resize(self):
    
    self.scene().resize()
    self.widget.resizeGL(self.widget.width(), self.widget.height())
#    self.updateSceneRect(rect)
    
  def addWidget(self, widget):
    
    return self.scene().addWidget(widget)

  def addItem(self, item):
    
    self.scene().addItem(item)
    item.setFlag(QtGui.QGraphicsItem.ItemIsMovable)

  def resizeEvent(self, event):
    
    self.resize()
    
#    QtGui.QGraphicsView.resizeEvent(self, event)
  
#  def event(self, event):
#    
#    if self.eventForwarding:
#      return False
#      
#    return self.forwardEvent(event)
    
  def forwardEvent(self, event):
    
    self.eventForwarding = True
    result = QtCore.QCoreApplication.instance().notify(self.widget, event)
    self.eventForwarding = False
    
    return result

  def viewportEvent(self, event):
    
    if self.eventForwarding:
      return False
    
    if event.type() == QtCore.QEvent.Wheel or event.type() == QtCore.QEvent.MouseButtonPress or event.type() == QtCore.QEvent.MouseButtonRelease:
      if self.handleMouseEvent(event, event.pos()):
        return True
      # Mouse release events need to be forwarded to the background. Otherwise selections in the spectra will not work properly if the mouse is released over e.g. a molecule.
      elif event.type() == QtCore.QEvent.MouseButtonRelease:
        self.forwardEvent(event)
      else:
        # Check if there are any items in the viewport that should take care of the event. Otherwise forward it to the background.
        state = QtGui.QGraphicsView.viewportEvent(self, event)
        if not event.isAccepted():
          state = self.forwardEvent(event)
        return state
    elif event.type() == QtCore.QEvent.DragEnter or event.type() == QtCore.QEvent.DragMove or event.type() == QtCore.QEvent.Drop:
      if self.handleMouseEvent(event, event.pos()):
        return True      
    else:
      self.forwardEvent(event)
      
    return QtGui.QGraphicsView.viewportEvent(self, event)

  def handleMouseEvent(self, event, pos):
    
    scene = self.scene()
    
    if scene.mouseGrabberItem():
      return False
      
    if not scene.itemAt(pos.x(), pos.y()):
      scene.setFocusItem(None)
      return self.forwardEvent(event)
      
    return False
    
