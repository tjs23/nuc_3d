from math import ceil

from PySide import QtCore, QtGui

from gui.qtgui.Base import Base


class WellPlateWidget(QtGui.QWidget, Base):
  
  wellSelected = QtCore.Signal(int)

  def __init__(self, parent = None, objects = None, colorKeyWord = None, baseColor = None, **kw):
    
    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.objects = objects
    self.colorKeyWord = colorKeyWord
    if(objects):
      self.numObjects = len(objects)
    else:
      self.numObjects = 96
    self.baseColor = baseColor
    
    if not (6 <= self.numObjects <= 384):
      raise Exception('Number of objects (%d) not in range 6 <= n <= 384.' % self.numObjects)
      
    if self.numObjects == 6:
      nCols = 3
    elif self.numObjects <= 24:
      nCols = 6
    elif self.numObjects <=96:
      nCols = 12
    else:
      nCols = 24
    nRows = int(ceil(self.numObjects / float(nCols)))
    
    if objects:
      maxColorValue = None
      for object in objects:
        value = getattr(object, colorKeyWord)
        if maxColorValue == None or value > maxColorValue:
          maxColorValue = value
          
    for n in xrange(self.numObjects):
      row = n/nCols
      col = n%nCols
      if objects:
        object = objects[n]
        value = getattr(object, colorKeyWord)
        if maxColorValue:
          colorIntensity = value*(235/maxColorValue)
        else:
          colorIntensity = 120

        text = 'Well %d' % (n+1)
        if colorKeyWord:
          text += ', %s: %d' % (colorKeyWord, value)
        Well(self, n, baseColor = self.baseColor, colorIntensity=colorIntensity, size=20, pos=(col*20, row*20), text=text)
      else:
        Well(self, n, baseColor = self.baseColor, size=20, pos=(col*20, row*20), text='Well %d' % (n+1))
      
    self.wellSelected.connect(wellSelectedReaction)

class Well(QtGui.QWidget):
  
  selected = QtCore.Signal(int)
  
  def __init__(self, parent, index, baseColor = 'blue', colorIntensity = 225, pos = None, size = None, text = None):
    
    QtGui.QWidget.__init__(self, parent)
    
    self.index = index
    if not size:
      size = 20
    if not pos:
      pos = (0, 0)
    self.setGeometry(pos[0], pos[1], size, size)
    
    self.setColor(baseColor or 'blue', colorIntensity)
    
    if not text:
      text = 'Index', index 
        
    self.setToolTip(text)
    self.selected.connect(parent.wellSelected)
    
  def setColor(self, baseColor, colorIntensity):
    
    highlight = max(200, colorIntensity)
    
    if baseColor == 'red':
      self.centralColor = QtGui.QColor(highlight, 175, 175)
      self.edgeColor = QtGui.QColor(colorIntensity, 0, 0)
    elif baseColor == 'green':
      self.centralColor = QtGui.QColor(175, highlight, 175)
      self.edgeColor = QtGui.QColor(0, colorIntensity, 0)
    elif baseColor == 'yellow':
      self.centralColor = QtGui.QColor(highlight, highlight, 175)
      self.edgeColor = QtGui.QColor(colorIntensity, colorIntensity, 0)
    elif baseColor in ('grey', 'gray'):
      self.centralColor = QtGui.QColor(highlight, highlight, highlight)
      self.edgeColor = QtGui.QColor(colorIntensity, colorIntensity, colorIntensity)
    else:
      self.centralColor = QtGui.QColor(175, 175, highlight)
      self.edgeColor = QtGui.QColor(0, 0, colorIntensity)
    
  def paintEvent(self, event):
    
    painter = QtGui.QPainter(self)
    painter.setRenderHints(QtGui.QPainter.Antialiasing)
    
    rect = QtCore.QRectF(1, 1, 18, 18)
    
    gradient = QtGui.QRadialGradient(7, 7, 4)
    gradient.setColorAt(0, self.centralColor)
    gradient.setColorAt(1, self.edgeColor)
    
    painter.setBrush(gradient)
    painter.drawEllipse(rect)
    
  def mousePressEvent(self, event):
    
    self.selected.emit(self.index)

def wellSelectedReaction(index):
  
  from gui.qtgui.MessageDialog import showInfo
  
  print('Well number %d was selected' % (index+1))

class DemoClass:
  
  def __init__(self, name = None, score = 0):
    
    self.name = name
    self.score = score

if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  import random
  
  app = Application()
  popup = BasePopup(title='Test Frame')
  popup.setSize(500, 330)
  
  objects = []
  
  for i in xrange(96):
    objects.append(DemoClass(score=random.randrange(0, 4)))
  
  WellPlateWidget(popup, objects, "score", baseColor = 'blue')

  app.start()
