from PySide import QtCore, QtGui
from numpy import ndarray

QColor = QtGui.QColor

BLANK = QColor()

def inverseGrey(color):

  r, g, b, a = color.getRgb()
  
  m = (11*r + 16*g + 5*b)/32
  
  if (m > 192) or (m < 64):
    m = 255-m
  elif m<128:
    m += 128
  elif m<192:
    m -= 128
  
  return QColor(m, m, m)
 
def clamp(x, minval, maxval):
  return min(max(x, minval), maxval)

def floatColorToHex(r, g, b, a=None):
  if a is None:
    colors = r, g, b
  else:
    colors = r, g, b, a
  color = (clamp(int(c * 255), 0, 255) for c in colors)
  return "#{0:02x}{1:02x}{2:02x}".format(*color)

class ColorDialog(QtGui.QColorDialog):

  def __init__(self, parent=None, doAlpha=False, **kw):
    
    QtGui.QColorDialog.__init__(self, parent)
       
    self.setOption(self.ShowAlphaChannel, doAlpha)
    self.setOption(QtGui.QColorDialog.DontUseNativeDialog,  True)
    self.aborted = False
    self.rejected.connect(self.quit)
  
  def set(self, color):
  
    self.setColor(color)
    
  def setColor(self, color):
    # color can be name, #hex, (r,g,b) or (r,g,b,a)
    
    if isinstance(color, (list, tuple, ndarray)) and len(color):

      if isinstance(color[0], float):
        color = [int(255*c) for c in color]
      
      qColor = QtGui.QColor(*color)

    elif isinstance(color, QtGui.QColor):
      qColor = QtGui.QColor(color)
   
    elif color[0] == '#':
      color = color.upper()
    
      if len(color) == 9:
        r = int(color[1:3], 16)
        g = int(color[3:5], 16)
        b = int(color[5:7], 16)
        a = int(color[7:9], 16)
        color = (r, g, b, a)
        
      else:
        r = int(color[1:3], 16)
        g = int(color[3:5], 16)
        b = int(color[5:7], 16)
        color = (r, g, b)
      
      qColor = QtGui.QColor(*color)
    
    self.setCurrentColor(qColor)
  
  def quit(self):
    
    self.aborted = True
    
  def get(self):
  
    color = self.currentColor()
    
    if self.aborted:
      return None
    else:
      return color
    
  def getColor(self, initialColor=None):
    
    if initialColor is not None:
      self.setColor(initialColor)
    
    self.exec_()
    
    color = self.currentColor()
    
    if self.aborted:
      return None
    else:
      return color

  def getHexColor(self, initialColor=None):
    
    if initialColor:
      self.setColor(initialColor)
    
    self.exec_()
    
    color = self.currentColor()
    
    if self.aborted:
      return None
    else:
      return '#%02X%02X%02X%02X' % color.getRgb()

from gui.qtgui.Base import Base, Icon
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.Frame import Frame
from gui.qtgui.Button import Button

import colorsys

DEFAULT_COLORS = []

for h in range(12):
  hue = h / 12.0
 
  for val, sat in ((0.5, 1.0), (0.75, 1.0), (1.0, 1.0), (1.0, 0.6), (1.0, 0.2)):
      hexCol = floatColorToHex(*colorsys.hsv_to_rgb(hue, sat, val))
      DEFAULT_COLORS.append(hexCol)

  hexCol = floatColorToHex(*colorsys.hsv_to_rgb(1.0, 0.0, hue))
  DEFAULT_COLORS.append(hexCol)
  

class ColorPulldown(PulldownList):

  def __init__(self, parent, colors=None, callback=None,
              index=0, numRows=6, **kw):
  
    PulldownList.__init__(self, parent, texts=None, objects=None,
                          icons=None, callback=callback, index=0,
                          **kw)
    
    if not colors:
      colors = DEFAULT_COLORS[:]
    
    self.view = None
    self.dialogItem = None
    self.numRows = numRows
    self.colors = colors
    self.objects = [None] * len(colors)
    
    self.setData(colors, index)
    
    self.object = self.objects[index]
    self.disconnect(self, QtCore.SIGNAL('currentIndexChanged(int)'), self._callback)
    self.connect(self, QtCore.SIGNAL('activated(int)'), self._callback)
        
  def addItem(self, text, object=None, icon=None):
    
    i = len(self.objects)
    row = int( i % self.numRows)
    col = int( i // self.numRows)
    model = self.model()
    
    if icon:
      item = QtGui.QStandardItem(icon, text)
    else:
      item = QtGui.QStandardItem(text)
    
    model.setItem(row, col, item)
    
    if object is None:
      object = text
    
    self.texts.append(text)  
    self.objects.append(object)
    
    return item
  
  def setColor(self, color):
    
    if isinstance(color, QtGui.QColor):
      color = '#%02X%02X%02X' % color.getRgb()[:3]
    
    else:
      color = color.upper()
  
    self.object = color
    
    if color in self.colors:
      i = self.colors.index(color)
      row = i % self.numRows
      col = int(i // self.numRows)
      self.setModelColumn(col)
      self.setCurrentIndex(row)        
 
    
  def setData(self, colors, index=0):
    
    texts = [''] * len(colors)
    objects = [QtGui.QColor(c) for c in colors]
    icons = [Icon(color=c) for c in colors]
    
    nCols = int(len(objects) // self.numRows) + 1

    self.view = view = QtGui.QTableView()
    
    view.setModel(self.model())
    view.horizontalHeader().setVisible(False)
    view.verticalHeader().setVisible(False)

    self.setView(view)
    
    PulldownList.setData(self, texts, objects, index, icons)
    
    view.resizeColumnsToContents()
    view.resizeRowsToContents()
    view.setMinimumWidth(view.horizontalHeader().length())

  def currentObject(self):
    
    if self.objects:
      itemIndex = self.view.selectionModel().currentIndex()
      index = itemIndex.row() + (itemIndex.column() * self.numRows)
      if index >= 0:
        return self.objects[index]
  
  def _callback(self, index):
    
    if index < 0:
      return
    
    itemIndex = self.view.selectionModel().currentIndex()
    index = itemIndex.row() + (itemIndex.column() * self.numRows)
    
    self.setModelColumn(itemIndex.column())
    self.setCurrentIndex(itemIndex.row())        
    
    if self.objects:
      object = self.objects[index]
    else:
      object = None        
       
    if object is self.object:
      return
     
    self.object = object
    if self.callback:
      self.callback(self.object)
    
    return True
    
Qt = QtCore.Qt

class GradientWidget(QtGui.QWidget, Base):

  def __init__(self, parent, colors=None, **kw):
  
    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
  
    if not colors:
      colors = ['#FF0000','#00FF00','#0000FF']
  
    self.colors = colors
    
    self.setMinimumHeight(22)

  
  def paintEvent(self, event):
    
    painter = QtGui.QPainter(self)
    
    w = self.width()
    h = self.height()
    
    p1 = QtCore.QPointF(0.0,h/2.0)
    p2 = QtCore.QPointF(float(w), h/2.0)
    gradient = QtGui.QLinearGradient(p1, p2)
    
    n = float(len(self.colors)) - 1.0
    for i, color in enumerate(self.colors):
      gradient.setColorAt(i/n, QtGui.QColor(color))
    
    painter.setBrush(gradient)
    painter.setPen('#000000')
    painter.drawRect(0, 0, w-1, h-1)
    
    return QtGui.QWidget.paintEvent(self, event)

  
  def setColors(self, colors):
  
    self.colors = colors
    self.update()
    
  
class GradientEditor(Frame):

  def __init__(self, parent, colors=None, callback=None, **kw):
   
    Frame.__init__(self, parent, **kw)

    self.gradientWidget = GradientWidget(self, grid=(0,0))
    self.frame = Frame(self, grid=(1,0))
    
    self.addButton = Button(self, grid=(0,1), icon=Icon('gui/qtgui/icons/list-add.png'),
                            callback=self.addColor)
    self.remButton = Button(self, grid=(1,1), icon=Icon('gui/qtgui/icons/list-remove.png'),
                            callback=self.removeColor)
                            
    self.layout().setRowStretch(0, 2)
    self.layout().setColumnStretch(0, 2)
        
    size = QtCore.QSize(12,12)                       
    self.addButton.setIconSize(size)
    self.remButton.setIconSize(size)
    
    self.callback = callback
    self.colorPulldowns = []
    self.colors = []
    
    if not colors:
      colors = self.gradientWidget.colors[:]
    
    self.setColors(colors) 
  
  
  def addColor(self):
  
    colors = self.colors
    colors = colors + [colors[-1]]
    self.setColors(colors)
  
    if self.callback:
      self.callback(self.colors)
  
  def removeColor(self, index=None):
    
    colors = self.colors
    if len(colors) > 2:
      if index is None:
        index = len(colors)-1

      colors.pop(index)
      self.setColors(colors)  
    
    
      if self.callback:
        self.callback(self.colors)
    
  def setColors(self, colors):
    
    self.colors = colors
    nColors = len(colors)
    nWidgets = len(self.colorPulldowns)
    
    while nWidgets > nColors:
      widget = self.colorPulldowns.pop()
      widget.hide()
      del widget
      nWidgets -= 1
    
    while nColors > nWidgets:
      pulldown = ColorPulldown(self.frame, grid=(0,nWidgets),
                               callback=lambda c, i=nWidgets:self.selectColor(c, i))
      self.colorPulldowns.append(pulldown)
      nWidgets += 1
    
    for i, color in enumerate(colors):
      self.colorPulldowns[i].setColor(color)
      
    self.gradientWidget.setColors(self.colors)
   
      
  def selectColor(self, color, index):
    
    if isinstance(color, QtGui.QColor):
      color = '#%02X%02X%02X' % color.getRgb()[:3]
      
    self.colors[index] = color
    self.gradientWidget.setColors(self.colors)
    
    if self.callback:
      self.callback(self.colors)
    
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.Button import Button
  
  def getColor():
  
    dialog = ColorDialog()
    dialog.setColor('#008000')
    
    print('Colour chosen:', dialog.getHexColor())
    
  def selectColor(color):
  
    print('Selected', color)
  
  app = Application()
  
  window = QtGui.QWidget()
  
  button = Button(window, 'Choose color...', callback=getColor, grid=(0,0))
  
  pulldown = ColorPulldown(window, callback=selectColor, grid=(0,1))
  
  pulldown.setColor('#FF8040')
  
  widget = GradientEditor(window, grid=(1,0), gridSpan=(1,2))
  
  window.show()
  
   # opens the dialog
  
  # You don't need to start if the dialog is used alone
  app.start()
