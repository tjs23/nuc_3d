
# TBC: Set icons, tree selection

from PySide import QtCore, QtGui

from gui.qtgui.Base import Base
from gui.qtgui.List import List
from gui.qtgui.LabelFrame import LabelFrame

class PageWidget(LabelFrame):
  """
  Subclass this to add custom panels to PageView
  Re-implement setObject(obj)
  """

  def __init__(self, parent=None, **kw):

    LabelFrame.__init__(self, parent, 'Default Page', **kw)
    self.object = None
  
  def setObject(self, obj):
    """
    Overwrite in subclass
    """
    
    self.setText(str(obj))
    self.object = obj
  

class PageView(QtGui.QWidget, Base):

  def __init__(self, parent, pageWidgetClass,
               getTextFunc, objects=None, **kw):

    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    
    self.pageWidgetClass = pageWidgetClass
    self.getTextFunc = getTextFunc
    
    self.listBox = List(self, grid=(0,0), multiSelect=False,
                        callback=self._choose, stretch=(0,0))
                        
    self.pageWidget = pageWidgetClass(self, grid=(0,1), stretch=(0,2))
    self.pageWidget.show()
    
    if objects:
      self.setObjects(objects)
      

  def _choose(self, obj):
  
    self.pageWidget.setObject(obj)
  
  
  def setObjects(self, objects):
    
    self.objects = objects or []
    
    texts = [self.getTextFunc(obj) for obj in self.objects]
    
    self.listBox.setItems(texts, self.objects) #, icons, tipTexts)


if __name__ == '__main__':

  import sys

  from math import exp

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  
  from gui.qtgui.Label import Label
  
  app = Application(sys.argv)
  popup = BasePopup(title='Test PageView')
  
  class NumberPage(PageWidget):
  
    def __init__(self, parent=None, **kw):
      
      PageWidget.__init__(self, parent, **kw)
      
      self.labels = []
      
      for i in range(4):
        self.labels.append(Label(self, "Row %d" % (i+1), grid=(i, 0)))
      
      if self.layout():
        self.layout().setRowStretch(i+1, 2)
      
    def setObject(self, x):
      
      x = x[0]
         
      self.labels[0].set('x: %d'     % x)
      self.labels[1].set('x^2: %d'   % (x*x))
      self.labels[2].set('1/x: %.3f' % (1.0/x))
      self.labels[3].set('e^x: %.3f' % exp(x))
    
      self.object = x
      self.setText('Number: %d' % x)

  pageView = PageView(popup, NumberPage, str, list(range(1,8)))
  
  app.start()
