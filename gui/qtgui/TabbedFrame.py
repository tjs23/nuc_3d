from PySide import QtCore, QtGui

from gui.qtgui.Base import Base, Icon

class TabbedFrame(QtGui.QTabWidget, Base):

  def __init__(self, parent=None, texts=None, icons=None,
               callback=None, widgets=None, tipTexts=None, **kw):

    QtGui.QTabWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
 
    if not widgets:
      widgets = []
    
    if not texts:
      texts = []
      
    if not icons:
      icons = []
    
    n = max(len(icons), len(texts))    
    
    if not tipTexts:
      tipTexts = [None] * n
    
    while len(widgets) < n:
      widgets.append(None)
    
    while len(texts) < n:
      texts.append(None)
      
    while len(icons) < n:
      icons.append(None)      

    self.widgets = []
    self.frames = self.widgets # Deprecate
    
    for i, widget in enumerate(widgets):
      self.appendTab(texts[i], icons[i], widget, tipTexts[i])

    self.callback = callback
    
    if callback:
      self.connect(self, QtCore.SIGNAL('currentChanged(QWidget *)'), self._callback)

  def clear(self):
  
    self.widgets = []
    QtGui.QTabWidget.clear(self)
  
  def appendTab(self, text=None, icon=None, widget=None, tipText=None):

    if not widget:
      widget = QtGui.QWidget(self)
      layout = QtGui.QGridLayout(widget)
      layout.setSpacing(2)
      layout.setContentsMargins(2,2,2,2)
      widget.setLayout( layout )
      widget.setParent(self)
    
    if icon:
      self.addTab(widget, Icon(icon), text)
    else:
      self.addTab(widget, text)
    
    if tipText:
      i = self.count()
      self.setTabToolTip(i, tipText)
    
    self.widgets.append(widget)
    widget.parent = self
    widget.name = text
    
    return widget

  def _callback(self, widget):
    
    if self.callback:
      if widget in self.widgets:
        i = self.widgets.index(widget)
        self.callback(i)
  
  #setCurrentIndex inbuilt
  
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.Button import Button
  from gui.qtgui.Frame import Frame

  def callback(widget):
    print('callback', widget)

  def callback1():
    print('callback1')

  def callback2():
    print('callback2')

  app = Application()
  popup = BasePopup(title='Test TabbedFrame')
  
  texts = ['Option A','Option B','Option C','Option D']
  tipTexts = ['Tip A','Tip B','Tip C','Tip D']
  
  tabbedFrame = TabbedFrame(popup, texts, callback=callback, tipTexts=tipTexts)
  frame1, frame2, frame3, frame4 = tabbedFrame.widgets
  
  button = Button(frame1, text='Hit me 1', callback=callback1)
  button = Button(frame2, text='Hit me 2', callback=callback2)

  app.start()

