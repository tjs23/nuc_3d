
from PySide import QtGui

from gui.qtgui.Base import Base

class LabelFrame(QtGui.QGroupBox, Base):

  def __init__(self, parent=None, text='', **kw):

    QtGui.QGroupBox.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.setText(text)
    
  def setText(self, text):

    self.setTitle(' ' + text + ' ')
  
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.Button import Button
  
  app = Application()
  popup = BasePopup(title='Test Frame')
  popup.resize(400, 120)
  
  frame = LabelFrame(popup, 'Frame Title 1')
  
  Button(frame, 'Button A', grid=(0,0))
  Button(frame, 'Button B', grid=(0,1))
  
  frame = LabelFrame(popup, 'Frame Title 2')
  
  Button(frame, 'Button C', grid=(0,0))
  Button(frame, 'Button D', grid=(0,1))
  
  app.start()
