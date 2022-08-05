import sys

from PySide import QtGui

from gui.qtgui.Base import Base

class ScrolledFrame(QtGui.QScrollArea, Base):

  def __init__(self, parent=None, minSize=(300,600), **kw):

    QtGui.QScrollArea.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    frame = QtGui.QWidget(self)
    frame.setMinimumSize(*minSize)
    frame.setSizePolicy(QtGui.QSizePolicy.Expanding,
                        QtGui.QSizePolicy.Expanding)

    self.setWidget(frame)
    self.frame = frame
    
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.ButtonList import UtilityButtonList
  from gui.qtgui.Label import Label

  app = Application()
  popup = BasePopup(title='Test Frame')
  popup.resize(200, 400)
  frame = ScrolledFrame(parent=popup).frame

  for i in range(30):
    t = 'A text label to fill some space; %s' % i 
    label = Label(frame, text=t)

  buttons = UtilityButtonList(frame)
  
  app.start()

