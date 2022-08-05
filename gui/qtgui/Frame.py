import sys

from PySide import QtGui

from gui.qtgui.Base import Base
from gui.qtgui.Layout import FlowLayout

class Frame(QtGui.QFrame, Base):

  def __init__(self, parent=None, **kw):

    QtGui.QFrame.__init__(self, parent)
    Base.__init__(self, parent, **kw)


class FlowFrame(QtGui.QFrame, Base):

  def __init__(self, parent=None, **kw):

    QtGui.QFrame.__init__(self, parent)
    
    layout = FlowLayout(self)
    self.setLayout(layout)
    
    Base.__init__(self, parent, **kw)


if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup

  app = Application()
  popup = BasePopup(title='Test Frame')
  popup.resize(400, 400)
  frame = Frame(parent=popup)
  app.start()

