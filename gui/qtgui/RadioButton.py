from PySide import QtCore, QtGui
from gui.qtgui.Base import Base

class RadioButton(QtGui.QRadioButton, Base):

  def __init__(self, parent, text='', *args, **kw):
      
    QtGui.QRadioButton.__init__(self, text, parent)
    Base.__init__(self, parent, *args, **kw)

