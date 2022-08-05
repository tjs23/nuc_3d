from PySide import QtCore, QtGui
from gui.qtgui.Base import Base

CHECKED = QtCore.Qt.Checked
UNCHECKED = QtCore.Qt.Unchecked

class CheckButton(QtGui.QCheckBox, Base):

  def __init__(self, parent, text='', selected=False, callback=None, **kw):
      
    QtGui.QCheckBox.__init__(self, text, parent)
    Base.__init__(self, parent, **kw)

    self.callback = callback
   
    self.setSelected(selected)

    if callback:
      self.connect(self, QtCore.SIGNAL('stateChanged(int)'), self._callback)
      
  def _callback(self, state):
    
    if self.callback:
      self.callback(bool(state))

  def isSelected(self):

    return bool(self.checkState())

  def get(self):

    return self.isSelected()

  def set(self, selected, doCallback=True):

    if not doCallback:
      callback = self.callback
      self.callback = None
    
    self.setSelected(selected)

    if not doCallback:
      self.callback = callback

  def setSelected(self, selected):

    if selected:
      self.setCheckState(CHECKED)
    else:
      self.setCheckState(UNCHECKED)

  def toggle(self):

    self.setSelected(not self.isSelected())


if __name__ == '__main__':

  from gui.qtgui.Button import Button
  from gui.qtgui.Application import Application

  app = Application()
  
  def get_me():
    print('get_me:', c.isSelected())

  def toggle_me():
    c.toggle()

  def my_callback(selected):
    print('my_callback:', selected)

  root = QtGui.QWidget()

  c = CheckButton(root, 'Check button', callback=my_callback,
                  tipText='What this does')
  b = Button(root, text='Query state', callback=get_me)
  b = Button(root, text='Toggle state', callback=toggle_me)
  
  root.show()
  app.start()
