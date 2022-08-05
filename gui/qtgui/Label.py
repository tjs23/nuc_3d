from PySide import QtGui, QtCore
Qt = QtCore.Qt

from gui.qtgui.Base import Base


class Label(QtGui.QLabel, Base):

  def __init__(self, parent, text='', hTextAlign=None, vTextAlign=None,
              textColor=None, **kw):

    QtGui.QLabel.__init__(self, text, parent)
    Base.__init__(self, parent, **kw)
    
    if hTextAlign or vTextAlign:
      alignment = self._getAlignment(hTextAlign, vTextAlign)
      self.setAlignment(alignment)
    
    if textColor:
      self.setStyleSheet('QLabel {color: %s;}' % textColor)
    
  def get(self):

    return self.text()

  def set(self, text=''):

    self.setText(text)


if __name__ == '__main__':

  from gui.qtgui.Button import Button
  from gui.qtgui.Application import Application

  msg = 'Hello world'
  count = 0

  def func():

    global count

    count += 1
    label.set(msg + ' ' + str(count))
    print(label.get())

  import sys
  app = Application()
  
  window = QtGui.QWidget()
 
  label = Label(window, text='Hello world', textColor='red', grid=(0,0))
  button = Button(window, text='Click me', command=func, grid=(0,1))

  window.show()

  app.start()
