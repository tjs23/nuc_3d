from PySide import QtCore, QtGui

from gui.qtgui.Base import Base, Icon

CHECKED = QtCore.Qt.Checked
UNCHECKED = QtCore.Qt.Unchecked

class Button(QtGui.QPushButton, Base):

  def __init__(self, parent, text='', callback=None, icon=None,
               toggle=None, command=None, iconSize=22, **kw):
    
    
    QtGui.QPushButton.__init__(self, text, parent)
    Base.__init__(self, parent, **kw)

    if icon: # filename or pixmap
      self.setIcon(Icon(icon))
      self.setIconSize(QtCore.QSize(iconSize,iconSize))
    
    if command:
      print("Use of qtgui.Button.command is deprecated")
      callback = command
    
    if toggle is not None:
      self.setCheckable(True)
      self.setSelected(toggle)
    
    self.textBox = QtGui.QFontMetricsF(self.font()).boundingRect
    self.callback = None
    self.setCallback(callback)
  
  def sizeHint(self):
    
    text = self.text()
    box = self.textBox(text)
    w = box.width() + 8
    
    if self.icon():
      iSize = self.iconSize()
    
      w += iSize.width()
      h = max(iSize.height(), box.height()) + 8
    
    else:
      h = box.height() + 8
      
    return QtCore.QSize(max(w, 24), max(h, 24))

  def setSelected(self, selected):
    
    if self.isCheckable(): 
      if selected:
        self.setChecked(CHECKED)
      else:
        self.setChecked(UNCHECKED)

  def setCallback(self, callback):
  
    if self.callback:
      self.disconnect(self, QtCore.SIGNAL('clicked()'), self.callback)
    
    if callback:
      self.connect(self, QtCore.SIGNAL('clicked()'), callback)
      # self.clicked.connect doesn't work with lambda, yet...
    
    self.callback = callback

  def disable(self):

    self.setDisabled(True)

  def enable(self):

    self.setEnabled(True)

  def setText(self, text):

    QtGui.QPushButton.setText(self, text)

  def setState(self, state):

    self.setEnabled(state)

  def configure(self, **options):

    print("qtgui.Button has no configure()")

  def config(self, **options):
  
    print("qtgui.Button has no config()")

class ButtonMenu(Button):

  def __init__(self, parent, menu, text='Select...',
               callback=None, **kw):
    
    Button.__init__(self, parent, text, **kw)
    
    menu.callback = callback
    
    self.setMenu(menu)   
    self.connect(menu, QtCore.SIGNAL('triggered(QAction *)'), self._setText)
    
  def _setText(self, action):
    
    text = action.text()
    self.setText(text)


if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.Menu import Menu
  from gui.qtgui.Base import Align

  app = Application()

  window = QtGui.QWidget()
  
  def click():
    print("Clicked")

  def clickObj(obj):
    print("Selected", obj)
  
  b1 = Button(window, text='Click Me', callback=click,
             tipText='Click for action',
             grid=(0, 0))

  b2 = Button(window, text='I am inactive', callback=click,
             tipText='Cannot click',
             grid=(0, 1))
  
  b2.disable()

  b3 = Button(window, text='I am green', callback=click,
             tipText='Mmm, green', bgColor='#80FF80',
             grid=(0, 2))

  b4 = Button(window, icon='icons/system-help.png', callback=click,
             tipText='A toggled icon button', toggle=True, 
             grid=(0, 3))

  menu = Menu()
  menu.addItem('One', object=1)
  menu.addItem('Five', object=5)
  menu.addItem('Eleven', object=11)
  
  b5 = ButtonMenu(window, menu, callback=clickObj,
                  tipText='A menu button', 
                  grid=(1, 1), vAlign=Align.top, stretch=(1,1))
  
  window.show()
  
  app.start()

