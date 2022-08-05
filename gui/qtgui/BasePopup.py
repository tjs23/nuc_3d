from PySide import QtGui, QtCore


    
class _PopupCore(object):

  def __init__(self, parent=None, title='', location=None, hide=False,
               modal=False, transient=False, quitFunc=None,
               tipText=None):
    
    self.location = location
    self.isTransient = transient
    self.modal = modal
    self.quitFunc = quitFunc
    
    if modal: # Set before visible
      modality = QtCore.Qt.ApplicationModal
      self.setWindowModality(modality)
    
    if tipText:
      self.setToolTip(tipText)


    if parent and not location:
      x = parent.x() + 50
      y = parent.y() + 50
      
      rect = self.rect()
      w = rect.width()
      h = rect.height()
      
      location = (x, y, w, h)
      
    if location:
      self.show()
      self.setGeometry(*location)

    if hasattr(parent, 'top'):
      self.top = parent.top
    else:
      self.top = self

    """
    if (parent and transient):
      self.transient(parent)
    """

    self.setWindowTitle(title)

    self.body(self)

    if hide:
      self.hide()
    else:
      self.show()
      self.raise_()
      
  def closeEvent(self, event):
  
    if self.quitFunc:
      self.quitFunc()
    
    QtGui.QWidget.closeEvent(self, event)
  
      
  def waitCursor(self):
    
    print("qt.BasePopup.waitCursor() not implemented")
    
    return
    
    #cursor = self.cursor()
    #self.setCursor(QtCore.Qt.WaitCursor)
    #self.setCursor(cursor)
 
  def open(self):

    self.showNormal()
    self.lift()
    self.activateWindow()

  def close(self, *event):
    
    self.updateLocation()
    self.hide()

  def lift(self):

    self.raise_()

  def lower(self, below=None):

    if below:
      self.stackUnder(below)
    else:
      self.lower()  
  
  def setSize(self, w, h):
  
    self.setGeometry(self.x(), self.y(), w, h)
  
  def setGeometry(self, x, y, w=None, h=None):
    
    if w is None:
      w = self.rect().width()
    
    if h is None:
      h = self.rect().width()  
     
    screenRect = QtGui.QApplication.desktop()
       
    sWidth = screenRect.width()
    sHeight = screenRect.height()
    
    w = min(sWidth, w)
    h = min(sHeight, h)
    
    if (x+w) > sWidth:
      if w == 1:
        x = sWidth // 2
      else:
        x = sWidth - w
    elif x < 0:
      x = 0
    
    if (y+h) > sHeight:
      if h == 1:
        y = sHeight // 2
      else:
        y = sHeight - h  
    elif y < 0:
      y = 0
          
    QtGui.QWidget.setGeometry(self, x, y, w, h)


  def updateLocation(self):

    self.location = self.getGeometry()

  def getGeometry(self):

    rect = self.geometry()
 
    return rect.x(), rect.y(), rect.width(), rect.height()

  def iconify(self):

    if self.isHidden():
      self.showNormal()
    
    self.showMinimized()
    
  def deiconify(self):

    self.showNormal()
    
  def withdraw(self):

    self.hide()

  def destroy(self):

    QtGui.QWidget.destroy(self)

  def setTitle(self, title=''):

    self.setWindowTitle(title)


  def body(self, master):
    
    pass # this method should be overridden by subclass

  def apply(self):
    
    return True # this method can be overridden by subclass

class BasePopup(QtGui.QWidget, _PopupCore):

  def __init__(self, parent=None, title='', location=None, hide=False,
               modal=False, transient=False, quitFunc=None,
               tipText=None):
    

    QtGui.QWidget.__init__(self, parent=None)
    _PopupCore.__init__(self, parent, title, location, hide,
                        modal, transient, quitFunc, tipText)

class DockPopup(QtGui.QDockWidget, _PopupCore):

  def __init__(self, parent=None, title='', location=None, hide=False,
               modal=False, transient=False, quitFunc=None,
               tipText=None):
    

    QtGui.QDockWidget.__init__(self, parent=None)
    _PopupCore.__init__(self, parent, title, location, hide,
                        modal, transient, quitFunc, tipText)
    
    self.frame = QtGui.QWidget()
    self.setWidget(self.frame)
    
if __name__ == '__main__':

  import sys
  from gui.qtgui.Button import Button
  from gui.qtgui.Label import Label
  from gui.qtgui.Application import Application
  
  app = Application()

  class TestPopup(BasePopup):

    def __init__(self, parent):

      BasePopup.__init__(self, parent, 'Test Popup', modal=False)

      self.result = None

    def body(self, master):

      self.setGeometry(600, 400, 50, 50)
      
      label = Label(master, text='label 1', grid=(0,0))
      label = Label(master, text='label 2', grid=(1,0))
      button = Button(master, text='ok', callback=self.apply, grid=(2,0))
      button = Button(master, text='cancel', callback=self.close, grid=(3,0))

    def apply(self):
      self.result = 77
      return True

  popup = None
  window = QtGui.QWidget()

  def new():

    global popup, window
    popup = TestPopup(window)
    popup.show()
    print(popup.result)

  def lift():

    if (popup):
      popup.lift()

  def close():

    if (popup):
      popup.close()

  def open_():

    if (popup):
      popup.open()

  def iconify():

    if (popup):
      popup.iconify()

  def deiconify():

    if (popup):
      popup.deiconify()


  button = Button(window, text='new popup', callback=new)

  button = Button(window, text='lift popup', callback=lift)

  button = Button(window, text='open popup', callback=open_)

  button = Button(window, text='close popup', callback=close)

  button = Button(window, text='iconify popup', callback=iconify)

  button = Button(window, text='deiconify popup', callback=deiconify)

  button = Button(window, text='quit', callback=window.destroy)

  window.show()
  
  app.start()


