from PySide import QtCore, QtGui
from gui.qtgui.Base import Icon

class Menu(QtGui.QMenu):

  def __init__(self, parent=None, text='', key=None, icon=None,
               callback=None, tearoff=False, persistant=False,
               setupFunc=None):

    QtGui.QMenu.__init__(self, text, parent=parent)

    if not key:
      key = text
    
    if parent:
      if hasattr(parent, 'addMenu'):
        parent.addMenu(self)
        
      if isinstance(parent, Menu):
        parent.itemDict[key] = self
  
    self.key = key
    self.text = text
    self.setupFunc = setupFunc
    self.itemDict = {}
    self._toolTipAction = None
    self._toolTipTimerId = None
    self._toolTip = None
    self._groupDict = {}
    
    self.setTearOffEnabled(tearoff)
    
    if icon:
      self.setIcon(icon)
      self.menuAction().setIconVisibleInMenu(True)
    else:
      self.menuAction().setIconVisibleInMenu(False)
      
    self.callback = callback
    
    self.connect(self, QtCore.SIGNAL("hovered(QAction *)"), self._toolTipHover)
    self.connect(self, QtCore.SIGNAL('triggered(QAction *)'), self._callback)
    
    if persistant:
      self.connect(self, QtCore.SIGNAL('aboutToHide()'), self._stayUp)
    
    if self.setupFunc:
      self.connect(self, QtCore.SIGNAL('aboutToShow()'), self._setupFunc)
  
  def _setupFunc(self):
  
    self.setupFunc(self)
  
  def _stayUp(self):
     
    if self.underMouse():
      self.show()

  def _callback(self, action):
    
    if self.callback:
      obj = action.data()
      
      if obj:
        self.callback(obj)
      else:
        self.callback(action)
  
  def getActions(self):
  
    return self.actions()
  
  def _fetchAction(self, keyOrIndex):
    
    if keyOrIndex in self.itemDict:
      return self.itemDict[keyOrIndex]
  
    elif type(keyOrIndex) is type(1):
      return self.actions()[keyOrIndex]

  def setItemChecked(self, keyOrIndex, bool):
  
    action = self._fetchAction(keyOrIndex)
    
    if action and action.isCheckable():
      action.setChecked(bool)
  
  def enableItem(self, keyOrIndex):
  
    self.setItemEnabled(keyOrIndex, True)
  
  def disableItem(self, keyOrIndex):
  
    self.setItemEnabled(keyOrIndex, False)
  
  def setItemEnabled(self, keyOrIndex, bool):
    
    action = self._fetchAction(keyOrIndex)
    
    if action:
      action.setEnabled(bool)
  
  def deleteItem(self, keyOrIndex):
    
    action = self._fetchAction(keyOrIndex)
    
    if action:
      self.removeAction(action)
  
  def insertItem(self, index, text, callback=None, object=None,
                 icon=None, checked=None,  key=None, shortcut=None,
                 tipText=None, group=None):
  
    self.addItem(text, callback, object, icon, key, shortcut,
                 tipText, index, group)
                 
  # addMenu(), addSeparator() inbuilt
  
  def actionEvent(self, event):
  
    QtGui.QMenu.actionEvent(self, event)
    
    # keep itemDict up-to-date whenever
    # anything changes a menu's actions
  
    aware = set(self.itemDict.values())
    for action in self.actions():
      if action not in aware:
        self.itemDict[action.text()] = action
      
  
  def addItem(self, text, callback=None, object=None, icon=None,
              key=None, checked=None, shortcut=None, tipText=None,
              index=None, group=None, widget=None):
    
    if widget:
      action = QtGui.QWidgetAction(self.parent())
      
      if text:
        frame = QtGui.QWidget(self)
        layout = QtGui.QHBoxLayout(frame)
        layout.setStretch(0,1)
        layout.setStretch(1,0)
        layout.setSpacing(4)
        layout.setContentsMargins(4,4,4,4)
        frame.setLayout(layout)
        widget.setParent(frame)
        
        label = QtGui.QLabel(text, frame)
        layout.addWidget(label)
        layout.addWidget(widget)
        
        action.setDefaultWidget(frame)
        
      else:  
        action.setDefaultWidget(widget)
    
    elif icon:
      action = QtGui.QAction(Icon(icon), text, self.parent())
      action.setIconVisibleInMenu(True)
    
    else:
      action = QtGui.QAction(text, self.parent())
    
    if checked is not None:
      action.setCheckable(True)
      action.setChecked(checked)     
    
    if shortcut:
      action.setShortcut(shortcut)
      
    if tipText:
      action.setStatusTip(tipText)
      action.setToolTip(tipText)
      
    if callback:
      if object is not None:
        func = lambda x=object:callback(x)
      elif checked is not None:  
        func = lambda x=action:callback(x)
      else:
        func = callback
      
      action.connect(action, QtCore.SIGNAL("triggered()"), func)
    
    if object:
      action.setData(object)
    
    if not key:
      key = text

    self.itemDict[key] = action
    
    if index is None:
      self.addAction(action)
    else:
      beforeAction = self.actions()[index]
      self.insertAction(beforeAction, action)
    
    if group is not None:
      if group in self._groupDict:
        actionGroup = self._groupDict[group]
      else:
        actionGroup = QtGui.QActionGroup(self)
        self._groupDict[group] = actionGroup
    
      actionGroup.addAction(action)
    
    return action
    
  def timerEvent(self, event):
    
    if event.timerId() is self._toolTipTimerId:
      tipText = self._toolTipAction.toolTip()
      self._toolTip = QtGui.QToolTip
      self._toolTip.showText(QtGui.QCursor.pos(), tipText)
    
      self.killTimer(event.timerId())
      self._toolTipTimerId = None
  
    QtGui.QMenu.timerEvent(self, event)

  def leaveEvent(self, event):
  
    if self._toolTipTimerId is not None:
      self.killTimer(self._toolTipTimerId)
      self._toolTipAction = None
      self._toolTipTimerId = None
 
    QtGui.QMenu.leaveEvent(self, event)

  def _toolTipHover(self, action):
    
    if action is not self._toolTipAction:
      self._toolTipAction = action 
    
      if self._toolTipTimerId is not None:
        self.killTimer(self._toolTipTimerId)
      
      if self._toolTip:
        self._toolTip.hideText()
        self._toolTip = None
        
      self._toolTipTimerId = self.startTimer(2000)

if __name__ == '__main__':

  from Application import Application
  from MainWindow import MainWindow

  def newProject():
    print('newProject')

  def callbackM(obj):
    print("Menu Selected Obj:", obj)

  def callback(obj):
    print("Item Selected Obj:", obj)
  
  app = Application()
  mainWindow = MainWindow(title='Test Menu')
  menuBar = mainWindow.menuBar()
  menu = Menu(menuBar, mainWindow.tr('&Project'), 'Project', callback=callbackM)
  menu.addItem('New', key=mainWindow.tr('&New'), shortcut='Ctrl+N', tipText='Create New Project', callback=newProject)
  menu.addItem('Something Else', tipText='Another option')
  menu.addItem('Selected', tipText='A checkable option, starts on', checked=True)
  menu.addItem('Unselected', tipText='A checkable option, starts off', checked=False)
  menu.addItem('Unselected', tipText='A checkable option, starts off', checked=False)
  menu.addItem('Object', callback, object=app, tipText='Item with object')
  app.start()

