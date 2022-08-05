import sys

from PySide import QtGui, QtCore

from gui.qtgui.Base import Base

Qt = QtCore.Qt

AREA_DICT = {'t':Qt.TopToolBarArea,'b':Qt.BottomToolBarArea,
             'l':Qt.LeftToolBarArea,'r':Qt.RightToolBarArea}

class ToolBar(QtGui.QToolBar, Base):

  def __init__(self, parent, name, callbacks=None,
               icons=None, texts=None, shortcuts=None, 
               objName=None, areas='tblr', iconSize=None,
               floatable=False, movable=True,
               isVertical=False, orientCallback=None, **kw):   
    
    QtGui.QToolBar.__init__(self, name, parent)
    
    if kw: # Only put in a layout in specified
      Base.__init__(self, parent, **kw)
    
    self.orientCallback = orientCallback
    
    self.actions = []
    
    if callbacks:
      self.setActions(callbacks, icons, texts, shortcuts)
    
    if objName:
      self.setObjectName(objName)
    
    if iconSize:
      if isinstance(iconSize, (tuple, list)):
        w, h = iconSize
      else:
        w = h = iconSize
        
      self.setIconSize(QtCore.QSize(w,h))
    
    self.setAreas(areas)
    
    self.setFloatable(floatable)
    self.setMovable(movable)
    
    if isVertical:
      self.setOrientation(Qt.Vertical)
    else:
      self.setOrientation(Qt.Horizontal)
    
    if orientCallback:
      self.orientationChanged.connect(self._orientCallback)


  def getActions(self):
    
    return self.actions


  def setAreas(self, areas):
    
    if isinstance(self.parent(), QtGui.QMainWindow):
      self.parent().addToolBar(self)
      
    if areas:
      areas = [x[0] for x in sorted(areas)]
      
      area = 0
      for a in areas:
        area |= AREA_DICT[a]
     
      if ('t' not in areas) and ('b' not in areas):
        self.setOrientation(Qt.Vertical)
     
    else:
      area = 0
    
    self.setAllowedAreas(area)
  
  
  def _addActions(self, callbacks, icons=None, texts=None, shortcuts=None):
  
    n = len(callbacks)
    actions = []
    
    if not icons:
      icons = [None] * n
 
    if not texts:
      texts = [None] * n

    if not shortcuts:
      shortcuts = [None] * n
 
    for i, func in enumerate(callbacks):
      icon = icons[i]
      text = texts[i]
      shortcut = shortcuts[i]
 
      if icon is None:
        icon = QtGui.QIcon()
 
        if text is None:
          text = str(i+1)
 
      elif not isinstance(icon, QtGui.QIcon):
        icon = QtGui.QIcon(icon) # Assumed to be an image file path
 
      action = QtGui.QAction(icon, text, self,
                             shortcut=shortcut,
                             triggered=func)
      
      self.addAction(action)
      actions.append(action)
    
    self.actions += actions
    return actions
    
    
  def addActions(self, callbacks, icons=None, texts=None, shortcuts=None):
  
    return self._addActions(callbacks, icons, texts, shortcuts)
  
  
  def setActions(self, callbacks, icons=None, texts=None, shortcuts=None):
    
    if self.actions:
      self.clear()
      self.actions = []
    
    return self._addActions(callbacks, icons, texts, shortcuts)
       
  
  def _orientCallback(self, orientation):
    
    isVertical = orientation == Qt.Vertical
    self.orientCallback(isVertical)

if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.Button import Button
  from gui.qtgui.Base import Align
  from PySide import QtCore

  def buttonName(text):
    print("Toggled ", text)

  app = Application()
  popup = BasePopup(title='Test Frame')
  popup.resize(400, 400)
  toolBar = ToolBar(popup, 'Example', [], vAlign=Align.top)
  
  actionData = [('Button 1', QtGui.QKeySequence("1")),
		('Button 2', QtGui.QKeySequence("2")),
		('Button 3', QtGui.QKeySequence("3")),
		('Button 4', QtGui.QKeySequence("4")),
		('Button 5', QtGui.QKeySequence("5")),
		]

  signalMapper = QtCore.QSignalMapper(app)
  
  for text, keys in actionData:
    action = QtGui.QAction(text, popup, shortcut=keys)
    signalMapper.setMapping(action, text)
    action.triggered.connect(signalMapper.map)
    toolBar.addAction(action)
    
  signalMapper.mapped.connect(buttonName)

  app.start()

