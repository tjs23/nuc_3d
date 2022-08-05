import sys

from PySide import QtGui, QtCore

from gui.qtgui.Base import Base
from gui.qtgui.Button import Button

# Replaces PartitionedSelector
# TBD work with colors, e.g. for spectra
# More prominant visual when button is down

class ButtonArray(QtGui.QWidget, Base):

  def __init__(self, parent, texts=None, objects=None, icons=None, callback=None,
               toggled=True, radio=False, selected=None, tipTexts=None, 
               maxCols=None, **kw):

    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    buttonGroup = self.buttonGroup = QtGui.QButtonGroup(self)
    buttonGroup.setExclusive(radio)
    
    self.maxCols = maxCols
    self.buttons = []
    self.objects = []
    self.radio = radio
    self.toggled = toggled
    self.radio = radio
     
    if not texts:
      texts = []
    
    n = len(texts)
 
    if not objects:
      objects = texts[:]
 
    if not icons:
      icons = [None] * n

    if not tipTexts:
      tipTexts = [None] * n
 
    assert len(objects) == n
    assert len(tipTexts) == n
    assert len(icons) == n
    
    self.setItems(texts, objects, icons, tipTexts, selected)

    buttonGroup.connect(buttonGroup, QtCore.SIGNAL('buttonClicked(int)'), self._callback)

    self.setCallback(callback)
  
  def setItems(self, texts, objects, icons=None, tipTexts=None, selected=None):
  
    selected = set(selected or [])
    n = len(texts)
    
    if not icons:
      icons = [None] * n

    if not tipTexts:
      tipTexts = [None] * n
    
    self.objects = []
    for i, obj in enumerate(objects):
      isSelected = obj in selected
      self.addItem(texts[i], obj, icons[i], tipTexts[i], isSelected)
    
    nButtons = len(self.buttons)
    nObjects = len(objects)
    
    if nButtons > nObjects:
      for i in range(nObjects, nButtons):
        self.buttons[i].hide()
    
  def addItem(self, text, object, icon, tipText, isSelected=None):
    
    m = self.maxCols
    i = len(self.objects)
    
    if i < len(self.buttons):
      button = self.buttons[i]
      button.setText(text)
      button.setIcon(icon or QtGui.QIcon())
      button.setToolTip(tipText)
      button.show()
      
    elif m:
      row = i // m
      col = i % m
      button = Button(self, text, icon=icon, tipText=tipText, grid=(row, col))
      button.setMinimumWidth(20)
      self.buttons.append(button)
      self.buttonGroup.addButton(button)
      self.buttonGroup.setId(button, i)
    
    else:
      button = Button(self, text, icon=icon, tipText=tipText, grid=None)
      button.setMinimumWidth(20)
      self.buttons.append(button)
      self.buttonGroup.addButton(button)
      self.buttonGroup.setId(button, i)
 
    button.setCheckable(self.toggled)
    
    self.objects.append(object)
    
    if isSelected is not None:
      self.buttons[i].setChecked(isSelected)
  
  def getSelected(self):
    
    indices = [i for i, b in enumerate(self.buttons) if b.isChecked()]
    
    return [self.objects[i] for i in indices]
  
  def setSelected(self, objs):
  
    allObjs = set(self.objects)
    indices = set()
    
    for obj in objs:
      if obj in allObjs:
        i = self.objects.index(obj)
      elif type(obj) is type(1):    
        i = obj
      else:
        continue  
        
      indices.add(i)
    
    if not indices:
      self.buttonGroup.setExclusive(False)
    
    for i in range(len(self.buttons)):
      bool = i in indices
      self.buttons[i].setChecked(bool)  

    self.buttonGroup.setExclusive(self.radio)
  
  def setButtonSelected(self, i, bool):
  
    self.buttons[i].setChecked(bool)
  
  def selectButton(self, i, doCallback=True):
    
    self.buttons[i].setChecked(True)
    
    if doCallback:
      self.callback(self.objects[i])
  
  def setCallback(self, callback):

    self.callback = callback

  def _callback(self, ind):

    if self.callback and ind >= 0:
      button = self.buttonGroup.buttons()[ind]
      self.callback(self.objects[ind])

if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  
  def callback(obj):
    print("Toggled", obj)
  
  app = Application()
  popup = BasePopup(title='Test Frame')
  popup.resize(400, 60)
  
  toolbar = ButtonArray(popup, ['One','Two','Three'],
                        [1,2,3], callback=callback)
  
  toolbar.selectButton(2)
  toolbar.setButtonSelected(0, True)
  
  print(toolbar.getSelected())

  vals = range(10,35)
  toolbar = ButtonArray(popup, [str(x) for x in vals],
                        vals, callback=callback, maxCols=9,
                        selected=[10,], radio=True)
    
  app.start()

