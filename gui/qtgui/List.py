from PySide import QtCore, QtGui

from gui.qtgui.Base import Base, Icon

SELECTED = 'icons/selected-yes.png'
UNSELECTED = 'icons/selected-no.png'

class List(QtGui.QListWidget, Base):

  def __init__(self, parent, texts=None, callback=None, objects=None,
               icons=None, tipTexts=None, doubleCallback=None,
               multiSelect=True, **kw):
    
    QtGui.QListWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.iconYes = Icon(SELECTED)
    self.iconNo = Icon(UNSELECTED)
    
    if multiSelect:
      self.setSelectionMode(self.MultiSelection)
    
    self.callback = callback
    self.doubleCallback = doubleCallback
    
    if texts:
      self.setItems(texts, objects, icons, tipTexts)
    
    #self.itemClicked.connect(self._callback)
    self.itemDoubleClicked.connect(self._doubleCallback)
    self.itemSelectionChanged.connect(self._callback)
  
  def _refreshIcons(self):

    for i in range(self.count()):
      item = self.item(i)
  
      if not item.data(33):
        if item.isSelected():
          item.setIcon(self.iconYes)
        else:
          item.setIcon(self.iconNo)
    
  def _callback(self):
    
    self._refreshIcons()
    if self.callback:
      selection = [it.data(32) or it.text() for it in self.selectedItems()]
      self.callback(selection)

  def _doubleCallback(self):
  
    if self.doubleCallback:
      item = self.currentItem()
      self.doubleCallback(item.data(32) or item.text())
  
  
  def _makeItem(self, text, object=None, icon=None, tipText=None):
  
    if icon:
      item = QtGui.QListWidgetItem(QtGui.QIcon(icon), text, self)
    else:
      item = QtGui.QListWidgetItem(self.iconNo, text, self)
    
    if object:
      item.setData(32, object)
    
    if icon:
      item.setData(33, True)
    else:
      item.setData(33, False)  
    
    if tipText:
      item.setToolTip(tipText)
    
    return item
    
               
  def addItem(self, text, object=None, icon=None, tipText=None, selected=False):
    
    if isinstance(text, QtGui.QListWidgetItem):
      # preserves original use
      item = text
    else:      
      item = self._makeItem(text, object, icon, tipText)
    
    item.setSelected(selected)
    
    return QtGui.QListWidget.addItem(self, item)


  def insertItem(self, index, text, object=None, icon=None, tipText=None):
  
    item = self._makeItem(text, object, icon, tipText)
    QtGui.QListWidget.insertItem(self, item)


  def setItems(self, texts, objects=None, icons=None, tipTexts=None):
    
    self.clear()
    
    n = len(texts)
    
    if not objects:
      objects = [None] * n
      
    if not icons:
      icons = [None] * n
      
    if not tipTexts:
      tipTexts = [None] * n
    
    for i, text in enumerate(texts):
      item = self._makeItem(text, objects[i], icons[i], tipTexts[i])
      QtGui.QListWidget.addItem(self, item)


  def deleteSelected(self):
    
    for item in self.selectedItems():
      self.removeItemWidget(witemdget)


  def deleteIndex(self, index):
    
    item = self.itemFromIndex(index)
    self.removeItemWidget(item)


  def deleteText(self, text):
    
    for item in self.findItems(text, QtCore.Qt.MatchExactly):
      self.removeItemWidget(item)
    
  
  def getIndex(self):
    
    return self.currentIndex()
    
    
  def getText(self):
  
    item = self.currentItem()
    
    if item:
      return item.text()
  
  
  def getObject(self):
    
    item = self.currentItem()
    
    if item and item.data(1):
      return item.data(1)
    

  def getSelectedIndices(self):

    return [self.itemToIndex(item) for item in self.selectedItems()]
    
    
  def getSelectedTexts(self):
    
    return [item.text() for item in self.selectedItems()]


  def getSelectedObjects(self):

    return [item.data(1) for item in self.selectedItems()]
    
    
  def setSelectedIndices(self, indices):
    
    indices = set(indices)
    
    for i in range(self.count()):
      item = self.item(i)
    
      if i in indices:
        item.setSelected(True)
      else:
        item.setSelected(False)
  
  
  def selectAll(self):
  
    for i in range(self.count()):
      item = self.item(i)
      item.setSelected(True)
  
  
  def selectNone(self):

    for i in range(self.count()):
      item = self.item(i)
      item.setSelected(False)


if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.Button import Button
  
  app = Application()
  popup = BasePopup(title='Test Frame')
  popup.resize(400, 150)
  
  def callback(obj):
    print("Clicked object:", obj)
  
  texts = ['Alpha','Beta','Gamma','Delta','Epsilon']
  objects = [1, 2, 3, 4, 5]
  
  List(popup, texts, callback, objects)
  
  app.start()

