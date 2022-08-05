import os
from PySide import QtCore, QtGui
from h5py import File, Group, Dataset, AttributeManager
from gui.qtgui.Tree import Tree

USER_ROLE = QtCore.Qt.UserRole
EDIT_ROLE = QtCore.Qt.EditRole
DISPLAY_ROLE = QtCore.Qt.DisplayRole
TOOLTIP_ROLE = QtCore.Qt.ToolTipRole
STATUS_ROLE = QtCore.Qt.StatusTipRole
BG_ROLE = QtCore.Qt.BackgroundRole
FG_ROLE = QtCore.Qt.ForegroundRole
CHECK_ROLE = QtCore.Qt.CheckStateRole
ICON_ROLE = QtCore.Qt.DecorationRole
SIZE_ROLE = QtCore.Qt.SizeHintRole
ALIGNMENT_ROLE = QtCore.Qt.TextAlignmentRole

NO_PROPS = QtCore.Qt.NoItemFlags
CHECKABLE = QtCore.Qt.ItemIsUserCheckable
ENABLED = QtCore.Qt.ItemIsEnabled
SELECTABLE = QtCore.Qt.ItemIsSelectable
EDITABLE = QtCore.Qt.ItemIsEditable

CHECKED = QtCore.Qt.Checked
UNCHECKED = QtCore.Qt.Unchecked
HORIZONTAL = QtCore.Qt.Horizontal
VERTICAL = QtCore.Qt.Vertical

QColor = QtGui.QColor
QIcon = QtGui.QIcon
QSize = QtCore.QSize
Qt = QtCore.Qt
QModelIndex = QtCore.QModelIndex

HEADERS = ('Group hierarchy', 'size', 'dType', 'chunks', 'max size', 'compression')
ATTRS = '.attrs'

from time import time

def usingTimer(function):
  
  def timer(*args, **kw):
    start = time()
    output = function(*args, **kw)
    end = time()
    print(('%16s %.6f' % (function.__name__, end-start)))
    
    return output
  
  return timer
  
# TBD 
# 
# Add attrs as if they were a group
# Add table for selected data
  
class HdfTreeModel(QtCore.QAbstractItemModel):
  
  def __init__(self, parent, hdfFile=None):
    
    QtCore.QAbstractItemModel.__init__(self, parent)
    
    self.tree = parent
    self.file = None
    self.nameDict = {}
    self.childDict = {}
    
    if hdfFile:
      self.setFile(hdfFile)
    
    
  def setFile(self, hdfFile):
    
    if hdfFile:
      if self.file:
        self.removeRows(0, len(self.file.values()))
        self.reset()
        self.nameDict = {}
        self.childDict = {}

      self.file = hdfFile
       
  #@usingTimer     
  def index(self, row, column, parent):
  
    if not self.file:
      return QModelIndex()
    
    #if not self.hasIndex(row, column, parent):
    #  return QModelIndex()
    nameDict = self.nameDict
    
    hdfName = nameDict.get(parent.internalPointer(), '/')
       
    if hdfName:
    
      if hdfName.endswith(ATTRS):
        hdfObj = self.file[hdfName[:-len(ATTRS)]].attrs

      else:
        hdfObj = self.file[hdfName]
 
      if hasattr(hdfObj, 'keys'):
        if hdfName in self.childDict:
          children = self.childDict[hdfName]
        
        elif hdfName.endswith(ATTRS):
          children = sorted(list(hdfObj.keys()))
          self.childDict[hdfName] = children
        
        else:
          if hdfObj.attrs.keys():
            children = sorted([ATTRS,] + list(hdfObj.keys()))
          else:
            children = sorted(list(hdfObj.keys()))
          self.childDict[hdfName] = children
 
        #children = sorted([x for x in hdfObj.keys()])
 
        if row < len(children):
 
          key = children[row]
          
          if key == ATTRS:
            childName = hdfName + ATTRS
          else:
            childName = '%s/%s' % (hdfName, key)
 
          if childName in nameDict:
            i = nameDict[childName]
          else:
            i = str(len(nameDict))
            nameDict[i] = childName
            nameDict[childName] = i
 
          return self.createIndex(row, column, i)
 
     
    return QModelIndex()
      
  
  #@usingTimer     
  def flags(self, index):

    if not index.isValid():
      return 0

    return Qt.ItemIsEnabled | Qt.ItemIsSelectable

      
  #@usingTimer     
  def _getRowInParent(self, hdfName):
    
    if hdfName == '/':
      return 0
    
    elif hdfName.endswith(ATTRS):
      return 0
    
    elif ATTRS in hdfName:
      gParentName = hdfName.split(ATTRS)[0]
      parentName = gParentName + ATTRS
    
      if parentName in self.childDict:
        children = self.childDict[parentName]
 
      else:
        parent = self.file[gParentName].attrs
        children = sorted(list(parent.keys()))
        self.childDict[parentName] = children
        
    else:
      parentName = '/'.join(hdfName.split('/')[:-1]) or '/'
      
      if parentName in self.childDict:
        children = self.childDict[parentName]
 
      else:
        parent = self.file[hdfName].parent
        if parent.attrs.keys():
          children = sorted([ATTRS,] + list(parent.keys()))
        
        else:
          children = sorted(list(parent.keys()))
          
        self.childDict[parentName] = children
    
    if hdfName in children:
      return children.index(hdfName)
      
    else:
      return 0
      
          
  #@usingTimer     
  def parent(self, index):
    
    if not index.isValid():
      return QModelIndex()
  
    hdfName = self.nameDict.get(index.internalPointer())
    
    if hdfName.endswith(ATTRS):
      parentName = hdfName[:-len(ATTRS)]
    
    elif ATTRS in hdfName:
      parentName, attrName = hdfName.split(ATTRS)
      parentName = parentName + ATTRS
    
    elif not hdfName or (hdfName not in self.file):
      return QModelIndex()
    
    else:
      parentName = '/'.join(hdfName.split('/')[:-1]) or '/'
    
    if parentName == '/':
      return QModelIndex()
    
    row = self._getRowInParent(parentName)
    
    if parentName in self.nameDict:
      i = self.nameDict[parentName]
      
    else:
      i = str(len(self.nameDict))
      self.nameDict[parentName] = i 
      self.nameDict[i] = parentName
   
    return self.createIndex(row, 0, i)
 
   
  #@usingTimer     
  def rowCount(self, parent):
    
    # number of children
       
    if parent.isValid():
      hdfName = self.nameDict.get(parent.internalPointer())
    
    elif self.file:
      hdfName = '/'
    
    else:
      return 0
    
    if hdfName.endswith(ATTRS):
      hdfObj = self.file[hdfName[:-len(ATTRS)]].attrs
      nRows = len(hdfObj.keys())
          
    elif hdfName and hdfName in self.file:
      hdfObj = self.file[hdfName]
      
      if hasattr(hdfObj, 'keys'):
        nRows = len(hdfObj)
        
        if hdfObj.attrs.keys():
          nRows += 1
 
      else: # Dataset
      
        if hdfObj.attrs.keys():
          nRows = 1
        else:
          nRows = 0
        
    else:
      nRows = 0
    
    return nRows   
    
      
  #@usingTimer     
  def columnCount(self, parent):
    
    return 4
    
    """
    
    if parent.row() < 0:
      return 1 
    
    hdfName = self.nameDict.get(parent.internalPointer())
    
    if hdfName and hdfName in self.file:
      hdfObj = self.file[hdfName]
      
      if isinstance(hdfObj, Group):
        return 2 # Name, id
 
      elif isinstance(hdfObj, Dataset):
        return 6 # Name, dtype, shape, maxShape, chunks, compression,
 
      else:
        return 1

    return 1"""

  
  #@usingTimer     
  def headerData(self, i, orientation, role):
    
    if orientation == HORIZONTAL:
      if role == DISPLAY_ROLE:
        return HEADERS[i]

    #  elif role == TOOLTIP_ROLE:
    #   return

    #  elif role == SIZE_ROLE:
    #    return QSize(size, 4 + bbox.height() * 2)


  #@usingTimer     
  def data(self, index, role):    
    
    if not index.isValid():
      return
    
    if role == DISPLAY_ROLE:
 
      hdfName = self.nameDict.get(index.internalPointer())
      
      if hdfName.endswith(ATTRS):
        col = index.column()      
        if col == 0:
          return ATTRS
        
        elif col == 1:
          hdfObj = self.file[hdfName[:-len(ATTRS)]].attrs
          return len(hdfObj.keys())
        
        return
        
      elif ATTRS in hdfName:
        col = index.column()      
        hdfName, attrName = hdfName.split(ATTRS)
        attrName = attrName[1:]
        
        if col == 0:
          return attrName
          s
        elif col == 1:
          hdfObj = self.file[hdfName].attrs[attrName]
          return ', '.join([str(x) for x in hdfObj.shape]) or '1'
          
        elif col == 2:
          hdfObj = self.file[hdfName].attrs[attrName]
          return str(hdfObj.dtype)
        
        return
      
      elif not hdfName or (hdfName not in self.file):
        return
      
      col = index.column()      
      if col == 0:
        return hdfName.split('/')[-1]
      
      hdfObj = self.file[hdfName]
      if isinstance(hdfObj, Group):
          
        if col == 1:
          return len(hdfObj)
          
      elif isinstance(hdfObj, Dataset):
        # Name, dtype, shape, maxShape, chunks, compression, 
        
        if col == 1:
          return ', '.join([str(x) for x in hdfObj.shape])
        
        elif col == 2:
          return str(hdfObj.dtype)
         
        elif col == 3:
          return ', '.join([str(x) for x in hdfObj.chunks or []])
       
        elif col == 4:
          return ', '.join([str(x) for x in hdfObj.maxshape])
        
        elif col == 5:
          return str(hdfObj.compression)

      elif isinstance(hdfObj, File):
        
        pass

      
    #elif (role == ICON_ROLE) and self.getIcon:
    #  if index.column() == 0:
    #    hdfName = index.internalPointer()
    #    icon = self.getIcon(node.data)
    #
    #    if icon:
    #      return QIcon(icon)
      
    #elif role == USER_ROLE:
    #  hdfName = index.internalPointer()
    #  return node.data
    
    elif role == TOOLTIP_ROLE:
 
      hdfName = self.nameDict.get(index.internalPointer())
      if not hdfName or (hdfName not in self.file):
        return
      
      if hdfName.endswith(ATTRS):
        return 'Attribute Group'

      elif ATTRS in hdfName:
        return 'Attribute'
        
      else:
        hdfObj = self.file[hdfName]
 
        if isinstance(hdfObj, File):
          return 'File'

        elif isinstance(hdfObj, Group):
          return 'Group'
 
        elif isinstance(hdfObj, Dataset):
          return 'Dataset'


       
    elif role == STATUS_ROLE:
 
      hdfName = self.nameDict.get(index.internalPointer())
      if not hdfName or (hdfName not in self.file):
        return
      
      if hdfName.endswith(ATTRS):
        return 'Attribute Group'

      elif ATTRS in hdfName:
        return 'Attribute'
              
      else:
        hdfObj = self.file[hdfName]
 
        if isinstance(hdfObj, File):
          return 'File'
 
        elif isinstance(hdfObj, Group):
          return 'Group'
 
        elif isinstance(hdfObj, Dataset):
          return 'Dataset'
        

    #elif role == EDIT_ROLE:
    
    #elif role == FG_ROLE:
    #  hdfName = index.internalPointer()
    #  color = self.getColor(node.data)
    #  return inverseGrey(color)

    #elif role == BG_ROLE:
    #  hdfName = index.internalPointer()
    #  color = self.getColor(node.data)
    #  return color

    #elif role == CHECK_ROLE:
    #  hdfName = index.internalPointer()
    #  value = self.getValue(node.data, index.column())
    #  if isinstance(value, bool):
    #    if value:
    #      return CHECKED
    #    else:
    #      return UNCHECKED
    #      
    #  else:
    #    return None

    #elif role == ALIGNMENT_ROLE:


class HdfTreePanel(Tree):
  
  def __init__(self, parent, hdfFile=None, callback=None,
               doubleCallback=None, **kw):
  
    Tree.__init__(self, parent, None, callback, doubleCallback, **kw)
    
    model = HdfTreeModel(self, None)
    
    self.setModel(model)

    if hdfFile:
      self.setFile(hdfFile)

  def currentChanged(self, index, prev):
    
    
    if self.callback:
      nodeName = self.model().nameDict.get(index.internalPointer())
      
      if nodeName.endswith(ATTRS):
        nodeName = nodeName[:-len(ATTRS)]
        node = self.model().file[nodeName]
        node = node.attrs

      elif ATTRS in nodeName:
        parentName, attrName = nodeName.split(ATTRS)
        node = self.model().file[parentName]
        node = node.attrs[attrName[1:]]
        
      else:
        node = self.model().file[nodeName]
      
      if node is not None:
        self.callback(node)

  def mouseDoubleClickEvent(self, event):
  
    if self.doubleCallback:
      index = self.currentIndex()
      nodeName = self.model().nameDict.get(index.internalPointer())

      if nodeName.endswith(ATTRS):
        nodeName = nodeName[:-len(ATTRS)]
        node = self.model().file[nodeName]
        node = node.attrs
        
      elif ATTRS in nodeName:
        parentName, attrName = nodeName.split(ATTRS)
        node = self.model().file[parentName]
        node = node.attrs[attrName[1:]]
        
      else:
        node = self.model().file[nodeName]
      
      if node is not None:
        mods = event.modifiers()
        haveCtrl = mods & QtCore.Qt.CTRL
        haveShift = mods & QtCore.Qt.SHIFT
 
        if haveShift or haveCtrl:
          self.doubleCallback(node, True)
 
        else:
          self.doubleCallback(node, False)

  def selectObject(self, obj):
  
    model = self.model()
   
    indices = model.match(self.rootIndex(), DISPLAY_ROLE, model.getValue(obj, 0), flags=0)
   
    toOpen = set()
    for index in indices:
      while index != QModelIndex():
        toOpen.add(index)
        index = index.parent

    for index in toOpen:
      self.expand(index)   

  def setFile(self, filePath):
 
    if os.path.exists(filePath):
      hdfFile = File(filePath, 'r')
 
      self.model().setFile(hdfFile)
      self.resizeColumnToContents(0)
 
    

if __name__ == '__main__':

  import sys
  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.FileSelect import selectFile, FileType
  from gui.qtgui.Button import Button
  
  def callback(obj):
    print(obj)
  
  dirname =  os.path.dirname 
  argv = sys.argv
  app = Application(argv[0])
  
  popup = BasePopup(title='Test HDF Tree')
  
  tree = HdfTreePanel(popup, callback=None, grid=(0,0))
  
  def selectHdfFile():
    fileTypes = [FileType('All files', ['*.*']),]
    msg = 'Select HDF file'
    filePath = selectFile(popup, msg, dirname(dirname(dirname(__file__))), fileTypes)
    
    if filePath:
      tree.setFile(filePath)
    
  but = Button(popup, 'Select HDF file', callback=selectHdfFile, grid=(1,0))
  
  if len(argv) > 1:
    tree.setFile(argv[1])
  
  app.start()
  


