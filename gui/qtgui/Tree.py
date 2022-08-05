import os
from os import path

from PySide import QtCore, QtGui

from gui.qtgui.Base import Base
from gui.qtgui.Button import Button
from gui.qtgui.CheckButton import CheckButton
from gui.qtgui.FileSelect import FileType
from gui.qtgui.Frame import Frame
from gui.qtgui.Label import Label
from gui.qtgui.Menu import Menu
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.MessageDialog import showWarning

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
QModelIndex = QtCore.QModelIndex
HEAD_ADJUST = QSize(50, 0)

ICON_DIR = path.join(path.dirname(__file__),'icons') 

class Tree(QtGui.QTreeView, Base):

  def __init__(self, parent, model=None, callback=None,
               doubleCallback=None, **kw):
  
    QtGui.QTreeView.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    if model:
      self.setModel(model)
    
    self.callback = callback
    self.doubleCallback = doubleCallback
    self.setIndentation(20)
    self.setSortingEnabled(True)
    

  def currentChanged(self, index, prev):

    if self.callback:
      self.callback(index)

    return QtGui.QTreeView.currentChanged(self, index, prev)
    
  def setSelected(self, indices):
    
    selectionModel = self.selectionModel()
    selectionModel.clear()
    
    for index in indices:
      selectionModel.select(index, selectionModel.Select | selectionModel.Rows)  
        
  def mouseDoubleClickEvent(self, event):
    
    if self.doubleCallback:
      mods = event.modifiers()
      haveCtrl = mods & QtCore.Qt.CTRL
      haveShift = mods & QtCore.Qt.SHIFT
      
      if haveShift or haveCtrl:
        self.doubleCallback(self.currentIndex(), True)
        
      else:
        self.doubleCallback(self.currentIndex(), False)
    
    return QtGui.QTreeView.mouseDoubleClickEvent(self, event)

   
class TreeNode(object):

  def __init__(self, data, row, nCols, parent=None):
    
    self.data = data
    self.row = row
    self.nCols = nCols
    self.parent = parent
    self.children = []
        
    if parent:
      parent.children.append(self)
  
  def delete(self):
  
    for child in self.children:
      child.delete()
  
    if self.parent:
      self.parent.children.remove(self)


def namedListsToNodes(rootData, nCols):
  
  row = 0
  
  def processList(data, i, parentNode):
    
    if isinstance(data, (tuple, list)):
      node = TreeNode(data, i, nCols, parentNode)
    
      childData = data[1:]
      for j, datum in enumerate(childData):
        childNode = processList(datum, j, node)
    
    else:
      node = TreeNode(data, i, nCols, parentNode)
    
    return node
  
  rootNode = processList(rootData, row, None)
  
  return rootNode  


class ObjectTreeModel(QtCore.QAbstractItemModel):
  # Hierarchy stored as a simple list of lists
  # First item in the list is always the name
  
  # TBC : Add editing, boolean toggles, bgColor
  
  def __init__(self, parent, headings, getValue, 
               getIcon=None, getTipText=None,
               rootNode=None, tipTexts=None):
    
    QtCore.QAbstractItemModel.__init__(self, parent)
    
    
    if not tipTexts:
      tipTexts = [None] * len(headings)
    
    self.tree = parent
    self.headings = headings
    self.rootNode = None
    self.tipTexts = tipTexts # For headings/default
    self.getValue = getValue # Takes an object and col, returns the text
    self.getIcon = getIcon # Takes an object
    self.getTipText = getTipText
    
    if not rootNode:
      rootNode = TreeNode([], 0, 0, None)
    
    self.setDataNode(rootNode)
  
  def setDataNode(self, rootNode):
        
    if rootNode:
      if self.rootNode:
        self.removeRows(0, len(self.rootNode.children))
        self.reset()
    
      if isinstance(rootNode, (tuple, list)):
        self.rootNode = namedListsToNodes(rootNode, len(self.headings))
      else:
        self.rootNode = rootNode
       
  def index(self, row, column, parent):
    
    if not self.hasIndex(row, column, parent):
      return QModelIndex()
    
    if not parent.isValid():
      node = self.rootNode
    else:
      node = parent.internalPointer()
    
    if row < len(node.children):
      childNode = node.children[row]
      return self.createIndex(row, column, childNode)
    
    else:
      return QModelIndex()
  
  
  def parent(self, index):
  
    if not index.isValid():
      return QModelIndex()
  
    node = index.internalPointer()
    
    parentNode = node.parent
 
    if parentNode:
      return self.createIndex(parentNode.row, 0, parentNode)
    
    else:
      return QModelIndex()
 
   
  def rowCount(self, parent):
    
    if parent.column() > 0:
      return 0
    
    elif parent.isValid():
      parentNode = parent.internalPointer()
    
    else:
      parentNode = self.rootNode
    
    return len(parentNode.children)

  
  def columnCount(self, parent):
    
    if parent.isValid():
      parentNode = parent.internalPointer()
      return parentNode.nCols

    else:
      return self.rootNode.nCols

  
  def headerData(self, i, orientation, role):
    
    if orientation == HORIZONTAL:
      if role == DISPLAY_ROLE:
        return self.headings[i]

      elif role == TOOLTIP_ROLE:
        return self.tipTexts[i]

      elif role == SIZE_ROLE:
        texts = self.headings[i].split('\n')
        texts.sort(key=len)
        bbox = self.tree.bbox(texts[-1])
        size = max(20, bbox.width() + 32 )
        return QSize(size, 4 + bbox.height() * 2)


  def data(self, index, role):
    
    if not index.isValid():
      return
    
    if role == DISPLAY_ROLE:
      col = index.column()
      node = index.internalPointer()
      value = self.getValue(node.data, col)
      return value
    
    elif (role == ICON_ROLE) and self.getIcon:
      if index.column() == 0:
        node = index.internalPointer()
        icon = self.getIcon(node.data)

        if icon:
          return QIcon(icon)
      
    elif role == USER_ROLE:
      node = index.internalPointer()
      return node.data
    
    elif role == TOOLTIP_ROLE:
      node = index.internalPointer()
      
      if self.getTipText:
        return self.getTipText(node.data)
      
      elif not node.children:
        return self.tipTexts[index.column()]
      
    elif role == STATUS_ROLE:
      node = index.internalPointer()
      if not node.children:
        return self.tipTexts[index.column()]
        
    #elif role == EDIT_ROLE:
    
    #elif role == FG_ROLE:
    #  node = index.internalPointer()
    #  color = self.getColor(node.data)
    #  return inverseGrey(color)

    #elif role == BG_ROLE:
    #  node = index.internalPointer()
    #  color = self.getColor(node.data)
    #  return color

    elif role == CHECK_ROLE:
      node = index.internalPointer()
      value = self.getValue(node.data, index.column())
      if isinstance(value, bool):
        if value:
          return CHECKED
        else:
          return UNCHECKED
          
      else:
        return None

    #elif role == ALIGNMENT_ROLE:
    

class ObjectTree(Tree):
  
  def __init__(self, parent, headings, getValue, getIcon=None,
               getTipText=None, getDragId=None, dataNode=None, tipTexts=None,
               callback=None, doubleCallback=None, mimeType='application/x-ccpn',
               **kw):
  
    Tree.__init__(self, parent, None, callback, doubleCallback, **kw)
    
    model = ObjectTreeModel(self, headings, getValue, getIcon,
                            getTipText, dataNode, tipTexts)
    
    
    self.headings = headings
    self.getDragId = getDragId
    self.mimeType = mimeType
    self.fontMetric = QtGui.QFontMetricsF(self.font())
    self.bbox = self.fontMetric.boundingRect
    self.setModel(model)

  def currentChanged(self, index, prev):
    
    node = index.internalPointer()
    
    if self.callback and node and node.data:
      self.callback(node.data)

  def mouseDoubleClickEvent(self, event):
    
    if self.doubleCallback:
      index = self.currentIndex()
      node = index.internalPointer()
      
      if node and node.data:
        mods = event.modifiers()
        haveCtrl = mods & QtCore.Qt.CTRL
        haveShift = mods & QtCore.Qt.SHIFT
 
        if haveShift or haveCtrl:
          self.doubleCallback(node.data, True)
 
        else:
          self.doubleCallback(node.data, False)

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

  def setData(self, dataList):
    
    """
    old contents were lingering
    oldModel = self.model()
    
    getValue = oldModel.getValue
    getIcon = oldModel.getIcon
    getTipText = oldModel.getTipText
    tipTexts = oldModel.tipTexts
    
    model = ObjectTreeModel(self, self.headings, getValue, getIcon,
                            getTipText, dataList, tipTexts)
    
    self.setModel(model)"""
    
    self.model().setDataNode(dataList)
   
    for i in range(len(self.headings)):
      self.resizeColumnToContents(i)
    
  def dragEnterEvent(self, event):
    
    event.ignore()

  def dragMoveEvent(self, event):
    
    event.ignore()
  
  def mousePressEvent(self, event):

    Tree.mousePressEvent(self, event)
    
    if not self.getDragId:
      return
      
    index = self.indexAt(event.pos())

    if index:
      node = index.internalPointer()
      
      if node is None:
        return
      
      dragId = self.getDragId(node.data)
      
      if not dragId:
        return
      
    else:
      return 
    
    model = self.model()
    
    if model.getIcon:
      icon = model.getIcon(node)
      pixmap = icon.pixmap(22,22)
    else:
      icon = 'icons/list-add.png'
      pixmap = QtGui.QPixmap(icon)
    
    pixmap.setMask(pixmap.createHeuristicMask())
 
    
    anchor = QtCore.QPoint(11,11)
    
    itemData = QtCore.QByteArray()
    dataStream = QtCore.QDataStream(itemData, QtCore.QIODevice.WriteOnly)
    dataStream << pixmap << anchor

    mimeData = QtCore.QMimeData()
    mimeData.setText(dragId)
    mimeData.setData(self.mimeType, itemData)

    drag = QtGui.QDrag(self)
    drag.setMimeData(mimeData)
    drag.setPixmap(pixmap)
    drag.setHotSpot(anchor)
    
    drag.exec_(QtCore.Qt.CopyAction | QtCore.Qt.MoveAction, QtCore.Qt.CopyAction)
    
    
class FileSystemTreePanel(QtGui.QWidget, Base):

  def __init__(self, parent, fileTypes=None, callback=None,
               iconProvider=None, showHiddenFiles=False, **kw):
  
    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.dirPath = None
    self.callback = callback
    self.fileTypes = fileTypes or []
    self.model = QtGui.QFileSystemModel()
    self.model.setNameFilterDisables(False)
    if showHiddenFiles:
      filters = self.model.filter()
      self.model.setFilter(filters | QtCore.QDir.Hidden)
    
    if iconProvider:
      self.model.setIconProvider(iconProvider)
    
    row = 0
    frame2 = Frame(self, grid=(row,0))
    texts = [ft.getFilterText() for ft in self.fileTypes]
    self.fileTypePulldown = PulldownList(frame2, texts=texts, objects=self.fileTypes,
                                         callback=self.setFileType, grid=(0,0))
    
    self.findButton = Button(frame2, '', grid=(0,2),
                             callback=self.searchFiles,
                             tipText='Find files',
                             icon=path.join(ICON_DIR,'edit-find.png'))

    self.homeButton = Button(frame2, '', grid=(0,3),
                             callback=self.goHome,
                             tipText='Go to home directiry',
                             icon=path.join(ICON_DIR,'go-home.png'))

    self.rootButton = Button(frame2, '', grid=(0,4),
                             callback=self.goRoot,
                             tipText='Go to root directory',
                             icon=path.join(ICON_DIR,'drive.png'))


    menu = Menu(frame2, callback=self.toggleColumn)
    menu.addItem('Show file size', object=1, icon=None,
                 key=None, checked=True, shortcut=None, tipText=None,
                 index=None)
    menu.addItem('Show file type', object=2, icon=None,
                 key=None, checked=False, shortcut=None, tipText=None,
                 index=None)
    menu.addItem('Show date modified', object=3, icon=None,
                 key=None, checked=False, shortcut=None, tipText=None,
                 index=None)
    self.HIDDEN_FILES_COL = 4
    menu.addItem('Show hidden files', object=self.HIDDEN_FILES_COL, icon=None,
                 key=None, checked=showHiddenFiles, shortcut=None, tipText=None,
                 index=None)
    
    self.configButton = Button(frame2, text='', grid=(0,6),
                             tipText='Configure file tree',
                             icon=path.join(ICON_DIR,'configure.png'))
    
    self.configButton.setMenu(menu)
    
    row += 1                    
    self.treeView = Tree(self, self.model, callback=self._selectPathIndex,
                         doubleCallback=self._callback)
    
    self.showCols = set() # wb104: not sure this is the best name for this variable
    self.treeView.setDragEnabled(True)
    self.treeView.setDragDropMode(self.treeView.DragOnly)
    self.toggleColumn(2)
    self.toggleColumn(3)
    if not showHiddenFiles:
      # do not use self.toggleColumn() here because filter dealt with above
      self.showCols.add(self.HIDDEN_FILES_COL)
    
    self.layout().addWidget(self.treeView, row, 0)
    self.layout().setRowStretch(row, 2)
                                         
    frame2.layout().setSpacing(1)
    frame2.layout().setContentsMargins(1,1,1,1)
    frame2.layout().setColumnStretch(1, 2)
    frame2.layout().setColumnStretch(0, 0)
    
    if fileTypes:
      self.setFileType(fileTypes[0])

    self.model.setRootPath(None)
    self.setShowDirs(True)
  
  def _callback(self, index, haveModKey):
  
    if self.callback:
      filePath = self.model.filePath(index)
      self.callback(filePath, haveModKey)
  
  def _selectPathIndex(self, index):
    
    filePath = self.model.filePath(index)
    
    if not self.model.isDir(index):
      self.dirPath = path.split(filePath)[0]
      
    else:
      self.dirPath = filePath
  
  def toggleColumn(self, col):
    
    if col in self.showCols:
      self.showCols.remove(col)
      if col < self.HIDDEN_FILES_COL:
        self.treeView.setColumnHidden(col, False)
      else:
        filters = self.model.filter()
        self.model.setFilter(filters | QtCore.QDir.Hidden)
   
    else:
      self.showCols.add(col)
      if col < self.HIDDEN_FILES_COL:
        self.treeView.setColumnHidden(col, True)    
      else:
        filters = self.model.filter()
        self.model.setFilter(filters & (~QtCore.QDir.Hidden))

    self.treeView.setColumnWidth(0, 250)
    self.treeView.setColumnWidth(1, 70)
    self.treeView.setColumnWidth(2, 70)
    self.treeView.setColumnWidth(3, 70)    

  def setFileType(self, fileType):
    
    text = fileType.getFilterText()
    self.model.setNameFilters(fileType.extensions)
  
  def setShowDirs(self, showDirs):
  
    if showDirs or not self.dirPath:
      self.treeView.setRootIndex(self.model.index('/'))
    else:
      self.treeView.setRootIndex(self.model.index(self.dirPath))
  
  def goHome(self):
  
    homeDir = self.getHomeDir()
    if homeDir:
      self.openDir(homeDir)
      self.treeView.setRootIndex(self.model.index(self.dirPath))

  def goRoot(self):
  
    rootDir = '/'
    self.openDir(rootDir)
    self.treeView.setRootIndex(self.model.index(rootDir))

  def getHomeDir(self):
  
    return os.environ.get('HOME') or os.environ.get('HOMEPATH')
    
  def searchFiles(self):
    
    from gui.qtgui.InputDialog import askString
    from fnmatch import fnmatch
    from os import path, walk
    
    if not self.dirPath:
      msg = 'No directory selected to search within'
      showWarning('Cannot continue', msg, parent=self)
      return
      
    text = askString('Search for files', 'Names matching (Use "*", "?" widcards):', '*', parent=self)
    if not text:
      return
    
    self.window().setCursor(QtCore.Qt.WaitCursor)
    
    dirPath = self.dirPath
    filePaths = []
    dirPaths = set()
    for parentDir, subFolders, fileNames in walk(str(dirPath)):
      for fileName in fileNames:
        if fnmatch(fileName, text):
          filePaths.append(path.join(parentDir, fileName))
          dirPaths.add(parentDir)
   
    if not filePaths:
      msg = 'No files found matching "%s"' % text
      showWarning('Failure', msg, parent=self)
      self.window().unsetCursor()
      return
          
    uniqPaths = set(dirPaths)
    for dirPath in dirPaths:
      
      prevPath = None
      while dirPath != prevPath:
        prevPath = dirPath
        dirPath, null = path.split(dirPath)
        uniqPaths.add(dirPath)
   
    self.treeView.collapseAll()
    #self.model.setNameFilters([text,])
    self.model.setRootPath(dirPath)
    expand = self.treeView.expand
    getIndex = self.model.index
    for dirPath in uniqPaths:
      expand(getIndex(dirPath))
    
    indices = [getIndex(filePath) for filePath in filePaths]    
    self.treeView.setSelected(indices)
    self.treeView.scrollTo(indices[0])
    self.window().unsetCursor()

  def openDir(self, dirPath):
    
    self.dirPath = dirPath
    self.treeView.collapseAll()
    
    home = self.getHomeDir()
    if dirPath and dirPath.startswith(home):
      self.treeView.setRootIndex(self.model.index(home))
    else:
      self.treeView.setRootIndex(self.model.index('/'))
        
    prevPath = None
    while dirPath != prevPath:
      self.treeView.expand(self.model.index(dirPath))
      prevPath = dirPath
      dirPath, null = path.split(dirPath)
    
    self.treeView.scrollTo(self.model.index(self.dirPath))

if __name__ == '__main__':

  import sys
  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  
  def callback(filePath):
    print(filePath)
    
  argv = sys.argv
  app = Application(argv)
  popup = BasePopup(title='Test File Tree')
  
  FileSystemTreePanel(popup, callback=callback)
  
  app.start()
  
