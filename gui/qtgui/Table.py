import re

from numpy import int8, int16, int32, int64, uint8, uint16, uint32, uint64
from numpy import float32, float64

from PySide import QtCore, QtGui

Qt = QtCore.Qt
QPainter = QtGui.QPainter
QColor = QtGui.QColor
QRect = QtCore.QRect

from datetime import datetime

import io

from gui.qtgui.Base import Base
from gui.qtgui.BasePopup import BasePopup
from gui.qtgui.ButtonList import ButtonList, UtilityButtonList
from gui.qtgui.CheckButton import CheckButton
from gui.qtgui.Entry import Entry
from gui.qtgui.FileSelect import selectSaveFile, FileType
from gui.qtgui.Graph import Graph, GraphAxis, GraphDataSet
from gui.qtgui.Label import Label
from gui.qtgui.Menu import Menu
from gui.qtgui.MultiWidget import MultiWidget
from gui.qtgui.Print import PrintDialog
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.RadioButtons import RadioButtons
from gui.qtgui.Splitter import Splitter
from gui.qtgui.TableModel import ObjectTableModel, HORIZONTAL

NUMPY_INT_TYPES = (int8, int16, int32, int64, uint8, uint16, uint32, uint64)
NUMPY_FLOAT_TYPES = (float32, float64)
INT_TYPES = (int,) + NUMPY_INT_TYPES
FLOAT_TYPES = (float,) + NUMPY_FLOAT_TYPES
NUMBER_TYPES = INT_TYPES + FLOAT_TYPES + (datetime,)
STRING_TYPES = (type(''), type(u''))
NUMBER_RE = re.compile('\s*(-?\d+\.*\d*|\d*\.*\d+)')

NULL_INDEX = QtCore.QModelIndex()

# ? Cleanup menu slots after use 
# ? Need text auto format ? - Maybe not
# ? Need separate highlight funcs ? - Maybe not
# ? Fill whitespace

# ? Scroll focus follows selected object
# Justification

class Table(QtGui.QTableView, Base):

  def __init__(self, parent=None, model=None, **kw):

    QtGui.QTableView.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    self.model = model

    self.resizeColumnsToContents()

    if model:
      self.setModel(model)
    self.setSortingEnabled(True)
    self.sortByColumn(0, QtCore.Qt.AscendingOrder)
  
  def setHorizontalHeaderVisible(self, isVisible):

    self.horizontalHeader().setVisible(isVisible)

  def setVerticalHeaderVisible(self, isVisible):

    self.verticalHeader().setVisible(isVisible)

  def setAllColumnResize(self, canResize=True):

    for column in self.model.numberCols():
      self.setColumnResize(column, canResize)

  def setColumnResize(self, column, canResize=True):
    
    header = self.horizontalHeader()
    
    if canResize:
      mode = header.ResizeToContents
      
    else:
      mode = header.Fixed
    
    header.setResizeMode(column, mode)

  ############################################################
  # other functions, which are inherited
  ############################################################

  # setFont(font)
  # setSortingEnabled(isEnabled)
  # setShowGrid(isShown)
  # setColumnWidth(column, width)
  # setRowHeight(row, height)
  # setMinimumSize(width, height)

BG_COLOR = QColor('#E0E0E0')

ALIGN_OPTS = {'C':Qt.AlignHCenter,
              'c':Qt.AlignHCenter,
              'Center':Qt.AlignHCenter,
              'center':Qt.AlignHCenter,
              'centre':Qt.AlignHCenter,
              'Centre':Qt.AlignHCenter,
              'L':Qt.AlignLeft,
              'l':Qt.AlignLeft,
              'left':Qt.AlignLeft,
              'Left':Qt.AlignLeft,
              'R':Qt.AlignRight,
              'r':Qt.AlignRight,
              'Right':Qt.AlignRight,
              'right':Qt.AlignRight,
              }


class Column:

  def __init__(self, heading, getValue, getEditValue=None, setEditValue=None,
               editClass=None, editArgs=None, editKw=None, tipText=None,
               getColor=None, getIcon=None, stretch=False, format=None,
               editDecimals=None, editStep=None, editMax=None, editMin=None,
               alignment=None, orderFunc=None):
    
    self.heading = heading
    self.getValue = getValue or self._defaultText
    self.getEditValue = getEditValue or getValue
    self.setEditValue = setEditValue
    self.editClass = editClass
    self.editArgs = editArgs or []
    self.editKw = editKw or {}
    self.stretch = stretch
    self.format = format
    self.editDecimals = editDecimals
    self.editStep = editStep
    self.editMin = editMin
    self.editMax = editMax
    self.defaultIcon = None
    #self.alignment = ALIGN_OPTS.get(alignment, alignment) | Qt.AlignVCenter
    # Alignment combinations broken in PySide v1.1.1
    # Use better default than top left
    self.alignment = alignment or Qt.AlignCenter
    self.orderFunc = orderFunc
    
    self.getIcon = getIcon or self._defaultIcon
    self.getColor = getColor or self._defaultColor
    self.tipText = tipText
    
    self._checkTextAttrs()
  
  def getFormatValue(self, obj):
    
    value = self.getValue(obj)
    format = self.format
    if isinstance(value, datetime):
      if not format:
        format = '%c'
      return value.strftime(format)
    elif format and (value is not None):
      try:
        return format % value
      except TypeError:
        return format.format(value)
      
    else:
      return value
    
  def _checkTextAttrs(self):
    
    if isinstance(self.getValue, str):
      attr = self.getValue
      self.getValue = lambda obj: getattr(obj,attr)

    if isinstance(self.getEditValue, str):
      attr = self.getEditValue
      self.getEditValue = lambda obj: getattr(obj,attr)    

    if isinstance(self.setEditValue, str):
      attr = self.setEditValue
      self.setEditValue = lambda obj, value: setattr(obj,attr,value)

    if isinstance(self.getIcon, QtGui.QIcon):
      self.defaultIcon = self.getIcon
      self.getIcon = self._defaultIcon
  
  def _defaultText(self, obj):
  
    return None
  
  def _defaultColor(self, obj):
    
    return BG_COLOR

  def _defaultIcon(self, obj):

    return self.defaultIcon

EDIT_ROLE = QtCore.Qt.EditRole

class ObjectTableItemDelegate(QtGui.QStyledItemDelegate):

  def __init__(self, parent):
    
    QtGui.QStyledItemDelegate.__init__(self, parent)
    self.customWidget = None
    self.parent = parent

  def createEditor(self, parentWidget, itemStyle, index): # returns the edit widget
    
    col = index.column()
    objCol = self.parent.columns[col]
    
    if objCol.editClass:
      widget = objCol.editClass(None, *objCol.editArgs, **objCol.editKw)
      widget.setParent(parentWidget)
      self.customWidget = True
      return widget
    
    else:
      obj = self.parent.objects[index.row()]
      editValue = objCol.getEditValue(obj)
    
      if isinstance(editValue, (list, tuple)):
        widget = PulldownList(None)
        widget.setParent(parentWidget)
        self.customWidget = True
        return widget
      
      else:
        # Use the default, type-dependant factory
        # Deals with strings, bools, date time etc.
        self.customWidget = None
        editor = QtGui.QStyledItemDelegate.createEditor(self, parentWidget, itemStyle, index)
        
        if isinstance(editor, QtGui.QDoubleSpinBox):
          numDecimals = objCol.editDecimals
          
          if numDecimals is not None:
            editor.setDecimals(numDecimals)
            
            if objCol.editStep:
              editor.setSingleStep(objCol.editStep)
            else:
              editor.setSingleStep(10**-numDecimals)
              
          elif objCol.editStep:
            editor.setSingleStep(objCol.editStep)
          
          if objCol.editMin is not None:
            editor.setMinimum(objCol.editMin)
            
          if objCol.editMax is not None:
            editor.setMaximum(objCol.editMax)
          
        if isinstance(editor, QtGui.QSpinBox):
          if objCol.editStep:
            editor.setSingleStep(objCol.editStep)
        
        return editor
      
    
  def setEditorData(self, widget, index): # provides the widget with data
    
    if self.customWidget:
      model = index.model()
      value = model.data(index, EDIT_ROLE)
      
      if not isinstance(value, (list, tuple)):
        value = (value,)
      
      if hasattr(widget, 'setColor'):
        widget.setColor(*value)
        
      elif hasattr(widget, 'setData'):
        widget.setData(*value)
      
      elif hasattr(widget, 'set'):
        widget.set(*value)

      elif hasattr(widget, 'setValue'):
        widget.setValue(*value)

      elif hasattr(widget, 'setFile'):
        widget.setFile(*value)
      
      else:
        msg = 'Widget %s does not expose "setData", "set" or "setValue" method; ' % widget
        msg += 'required for table proxy editing'
        raise Exception(msg)
        
      
    else:
      return QtGui.QStyledItemDelegate.setEditorData(self, widget, index)
       
  def updateEditorGeometry(self, widget, itemStyle, index):# ensures that the editor is displayed correctly 
    
    if self.customWidget:
      cellRect = itemStyle.rect
      x = cellRect.x()
      y = cellRect.y()
      hint = widget.sizeHint()
      
      if hint.height() > cellRect.height():
        if isinstance(widget, QtGui.QComboBox): # has a popup anyway
          widget.move(cellRect.topLeft())
          
        else:
          pos = widget.mapToGlobal(cellRect.topLeft())
          widget.setParent(self.parent, QtCore.Qt.Popup) # popup so not confined
          widget.move(pos)
     
      else:
        width = max(hint.width(), cellRect.width())
        height = max(hint.height(), cellRect.height())
        widget.setGeometry(x, y, width, height)
        
     
    else:
      return QtGui.QStyledItemDelegate.updateEditorGeometry(self, widget, itemStyle, index)
      
  def setModelData(self, widget, mode, index):#returns updated data 
     
    if self.customWidget:
    
      if hasattr(widget, 'get'):
        value = widget.get()
        
      elif hasattr(widget, 'value'):
        value = widget.value()
        
      elif hasattr(widget, 'getFile'):
        files = widget.selectedFiles()
        
        if not files:
          return
          
        value = files[0]
       
      else:
        msg = 'Widget %s does not expose "get" or "value" method; ' % widget
        msg += 'required for table proxy editing'
        raise Exception(msg)
      
      del widget
      model = index.model()
      model.setData(index, value, EDIT_ROLE)
     
    else:
      return QtGui.QStyledItemDelegate.setModelData(self, widget, mode, index)
      
class ObjectTableViewport(QtGui.QWidget):      
  
  def __init__(self, parent):
    
    QtGui.QWidget.__init__(self, parent=parent)

  """ 
    self.setAutoFillBackground(True)
    self.parent = parent
  
  def paintEvent(self, event):
  
    header = self.parent.verticalHeader()
    pointA = header.rect().bottomLeft()
    pointB = self.rect().bottomRight()
    
    painter = QPainter(self)
    painter.setBrush(BG_COLOR)
    painter.drawRect(QRect(pointA, pointB))
    painter.end()
    
    
    return QtGui.QWidget.paintEvent(self, event)
  """   

LESS_THAN = QtGui.QSortFilterProxyModel.lessThan

class ObjectTableProxyModel(QtGui.QSortFilterProxyModel):      
  
  def __init__(self, parent):
    
    QtGui.QSortFilterProxyModel.__init__(self, parent=parent)
    self.table = parent
    
  def lessThan(self, leftIndex, rightIndex):
    
    table = self.table
    columns = table.columns
    objCol = columns[leftIndex.column()]
    
    if objCol.orderFunc:

      return objCol.orderFunc(objCol.getValue(objA), objCol.getValue(objB))
    
    elif objCol.format:
      objects = table.objects
      objA = objects[leftIndex.row()]
      objB = objects[rightIndex.row()]
      return objCol.getValue(objA) < objCol.getValue(objB)
      
    else:
      return LESS_THAN(self, leftIndex, rightIndex)

class ObjectHeaderView(QtGui.QHeaderView):

  def __init__(self, orient, parent):
  
    QtGui.QHeaderView.__init__(self, orient, parent)
    self.table = parent

  #def sizeHint(self):
  
  #  return QtCore.QSize(30*len(self.table.columns), self.table.bbox('A').height())

  #def minimumSize(self):
  #
  #  return QtCore.QSize(30*len(self.table.columns), self.table.bbox('A').height())
    
class ObjectTable(QtGui.QTableView, Base):

  def __init__(self, parent, columns, objects=None, callback=None,
               multiSelect=True, selectRows=True, numberRows=False,
               sortCol=0, **kw):
    
    QtGui.QTableView.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    styleSheet = """
    QHeaderView { font-weight:normal }
    QHeaderView::section {padding:2px;}
    QHeaderView::down-arrow { top: 10px; right: 1px; image:url(icons/sort-up.png) }
    QHeaderView::up-arrow { top: 10px; right: 1px; image:url(icons/sort-down.png) }
    """
    
    self.setStyleSheet(styleSheet)
    self.graphPanel = None
    self.filterPanel = None
    self.model = None
    self.columns = columns
    if objects:
      self.objects = list(objects)
    else:
      self.objects = []
    self.callback = callback
    self.fontMetric = QtGui.QFontMetricsF(self.font())
    self.bbox = self.fontMetric.boundingRect
    self._silenceCallback = False
    self.selectRows = selectRows
    
    self.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
    self.setHorizontalScrollMode(self.ScrollPerItem)
    self.setVerticalScrollMode(self.ScrollPerItem)
    self.setSortingEnabled(True)
    self.setAutoFillBackground(True)
    
    #self.setSizePolicy(QtGui.QSizePolicy.Preferred, 
    #                   QtGui.QSizePolicy.Preferred)
    
    if multiSelect:
      self.setSelectionMode(self.ExtendedSelection)
      # + Continuous etc possible
    else:
      self.setSelectionMode(self.SingleSelection)
    
    if selectRows:
      self.setSelectionBehavior(self.SelectRows)
    else:
      self.setSelectionBehavior(self.SelectItems)
      # + Columns possible
    
    self._setupModel()
    
    if sortCol is not None:
      self.sortByColumn(sortCol, QtCore.Qt.AscendingOrder)
    
    delegate = ObjectTableItemDelegate(self)
    self.setItemDelegate(delegate)
    
    model = self.selectionModel()
    model.selectionChanged.connect(self._callback)
    #model.currentRowChanged.connect(self._callback)
    
    header = self.verticalHeader()
    header.setResizeMode(header.Interactive)
    header.setStretchLastSection(False)
    
    rowHeight = self.bbox('A').height() + 4
    header.setMinimumSectionSize(rowHeight)
    header.setDefaultSectionSize(rowHeight)
    
    if numberRows:
      header.setVisible(True)
    else:
      header.setVisible(False)
  
    #header = ObjectHeaderView(Qt.Horizontal, self)
    header = self.horizontalHeader()
    header.setMovable(True)
    header.setMinimumSectionSize(30)
    header.setDefaultSectionSize(30)
    #header.setSortIndicatorShown(False)
    #header.setStyleSheet('QHeaderView::down-arrow { image: url(icons/sort-up.png);} QHeaderView::up-arrow { image: url(icons/sort-down.png);}')
    
    self.setupHeaderStretch()
  
  def _setupModel(self):
  
    
    if self.model:
      sortDetails = (self.model.sortColumn(), self.model.sortOrder())
      
    else:
      sortDetails = None  
    
    objModel = ObjectTableModel(self)
    
    model = ObjectTableProxyModel(self)
    model.setSourceModel(objModel)
    
    if sortDetails is not None:
      column, order = sortDetails
      model.sort(column, order)
    
    self.setModel(model)
    self.model = model
    
    return model
    
  def resizeEvent(self, event):
    
    if self.graphPanel and self.graphPanel.isVisible():
      pos = self.graphPanel.pos()
      x = pos.x()
      y = pos.y()
      w = self.width()-x
      h = self.height()
      self.graphPanel.setGeometry(x, y, w, h)
      
    # If the table is connected to a qtgui Splitter it
    # should not resize unless it is specifically asked to.
    # This helps avoiding infinite repaint loops.
    if not (isinstance(self.parent, Splitter) or self.parent.__class__.__name__ == Splitter.__name__) or \
           self.parent.doResize == True:
      return QtGui.QTableView.resizeEvent(self, event)
  
  def clearSelection(self):
  
    model = self.selectionModel()
    model.clear()
  
  def _callback(self, itemSelection):
    
    if self._silenceCallback:
      return
        
    if self.graphPanel and self.graphPanel.isVisible():
      graph = self.graphPanel.graph
      graph.coordsOff()
      rows = self.getSelectedRows()
      vLines = []
      hLines = []
      
      for row in rows:
        for dataSet in graph.dataSets:
          x, y = dataSet.dataPoints[row][:2]
          vLines.append(x)
          hLines.append(y)
      
      graph.drawVerticalLines(vLines)
      graph.drawHorizontalLines(hLines)
      
    elif self.callback:
      index = self.getCurrentIndex()
      row = index.row()
      col = index.column()
       
      model = self.selectionModel()
      
      if self.selectRows:
        selection = model.selectedRows(column=0)
      else:
        selection = model.selectedIndexes()
      
      if selection:
        rows = [i.row() for i in selection]
        rows.sort()
        if row not in rows:
          row = rows[0]
 
        index = self.model.index(row, 0)
        row = self.model.mapToSource(index).row()
 
        if row >= 0:
          obj = self.objects[row]
          self.callback(obj, row, col)
      
      else:
        self.callback(None, row, col)
  
  def getCurrentIndex(self):
    
    model = self.selectionModel()
    index = model.currentIndex()
    index = self.model.mapToSource(index)
    
    return index
      
  def getCurrentObject(self):
    
    selectionModel = self.selectionModel()
    index = selectionModel.currentIndex()
    index = self.model.mapToSource(index)
    
    row = index.row()
    
    if row < 0:
      return
    else:
      return self.objects[row]

  def getCurrentRow(self):
    
    selectionModel = self.selectionModel()
    index = selectionModel.currentIndex()
    index = self.model.mapToSource(index)
    
    return index.row()
            
  def getSelectedRows(self):
    
    model = self.selectionModel()
    
    if self.selectRows:
      selection = model.selectedRows(column=0)
    else:
      selection = model.selectedIndexes()
    
    rows = [i.row() for i in selection]
    #rows = list(set(rows))
    #rows.sort()
    
    return rows
  
  def getSelectedObjects(self):
     
    model = self.selectionModel()
    if self.selectRows:
      selection = model.selectedRows(column=0)
    else:
      selection = model.selectedIndexes()
      
    objects = self.objects
    selectedObjects = []
    
    for index in selection:
      row = self.model.mapToSource(index).row()
      selectedObjects.append(objects[row])
    
    return selectedObjects
  
  def setCurrentRow(self, row):
    
    selectionModel = self.selectionModel()
    index = self.model.index(row, 0)
    self._silenceCallback = True
    selectionModel.clearSelection()
    self._silenceCallback = False
    selectionModel.select(index, selectionModel.Select | selectionModel.Rows)
    self.setFocus(Qt.OtherFocusReason)

  def setCurrentObject(self, obj):
    
    if obj in self.objects:
      row = self.objects.index(obj)
      selectionModel = self.selectionModel()
      index = self.model.sourceModel().index(row, 0)
      index = self.model.mapFromSource(index)
      
      self._silenceCallback = True
      selectionModel.clearSelection()
      self._silenceCallback = False
      selectionModel.select(index, selectionModel.Select | selectionModel.Rows)
      self.setFocus(Qt.OtherFocusReason)
 
  def selectRow(self, row):
  
    selectionModel = self.selectionModel()
    index = self.model.index(row, 0)
    self._silenceCallback = True
    selectionModel.clearSelection()
    self._silenceCallback = False
    selectionModel.setCurrentIndex(index, selectionModel.SelectCurrent | selectionModel.Rows)
    self.setFocus(Qt.OtherFocusReason)
 
  def selectObject(self, obj):
  
    if obj in self.objects:
      row = self.objects.index(obj)
      selectionModel = self.selectionModel()
      index = self.model.sourceModel().index(row, 0)
      index = self.model.mapFromSource(index)
      
      self._silenceCallback = True
      selectionModel.clearSelection()
      self._silenceCallback = False
      selectionModel.setCurrentIndex(index, selectionModel.SelectCurrent | selectionModel.Rows)
      self.setFocus(Qt.OtherFocusReason)
 
  def setSelectedObjects(self, selection):
  
    return self.setCurrentObjects(selection)
 
  def setCurrentObjects(self, selection):
    
    objects = self.objects
    selectionModel = self.selectionModel()

    rows = []
    uniqObjs = set(selection)
    
    for row, obj in enumerate(objects):
      if obj in uniqObjs:
        rows.append(row)
    
    if rows:
      self._silenceCallback = True
      selectionModel.clearSelection()
      self.setUpdatesEnabled(False)
      
      selectMode = selectionModel.Select | selectionModel.Rows
      select = selectionModel.select
      getIndex = self.model.sourceModel().index
      mapFromSource = self.model.mapFromSource

      for row in rows:
        index = getIndex(row, 0)
        index = mapFromSource(index)
        select(index, selectMode)
      
      self._silenceCallback = False
      self._callback(None)
      self.setUpdatesEnabled(True)

    self.setFocus(Qt.OtherFocusReason)

  
  def contextMenuEvent(self, event):

    mouseMenu = Menu(None, 'Options')
    
    # Filters
    
    mouseMenu.addItem('Filter', self.filterRows)
    
    # Export text
    
    mouseMenu.addItem('Export Text', self.exportText)
    
    # Export image
    
    imageMenu = Menu(mouseMenu, 'Export Image')
    
    imageMenu.addItem('JPEG', self.exportJpeg)
    imageMenu.addItem('PNG', self.exportPng)
    imageMenu.addItem('PDF', self.exportPdf)
    imageMenu.addItem('SVG', self.exportSvg)
    
    # Graph menu
    
    graphMenu = Menu(mouseMenu, 'Graph')
    options = [('Row Number',None)]
    columns = self.columns
    objects = self.objects
    reMatch = re.match
    
    for i, column in enumerate(columns):
      
      for obj in objects:
        value = column.getValue(obj)
        
        if value is not None:
          if isinstance(value, bool):
            break
          
          elif isinstance(value, NUMBER_TYPES):
            heading = column.heading.replace('\n',' ')
            options.append((heading, i))
            break
            
          elif isinstance(value, STRING_TYPES):
            match = reMatch(NUMBER_RE,value)
            if match:
              heading = column.heading.replace('\n',' ')
              options.append((heading, i))
              break
              
    for headingA, i in options:
      menu = Menu(graphMenu, 'X:%s' % (headingA))
      
      for headingB, j in options:
        if i == j:
          continue
        
        menu.addItem('Y:%s' % (headingB), lambda x=i, y=j: self.graphColumns(x, y))
    
    # Printing
    
    printMenu = Menu(mouseMenu, 'Print')
    printMenu.addItem('Pixmap Image', self.printPixmap)
    
    # Table info
    
    infoMenu = Menu(mouseMenu, 'Table info')
    row = self.rowAt(event.y())+1
    col = self.columnAt(event.x())+1
    infoMenu.addItem('Total rows: %d' % len(self.objects))
    infoMenu.addItem('Selected rows: %d' %  len(self.getSelectedRows()))
    infoMenu.addItem('Current row, col: %d,%d' % (row, col))
    
    mouseMenu.exec_(event.globalPos())
  
  def setupHeaderStretch(self):
    
    columns = self.columns
    header = self.horizontalHeader()
    stretch = False
    objects = self.objects
    setColumnWidth = self.setColumnWidth
    bbox = self.bbox
    
    for i, column in enumerate(columns):
      if column.stretch:
        header.setResizeMode(i, header.Stretch)
        stretch = True
        
    if objects:
      minSize = bbox('MM').width()
      colSizes = [max(minSize, bbox(col.heading).width()) for col in columns]
      for i, objCol in enumerate(columns):
        for obj in objects[:35] + objects[-5:]:
          value = objCol.getFormatValue(obj)
          
          if isinstance(value, float):
            value = '%.5f' % (value,)
          else:
            value = '%s' % (value,)
          
          size = bbox(value).width()
 
          if size > colSizes[i]:
            colSizes[i] = size
 
      for i, width in enumerate(colSizes):
        setColumnWidth(i, width+12)
      
    if not stretch:
      header.setStretchLastSection(True)
 
  def closeEvent(self, event):
    
    self.hideGraph()
    self.hideFilter()
    QtGui.QTableView.closeEvent(self, event)

  def destroy(self, *args):
    
    if self.filterPanel:
      self.filterPanel.destroy()
    
    if self.graphPanel:
      self.graphPanel.destroy()  
      
    QtGui.QTableView.destroy(self, *args)
  
  def exportJpeg(self):
   
    fileTypes = [FileType('JPEG', ['*.jpg','.jpeg','.jpr']),]
    filePath = selectSaveFile(self, caption='Select image file name',
                              directory=None, fileTypes=fileTypes)
                   
    if filePath:
      widget = self
      pixmap = QtGui.QPixmap.grabWidget(widget, widget.rect())
      pixmap.save(filePath, 'JPEG')
    
  def exportPng(self):
  
    fileTypes = [FileType('PNG', ['*.png']),]
    filePath = selectSaveFile(self, caption='Select image file name',
                              directory=None, fileTypes=fileTypes)
                   
    if filePath:
      widget = self
      pixmap = QtGui.QPixmap.grabWidget(widget, widget.rect())
      pixmap.save(filePath, 'PNG')
    
  def exportText(self):
  
    popup = ObjectTableExport(self)
    popup.show()
    popup.move(self.window().pos())
  
  def filterRows(self):
   
    self.hideGraph()
 
    if not self.filterPanel:
      self.filterPanel = ObjectTableFilter(self)
 
    self.filterPanel.move(self.window().pos())
    self.filterPanel.show()
  
  def graphColumns(self, xCol, yCol):
   
    self.hideFilter()
 
    header = self.horizontalHeader()
    header.setStretchLastSection(False)
    header.setResizeMode(header.Fixed)
    x = self.verticalHeader().width()
 
    otherCols = set(range(len(self.columns)))
    if xCol is not None:
      otherCols.remove(xCol)
      self.setColumnHidden(xCol, False)
      x += header.sectionSize(xCol)
 
    if yCol is not None:
      otherCols.remove(yCol)
      self.setColumnHidden(yCol, False)
      x += header.sectionSize(yCol)
 
    for i in otherCols:
      self.setColumnHidden(i, True)
 
    w = self.width() - x
    y = 0
    h = self.height()
 
    if self.graphPanel:
      self.graphPanel.updateGraph(xCol, yCol)

    else:
      self.graphPanel = ObjectTableGraph(self, xCol, yCol)
 
    self.graphPanel.setGeometry(x, y, w, h)
    self.graphPanel.show()
 
  def hideGraph(self):
   
    if self.graphPanel:
      self.graphPanel.hide()
      for i in range(len(self.columns)):
        self.setColumnHidden(i, False)
 
      self.setupHeaderStretch()
  
  def hideFilter(self):
     
    if self.filterPanel:
      self.filterPanel.close()
  
  def unfilter(self):
  
    if self.filterPanel:
      self.filterPanel.unfilterTable()
    
  def printPixmap(self):
  
    dialog = PrintDialog(self) 
    dialog.printWidgetPixmap(self)
 
  def getObjects(self, filtered=False):
  
    if self.filterPanel:
      if filtered:
        return list(self.objects)

      else:
        return list(self.filterPanel.origObjects)
      
    else:
      return list(self.objects)    
  
  def getSortedObjects(self):
  
    model = self.model # .filterModel()
    n = len(self.objects)
    
    rows = []
    for i in range(n):
      indexA = model.index(i, 0) 
      index = model.mapToSource(indexA)# original, underlying
      rows.append(index.row())
      
    return [self.objects[i] for i in rows]
      
  
  def _syncFilterObjects(self, applyFilter=False):
    
    if self.filterPanel:
      if applyFilter is not None:
        self.filterPanel.origObjects = self.objects

      status = self.filterPanel.status
      if applyFilter and (status is not None):
        self.filterPanel.filterTable(status)
   
  def replaceCurrentObject(self, newObj):
    
    model = self.model
    selectionModel = self.selectionModel()
    row = selectionModel.currentIndex().row() # the visible row
    indexA = model.index(row, 0)
    indexB = model.index(row, len(self.columns)-1)
    index = model.mapToSource(indexA) # the underlying, original
    
    if index.row() < 0:
      return
    
    self.objects[index.row()] = newObj
    self._syncFilterObjects()
 
    model.dataChanged.emit(indexA,indexB)    

  
  def setObject(self, i, newObj):
    """Replaces an object in the _underlying_ list"""
    
    model = self.model
    self.objects[i] = newObj
    self._syncFilterObjects()
    
    index = model.sourceModel().index(i, 0) # the underlying, original
    indexA = model.mapFromSource(index) # the visible row
    indexB = model.index(indexA.row(), len(self.columns)-1)

    model.dataChanged.emit(indexA,indexB)
  
  
  def setRowObject(self, row, newObj):
    """Replaces an object in the _visible_ rows"""
    
    model = self.model
    indexA = model.index(row, 0) # the visible row
    indexB = model.index(row, len(self.columns)-1) # the visible row
    index = model.mapToSource(indexA) # the underlying, original
    
    self.objects[index.row()] = newObj
    self._syncFilterObjects()

    model.dataChanged.emit(indexA,indexB)
   
  
  def setObjectsAndColumns(self, objects, columns):
     
    self.setObjects([])
    self.setColumns(columns)
    self.setObjects(objects)
    
   
  def setObjects(self, objects, applyFilter=False):
    
    model = self.model
    sourceModel = model.sourceModel()
    selected = set(self.getSelectedObjects())
    current = self.getCurrentObject()
    
    n = len(objects)
    m = len(self.objects)
    c = len(self.columns) 
    
    self.objects = list(objects)
    sourceModel.objects = self.objects
    
    if n > m:
      sourceModel.beginInsertRows(QtCore.QModelIndex(), m, n-1)
    
    elif m > n:
      sourceModel.beginRemoveRows(QtCore.QModelIndex(), n, m-1)

    indexA = model.index(0, 0)
    indexB = model.index(n-1, c-1)
    model.dataChanged.emit(indexA,indexB) # the visible rows
    
    if n > m:
      sourceModel.endInsertRows()     
    
    elif m > n:
      sourceModel.endRemoveRows()
    
    if selected:
      self._silenceCallback = True
      selectionModel = self.selectionModel()
      selectionModel.clearSelection()
      selectMode = selectionModel.Select | selectionModel.Rows
      select = selectionModel.select
      getIndex = self.model.sourceModel().index
      mapFromSource = self.model.mapFromSource
      
      for row, obj in enumerate(objects):
        if obj in selected:
          index = getIndex(row, 0)
          index = mapFromSource(index)
 
          if obj is current:
            selectionModel.setCurrentIndex(index, selectionModel.SelectCurrent)
 
          select(index, selectMode)
      
      self._silenceCallback = False
          
    sortCol = model.sortColumn()
    sortOrder = model.sortOrder()
    
    if sortCol >= 0:
      model.sort(sortCol, sortOrder)

    self._syncFilterObjects(applyFilter)
    self.setupHeaderStretch()
 
  def setColumns(self, columns):
    
    model = self.model
    sourceModel = model.sourceModel()
    selected = self.getSelectedObjects()
    current = self.getCurrentObject()
    
    n = len(columns)
    m = len(self.columns)
    r = len(self.objects) 
       
    self.columns = columns
    sourceModel.columns = columns
    
    if n > m:
      sourceModel.beginInsertColumns(QtCore.QModelIndex(), m, n-1)
    
    elif m > n:
      sourceModel.beginRemoveColumns(QtCore.QModelIndex(), n, m-1)

    indexA = model.index(0, 0)
    indexB = model.index(r-1, n-1)
    model.dataChanged.emit(indexA,indexB) # the visible rows
    model.headerDataChanged.emit(HORIZONTAL, 0,n-1)
    
    if n > m:
      sourceModel.endInsertColumns()     
    
    elif m > n:
      sourceModel.endRemoveColumns()
     
    if selected:
      # may need to now highlight more columns
      selectionModel = self.selectionModel()
      for obj in selected:
        row = self.objects.index(obj)
        index = model.index(row, 0)
        
        if obj is current:
          selectionModel.select(index, selectionModel.SelectCurrent | selectionModel.Rows)
          selectionModel.setCurrentIndex(index, selectionModel.SelectCurrent | selectionModel.Rows)
 
        else:
          selectionModel.select(index, selectionModel.Select | selectionModel.Rows)

    self.setupHeaderStretch()

  def getObject(self, row):
  
    return self.objects[row]

SEARCH_MODES = [ 'Literal','Case Sensitive Literal','Regular Expression' ]

class ObjectTableFilter(BasePopup):

  def __init__(self, table):
  
    BasePopup.__init__(self, title='Filter Table')
    
    self.table = table
    self.status = None
    self.origObjects = table.objects

    label = Label(self, 'Filter Column', grid=(0,0))
    
    columns = table.columns
    
    texts = [c.heading for c in columns]
    objects = range(len(columns))
    
    tIndex = table.getCurrentIndex()
    if tIndex is None:
      index = 0
    else:
      index = tIndex.column()  
    
    self.colPulldown = PulldownList(self, texts, objects, index=index, grid=(0,1))

    label = Label(self, 'Objects to filter', grid=(0,2))
    
    self.filterObjRadio = RadioButtons(self, ['All','Selected'], 0, grid=(0,3))

    self.entry = Entry(self, grid=(1,0), gridSpan=(1,4))

    self.filterModeRadio = RadioButtons(self, SEARCH_MODES, 0, grid=(2,0), gridSpan=(1,4))
    
    
    texts = ['Reset','Filter\nInclude','Filter\nExclude',None]
    callbacks = [self.unfilterTable, self.filterInclude,
                 self.filterExclude, self.close]
    icons = ['icons/edit-undo.png', None, None, 'icons/window-close.png']
    buttons = ButtonList(self, texts, callbacks, icons, grid=(3,0), gridSpan=(1,4))    
    
    self.setWindowFlags(QtCore.Qt.Tool)
    self.setSize(300,100)
  
  def close(self):
  
    BasePopup.close(self)
    
  def unfilterTable(self):
  
    self.table.setObjects(self.origObjects)
    self.status = None
  
  def filterInclude(self, *event):
  
    self.filterTable(True)
  
  def filterExclude(self, *event):
  
    self.filterTable(False)
    
  def filterTable(self, includeMatches=True):
  
    if not self.origObjects:
      self.status = None
      return
    
    string = self.entry.get()
    if not string:
      self.status = None
      return
    
    self.status = includeMatches
    columns = self.table.columns
    objCol = columns[self.colPulldown.currentObject()]
    mode = self.filterModeRadio.get()
    flag = re.S
    
    def exclude(a,b,c):
      return not re.search(a,b,c)
    
    if includeMatches:
      find = re.search
    else:
      find = exclude
    
    if mode != SEARCH_MODES[2]:
      string = re.escape(string)
    
    if mode == SEARCH_MODES[0]:
      flag = re.I
    
    objects = []
    objectsAppend = objects.append
    
    if self.filterObjRadio.getIndex() == 0:
      filterObjs = self.origObjects
    else:
      filterObjs = self.table.getSelectedObjects()
    
    for  obj in filterObjs:
      value = u'%s' % (objCol.getValue(obj))
      match = find(string, value, flag)
      
      if match:
        objectsAppend(obj)
    
    self.table.clearSelection()    
    self.table.setObjects(objects, None)
  
class ObjectTableGraph(QtGui.QWidget):

  def __init__(self, table, xCol, yCol):
    
    QtGui.QWidget.__init__(self, parent=table)
    
    self.table = table
    self.xCol = xCol
    self.yCol = yCol
                       
    xAxis = GraphAxis('xName', labels=None, ticks=True)
    yAxis = GraphAxis('yName', labels=None, ticks=True)
   
    axes = (xAxis, yAxis)
    self.graph = graph = Graph(self, axes, [], title='Table Columns Graph',
                               callback=None, grid=(0,0))
                           
    texts = [None, None,'Fit',
             'Line','Scatter','Histogram']
    icons = ['icons/zoom-in.png',
             'icons/zoom-out.png',
             'icons/zoom-fit-best.png',
             None, None, None]         
    callbacks = [lambda : graph.setZoom(1.2),
                 lambda : graph.setZoom(0.8333333),
                 graph.resetView,
                 lambda : graph.setPlotType('line'),
                 lambda : graph.setPlotType('scatter'),
                 lambda : graph.setPlotType('histogram')]
    self.buttons = UtilityButtonList(self, texts=texts, icons=icons, doClone=False,
                                     doHelp=False, closeCmd=self.table.hideGraph,
                                     callbacks=callbacks, grid=(1,0))
   
    self.updateGraph(xCol, yCol)
    
  def updateGraph(self, xCol, yCol):
    
    self.xCol = xCol
    self.yCol = yCol
    
    columns = self.table.columns
    objects = self.table.objects
    
    if xCol is None:
      xName = 'Row'
    else:  
      xName = columns[xCol].heading
    
    if yCol is None:
      yName = 'Row'
    else:  
      yName = columns[yCol].heading
    
    self.graph.xAxis.name = xName
    self.graph.yAxis.name = yName
    
    values = [None] * len(objects)
    
    for i, obj in enumerate(self.table.objects):
      if xCol is None:
        x = i
      else:  
        x = columns[xCol].getValue(obj)
      
      if yCol is None:
        y = i
      else:  
        y = columns[yCol].getValue(obj)
      
      values[i] = (x,y)
    
    dataSet = GraphDataSet(values, None, '#008000')
    
    title = '%s vs %s' % (xName, yName)
    self.graph.updateData([dataSet], title)

TAB_FORMAT = 'Tab-separated'
COMMA_FORMAT = 'Comma-separated'
EXPORT_FORMATS = (TAB_FORMAT, COMMA_FORMAT)

class ObjectTableExport(BasePopup):

  def __init__(self, table):
  
    BasePopup.__init__(self, title='Export Table Text', transient=True, modal=True) 
    self.setWindowFlags(QtCore.Qt.Tool)

    self.table = table
    
    label = Label(self, 'Columns to export:', grid=(0,0), gridSpan=(1,2))
    
    labels = ['Row Number',] + [c.heading.replace('\n', ' ') for c in table.columns]
    values = [True] * len(labels)
    
    self.multi = MultiWidget(self, CheckButton, labels=labels, values=values,
                             minRows=None, maxRows=None, useImages=False,
                             maxRow=15, grid=(1,0), gridSpan=(1,2))
    self.multi.setFrameStyle(self.multi.Plain)
    
    label = Label(self, 'Export format:', grid=(3,0))
    
    self.formatPulldown = PulldownList(self, EXPORT_FORMATS, grid=(3,1))
    
    texts = ['Save File',]
    callbacks = [self.saveFile,]
    icons = ['icons/save.png',]
    buttons = UtilityButtonList(self, texts=texts, callbacks=callbacks,
                                icons=icons, doHelp=False,
                                doClone=False, closeCmd=self.close,
                                grid=(4,0), gridSpan=(1,2))
   
    self.setMaximumWidth(300)
    
  
  def saveFile(self):
    
    from gui.qtgui.FileSelect import selectSaveFile, FileType
    from gui.qtgui.MessageDialog import showWarning, showError

    values = self.multi.get()
    headings = self.multi.labels
    cols = [i for i, value in enumerate(values) if value]
    headings = [headings[i] for i in cols]
    #headings = [headings[i].encode('utf-8') for i in cols]
    
    if not headings:
      msg = 'No table columns selected for export.'
      showWarning('Failure', msg, parent=self)
      return
    
    if self.formatPulldown.get() == COMMA_FORMAT:
      fileTypes = [FileType('All files', ['*',]), FileType('CSV', ['*.csv',]),]
    else:
      fileTypes = [FileType('All files', ['*',]), FileType('Text', ['*.txt',]),]
    
    fileName = selectSaveFile(self, 'Select Export File', fileTypes=fileTypes)

    if not fileName:
      return
    
    if self.formatPulldown.get() == COMMA_FORMAT:
      sep = ','
    else:
      sep = '\t'

    try:
      fileObj = io.open(fileName, 'wb')
    except IOError as e:
      showError('File error', str(e), parent=self)
      return 

    fileObj.write(sep.join(headings) + '\n')
    
    row = 1
    indices = [i-1 for i in cols[1:]]
    objCols = self.table.columns
    for obj in self.table.objects:
      if values[0]:
        texts = ['%d' % row,]
      else:
        texts = []
      
      for i in indices:
        value = objCols[i].getValue(obj)
        
        if isinstance(value, float):
          if value == 0:
            text = '0'
          elif (abs(value) > 10000) or (abs(value) < 0.001):
            text = '%6.5e' % (value,)
          else:
            text = '%6.5f' % (value,)
          
        else:
          text = '%s' % (value,)
          
        texts.append(text)
      
      fileObj.write(sep.join(texts) + '\n')
      row += 1
  
    self.close()
    
if __name__ == '__main__':

  import sys

  from gui.qtgui.Application import Application
  from gui.qtgui.SpinBox import IntSpinBox
  from gui.qtgui.Colors import ColorPulldown, DEFAULT_COLORS
  from gui.qtgui.Base import Icon
  
  app = Application(sys.argv)
  popup = BasePopup(title='Test Table')
  
  class TestObj:
  
    def __init__(self, i, name, greek):
      self.i = i
      self.name = name
      self.greek = greek
      self.boolList = [False, True, False, True]
      self.isTicked = i % 2 == 1
      self.color = DEFAULT_COLORS[i % len(DEFAULT_COLORS)]
  
  objects = (TestObj(1, 'One', u'\u03B1'),
             TestObj(2, 'Two', u'\u03B2'),
             TestObj(3, 'Three',u'\u03B3'),
             TestObj(4, 'Four', u'\u03B4'),
             TestObj(5, 'Five', u'\u03B5'),
             TestObj(6, 'Six', u'\u03B6'),
             TestObj(7, 'Seven', u'\u03B7'),
             TestObj(8, 'Eight', u'\u03B8'),
             TestObj(9, 'Nine', u'\u03B9'),
             TestObj(10,'Ten', u'\u03BA'))
  
  def setBoolsEditValue(obj, value):
    obj.boolList = value

  def setBoolEditValue(obj, value):
    obj.isTicked = value
 
  def setEditValue(obj, value):
    print("Value Set", obj, value)
  
  def getName(obj):
    return obj.name
    
  def getGreek(obj):
    return obj.greek

  def getBgColor(obj):
    if obj.i % 2 == 0:
      return QtGui.QColor('#D0FFD0')

    else:
      return QtGui.QColor('#D0D0FF')
  
  def getBool(obj):
    return obj.isTicked
    
  def getEditValue(obj): # required to setup pulldown
    objects = [x*x for x in range(obj.i-5, obj.i+5)]
    texts = [str(x) for x in objects]
    index = objects.index(obj.i*obj.i)
    
    return texts, objects, index

  def getBoolListStr(obj):
    return ', '.join([x and 'Yes' or 'No' for x in obj.boolList]) or 'None'

  def getEditBoolList(obj): # required to setup mult-widget
    values = obj.boolList
    labels = ['label %d' % (x+1) for x in range(len(values))]
  
    return values, labels
  
  def getIcon(obj):
    
    color = obj.color
    
    return Icon(color=color)

  def setColor(obj, colorObj):
    
    if colorObj:
      r, g, b, a = colorObj.getRgb()
      obj.color = '#%02x%02x%02x' % (r,g,b) 
  
  def callback(obj, row, col):
    
    if not obj:
      print("Deselected")
    
    else:
      print("Callback", obj.name, row, col)
  
  def compareTextLen(objA, objB):
    
    return len(objA.name) < len(objB.name)
    
  def compareBool(objA, objB):
    
    return getBool(objA) and not getBool(objB)
               
  objectCols = [Column('x',   'i', setEditValue='i',
                       tipText='Tip C', alignment='C'),
                Column('Name', getName, setEditValue=setEditValue,
                       tipText='Tip A', orderFunc=compareTextLen),
                Column('Greek', 'greek', tipText='Tip B',
                       getColor=getBgColor, alignment='center'),
                Column('?', getBool, setEditValue=setBoolEditValue,
                       orderFunc=compareBool, tipText='Something to tick'),
                Column('x^2', lambda x:x.i*x.i, getEditValue=getEditValue,
                       setEditValue=setEditValue, tipText='Some int, editable',
                       alignment='r'),
                       #editClass=IntSpinBox, editArgs=(0, 0, 11))
                Column('1/x', lambda x:1.0/x.i, setEditValue=setEditValue,
                       tipText='Some float', format='%.4f',
                       editDecimals=4, editStep=0.1),
                Column('Color', 'color', setEditValue=setColor,
                       getIcon=getIcon, editClass=ColorPulldown),
                Column('Bools', getBoolListStr,
                       getEditValue=getEditBoolList, setEditValue=setBoolsEditValue,
                       editClass=MultiWidget, editArgs=(CheckButton,),
                       tipText='Some bools'),
                ]

  objTable = ObjectTable(popup, objectCols, objects,
                         callback=callback, multiSelect=True) 
  objTable.setObjects(objects[5:])
  
  def setColsA():
    objTable.setColumns(objectCols)

  def setColsB():
    objTable.setColumns(objectCols[:4])
  
  newObjs = (TestObj(1, 'Partridge in a pear tree', u'\u03B1'),
             TestObj(2, 'Turtle doves', u'\u03B2'),
             TestObj(3, 'French hens',u'\u03B3'),
             TestObj(4, 'Calling birds', u'\u03B4'),
             TestObj(5, 'Golden rigs', u'\u03B5'),
             TestObj(6, 'Geese a laying', u'\u03B6'))
  
  def setNewObjs():
    objTable.setObjects(objects)
    #objTable.setObjects(newObjs)
  
  def setSelection():
    objTable.setCurrentObjects(objects[7:10])
    #objTable.replaceCurrentObject(TestObj(100, 'Hundred', u'\u03B1'))
  
  def getSelection():
    for obj in objTable.getSelectedObjects():
      print(obj.i, obj.name)
    
  ButtonList(popup, ['Columns A', 'Columns B', 'New Objs', 'Select', 'Get Selected'],
                    [setColsA, setColsB, setNewObjs, setSelection, getSelection])
  
  app.start()


