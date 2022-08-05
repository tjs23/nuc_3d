from PySide import QtCore, QtGui

from gui.qtgui.Base import Base, Icon
from gui.qtgui.Frame import Frame
from gui.qtgui.Label import Label
from gui.qtgui.Button import Button

class MultiWidget(QtGui.QFrame, Base):

  def __init__(self, parent, widgetClass, labels=None, values=None,
               callback=None, minRows=3, maxRows=None, useImages=True,
               maxRow=15, widgetArgs=None, **kw):

    QtGui.QFrame.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    self.setAutoFillBackground(True)
    self.setLineWidth(1)
    self.setFrameStyle(QtGui.QFrame.StyledPanel)
 
    #self.setStyleSheet("background-color: rgb(%d, %d, %d);" %  rgb)

    if not values:
      values = []
      
    if not labels:
      labels = [] 

    # Widget types currently: Entries, Pulldown, Checkbutton

    if useImages:
      self.crossImage = Icon('icons/dialog-cancel.png')
      self.tickImage  = Icon('icons/dialog-ok-apply.png')
      self.plusImage  = Icon('icons/list-add.png')
      self.minusImage = Icon('icons/list-remove.png')
    else:
      self.crossImage = None
      self.tickImage  = None
      self.plusImage  = None
      self.minusImage = None
    
    module = widgetClass.__module__
    self.widgetType  = module.split('.')[-1]
    self.widgetClass = widgetClass
    self.widgetArgs  = widgetArgs or {}
    self.callback    = callback
    self.maxRow = maxRow
    
    self.numRows = 0
    self.values = values
    self.labels = labels
    self.minRows = minRows
    self.maxRows = maxRows
    self.widgets = []
    self.labelWidgets = []

    self.topFrame = Frame(self, grid=(0,0), stretch=(1,1))
    self.bottomFrame = Frame(self, grid=(1,0), stretch=(0,0))

    self.addButton    = None
    self.removeButton = None
    self.commitButton = None
    self.cancelButton = None
    self.allButton    = None
    self.noneButton   = None

    self.set(self.values, self.labels)
  
    self.updateButtons()
  
  def keyPressEvent(self, event):
    
    QtGui.QWidget.keyPressEvent(event)
    
    key = event.key()
    
    if key == QtCore.Qt.Key_Escape:
      self.cancel()
      
    elif key == QtCore.Qt.Key_Return:
      self.commit()
    
    elif key == QtCore.Qt.Key_Escape:
      self.commit()

    
  def get(self):
  
    values = []
    
    if self.widgetType == 'PulldownList':
      for i in range(self.numRows):
        value = self.widgets[i].currentObject()
        if value is not None:
          values.append( value )
          
    else:
      for i in range(self.numRows):
        value = self.widgets[i].get()
        if (value is not None) and (value != ''):
          values.append( value )
    
    return values

  def set(self, values=None, labels=None):

    if values is not None:
      self.values = values
      
    if labels is not None:
      self.labels = labels      
        
    N = max( len(self.values), self.minRows)
    if self.maxRows is not None:
      N = min(N, self.maxRows)

    if self.numRows < N:
      for i in range(self.numRows, N):
        self.addRow()

    elif self.numRows > N:
      for i in range(N, self.numRows):
        self.removeRow(doClear=True)

    for i in range(self.numRows):
      if i < len(self.values):
        value = self.values[i]
      else:
        value = None
        
      if i < len(self.labels):
        option = self.labels[i]
      else:
        option = None

      self.updateWidget(i, value, option)

      
  def setlabels(self, labels):
  
    self.set(self.values, labels)


  def setValues(self, values):
  
    self.set(values, self.labels)


  def updateWidget(self, row, value, option):
 
    widget = self.widgets[row]
    if 'Entry' in self.widgetType:
      widget.set(value)
      label = self.labelWidgets[row]
      label.set(option or '')

    elif 'SpinBox' in self.widgetType:
      widget.set(value or 0)
      label = self.labelWidgets[row]
      label.set(option or '')
 
    elif self.widgetType == 'PulldownList':
      index = 0
      if value in self.labels:
        index = self.labels.index(value or self.labels[0])
      widget.setData(self.labels, self.labels, index)

    elif self.widgetType == 'CheckButton':
      widget.set(value)
      label = self.labelWidgets[row]
      label.set(option or '')

    else:
      raise 'Widget type %s not supported in MultiWidget' % self.widgetType
    

  def updateButtons(self):

    row = 0
    col = 0
    if self.widgetType != 'CheckButton':
      if (self.maxRows is None) or (self.numRows < self.maxRows):
        if not self.addButton:
          if self.plusImage:
            self.addButton = Button(self.bottomFrame, icon=self.plusImage,
                                   tipText='Add extra row', grid=(row, col),
                                   callback=self.addRow)
          else:
            self.addButton = Button(self.bottomFrame, text='+',
                                    callback=self.addRow, 
                                    tipText='Add extra row', grid=(row, col))
            
        col += 1
        self.addButton.show()
        
      elif self.addButton:
        self.addButton.hide()

      if self.numRows > self.minRows:
        if not self.removeButton:    
          if self.minusImage:
            self.removeButton = Button(self.bottomFrame, icon=self.minusImage,
                                       tipText='Remove last row', 
                                       callback=self.removeRow, grid=(row, col))
          else:
            self.removeButton = Button(self.bottomFrame, text='-',
                                       callback=self.removeRow, 
                                       tipText='Remove last row', grid=(row, col))
        col += 1
        self.removeButton.show()
        
      elif self.removeButton:
        self.removeButton.hide()
    
    elif self.numRows > 5:
      
      if not self.allButton:
        self.allButton = Button(self.bottomFrame, text='All',
                                callback=self._selectAll,
                                tipText='Select all items', grid=(row, col))
      col += 1
      
      if not self.noneButton:
        self.noneButton = Button(self.bottomFrame, text='None',
                                 callback=self._selectNone,
                                 tipText='Select no items', grid=(row, col))
      col += 1

      
    if self.callback:
      if not self.commitButton:
        if self.tickImage:
          self.commitButton = Button(self.bottomFrame, icon=self.tickImage,
                                     tipText='Commit selection',
                                     callback=self.commit, grid=(row, col))
        else:
          self.commitButton = Button(self.bottomFrame, text='OK',
                                     callback=self.commit, 
                                     tipText='Commit selection', grid=(row, col))

      col += 1
 
      if not self.cancelButton:
        if self.crossImage:
          self.cancelButton = Button(self.bottomFrame, icon=self.crossImage,
                                     callback=self.cancel,  
                                     tipText='Abort selection', grid=(row, col))
        else:
          self.cancelButton = Button(self.bottomFrame, text='Cancel',
                                     callback=self.cancel, 
                                     tipText='Abort selection', grid=(row, col))

        
  def _selectAll(self):
  
    values = [True] * len(self.values)
    self.set(values, self.labels)

  def _selectNone(self):
  
    values = [False] * len(self.values)
    self.set(values, self.labels)


  def addWidget(self, i, value=None):
    
    row = i % self.maxRow
    col = i // self.maxRow
 
    if 'Entry' in self.widgetType:
      option = ''
      if i < len(self.labels):
        option = self.labels[i]

      label = Label(self.topFrame, text=option or '',
                    grid=(row, 2*col), stretch=(1,0))
      self.labelWidgets.append(label)
      
      widget = self.widgetClass(self.topFrame, text=value,
                                grid=(row, 2*col+1), stretch=(1,1),
                                **self.widgetArgs)

    elif 'SpinBox' in self.widgetType:
      option = ''
      if i < len(self.labels):
        option = self.labels[i]

      label = Label(self.topFrame, text=option or '',
                    grid=(row, 2*col), stretch=(1,0))
      self.labelWidgets.append(label)
      
      widget = self.widgetClass(self.topFrame, value=value or 0,
                                grid=(row, 2*col+1), stretch=(1,1),
                                **self.widgetArgs)
 
    elif self.widgetType == 'PulldownList':

      index = 0
      if value in self.labels:
        index = self.labels.index(value or self.labels[0])
      
      widget = self.widgetClass(self.topFrame, grid=(row, 2*col),
                                gridSpan=(1,2), stretch=(1,1) )
      widget.setData(self.labels, self.labels, index)
      

    elif self.widgetType == 'CheckButton':

      widget = self.widgetClass(self.topFrame, grid=(row, 2*col), stretch=(1,0))
      widget.set(value)
      
      option = ''
      if row < len(self.labels):
        option = self.labels[i]

      label = Label(self.topFrame, text=option or '',
                    grid=(row, 2*col+1), stretch=(1,1))
      self.labelWidgets.append(label)

    else:
      msg = 'Widget type %s not supported in MultiWidget' % self.widgetType
      raise Exception(msg)

    return widget


  def commit(self, event=None):

    values = self.get()
          
    if self.callback:
      self.callback(values)


  def cancel(self, event=None):
  
    if self.callback:
      self.callback(None)


  def addRow(self):

    if (self.maxRows is None) or (self.numRows < self.maxRows):
      i = self.numRows
      self.numRows += 1
      row = i % self.maxRow
      col = i // self.maxRow
      col2 = 2*col
      
      if self.numRows > len(self.widgets):
        widget = self.addWidget(i)
        self.widgets.append(widget)

      else:
        widget = self.widgets[i]
        if 'Entry' in self.widgetType:
          label = self.labelWidgets[i]
          label.show()
          widget.show()

        elif 'SpinBox' in self.widgetType:
          label = self.labelWidgets[i]
          label.show()
          widget.show()
            
        elif self.widgetType == 'CheckButton':
          widget.grid(row, col2)
          label = self.labelWidgets[i]
          label.show()
        
        else:
          widget.show()

      if row == 0:
        self.topFrame.show()


      self.updateButtons()

  def removeRow(self, doClear=False):
 
    if self.numRows > self.minRows:
      self.numRows -= 1
      row    = self.numRows
      widget = self.widgets[row]
      widget.hide()
      
      if self.widgetType == 'CheckButton':
        label = self.labelWidgets[row]
        label.hide()
        
      elif 'Entry' in self.widgetType:
        label = self.labelWidgets[row]
        label.hide()

      elif 'SpinBox' in self.widgetType:
        label = self.labelWidgets[row]
        label.hide()

      if doClear:
        self.updateWidget(row, None, '') 

      if self.numRows == 0:
        self.topFrame.hide()

      self.updateButtons()

if __name__ == '__main__':

  from gui.qtgui.Entry import Entry, FloatEntry
  from gui.qtgui.PulldownList import PulldownList
  from gui.qtgui.CheckButton import CheckButton
  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.Base import Align
  
  app = Application()
  popup = BasePopup(title='Test Frame')
  popup.setSize(350, 600)
  def fn(values):
    print(values)

  mw = MultiWidget(popup, CheckButton, callback=fn, minRows=0, hAlign=Align.left,
                   labels=['One','Two','Three'], values=[1,0,1],
                   tipText='Some advice')
  
  mw.set(labels=['Easy','As','One','Two','Three'],
         values=[1,0,1,0,1])

  mw2 = MultiWidget(popup, FloatEntry, callback=fn, hAlign=Align.left,
                    minRows=3, maxRows=3, useImages=False,
                    labels=['1H:','13C:','15N:'],
                    values=[1.0,0.1,0.15])
  
  mw2 = MultiWidget(popup, Entry, callback=fn,
                    minRows=0, useImages=False, 
                    values=['Easy','As','One','Two','Three'])

  mw3 = MultiWidget(popup, PulldownList, callback=None, hAlign=Align.left,
                    values=['Bjorn','Benny','Benny','Bjorn'],
                    labels=['Bjorn','Agnetha','Benny','Anni-Frid'])

  
  app.start()
