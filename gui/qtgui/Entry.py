from PySide import QtCore, QtGui

from gui.qtgui.Base import Base

from gui.qtgui.Frame import Frame
from gui.qtgui.Label import Label

# Width?
# Allow setting of max lengh based on data model?

import re

SPLIT_REG_EXP = re.compile(',?\s*')
SEPARATOR = ', '
MAXINT = 2**31-1
INFINITY = float('Inf')


class Entry(QtGui.QLineEdit, Base):

  def __init__(self, parent, text='', callback=None, maxLength=32, 
               listener=None, stripEndWhitespace=True, **kw):
    
    QtGui.QLineEdit.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.setText(self.valueToText(text))
    self.setMaxLength(maxLength)
    
    self._isAltered = False
    self._stripEndWhitespace = stripEndWhitespace
    self.callback = callback
    
    self.textChanged.connect(self._changed)
    
    self.connect(self, QtCore.SIGNAL('returnPressed()'), self._callback)
    self.connect(self, QtCore.SIGNAL('editingFinished()'), self._callback)
      
    if listener:
      if isinstance(listener, (set, list, tuple)):
        for signal in listener:
          signal.connect(self.set)
        
      else:
        listener.connect(self.set)
  
  def _callback(self):
    
    if self.callback and self._isAltered:
      self.callback(self.get())
      self._isAltered = False
  
  def _changed(self):
    
    self._isAltered = True
  
  def textToValue(self, text):
    # Overwritten in subclasses to make float, int etc
    
    if self._stripEndWhitespace:
      text = text.strip()
    
    return text or None

  def valueToText(self, value):
    # Overwritten in subclasses to convert float, int
    
    return value or ''
  
  def get(self):
    
    return self.textToValue(self.text())
    
  def set(self, value, doCallback=True):
    
    self.setText(self.valueToText(value))
    
    if doCallback:
      self._callback()

class IntEntry(Entry):

  def __init__(self, parent, text=0, callback=None,
               minValue=-MAXINT, maxValue=MAXINT,  **kw):  
    
    Entry.__init__(self, parent, text, callback, **kw)
   
    valid = QtGui.QIntValidator(minValue, maxValue, self)
    self.setValidator(valid)

  def textToValue(self, text):
    
    if not text:
      return None
    else:
      return int(text)

  def valueToText(self, value):
    
    if value is None:
      return ''
    else:
      return '%s' % value

  def setRange(self, minValue=-MAXINT, maxValue=MAXINT):
  
    valid = QtGui.QIntValidator(minValue, maxValue, self)
    self.setValidator(valid)
    

class FloatEntry(Entry):
  
  decimals = 4
  
  def __init__(self, parent, text=0.0, callback=None,
               minValue=-INFINITY, maxValue=INFINITY,
               decimals=4, **kw):  

    Entry.__init__(self, parent, text, callback, **kw)
    
    self.decimals = decimals
    self.setText(self.valueToText(text))
    
    valid = QtGui.QDoubleValidator(minValue, maxValue, decimals, self)
    self.setValidator(valid)

  def textToValue(self, text):
    
    if not text:
      return None
    else:
      return float(text)

  def valueToText(self, value):
    
    if value is None:
      text = ''
    elif value == 0:
      text = '0.0'
    elif abs(value) > 999999 or abs(value) < 0.00001:
      format = '%%.%de' % (self.decimals)
      text = format % value
    else:
      format = '%%.%df' % (self.decimals)
      text = format % value   
       
    return text

  def setRange(self, minValue=-INFINITY, maxValue=INFINITY):
  
    valid = QtGui.QIntValidator(minValue, maxValue, self)
    self.setValidator(valid)

class RegExpEntry(Entry):

  def __init__(self, parent, text='', callback=None, **kw):  
    
    Entry.__init__(self, parent, text, callback, **kw)
    
    self.setValidator(QtGui.QRegExpValidator)


class ArrayEntry(Entry):

  def __init__(self, parent, text='', callback=None, **kw):  
    
    Entry.__init__(self, parent, text, callback, **kw)
  
  def textToValue(self, text):
  
    return re.split(SPLIT_REG_EXP, text) or []
  
  def valueToText(self, array):
  
    return SEPARATOR.join(array)


class IntArrayEntry(IntEntry):

  def __init__(self, parent, text='', callback=None, **kw):  
    
    IntEntry.__init__(self, parent, text, callback, **kw)
  
  def textToValue(self, text):
  
    array = re.split(SPLIT_REG_EXP, text) or []
    return  [IntEntry.textToValue(self, x) for x in array]
  
  def valueToText(self, values):
    
    texts = [IntEntry.valueToText(self, x) for x in values]
    return SEPARATOR.join(texts)


class IntRangesEntry(IntEntry):

  def __init__(self, parent, text='', callback=None,
               minValue=0, maxValue=MAXINT, **kw):  
    
    Entry.__init__(self, parent, text, callback, **kw)
    self.maxValue = maxValue
    self.minValue = minValue
  
  def textToValue(self, text):
    
    array = []
    
    items = re.split('\s*[/,;\s]\s*', text) or []
    
    for item in items:
      
      try:
        x = int(item)
        array.append(x)
        
      except ValueError:
        valRange = re.split('\s*[-:]|\.{2}\s*', item) or []
        
        if len(valRange) > 1:
          x = valRange[0]
          y = valRange[-1]
        
          try:
             x = int(x)
             y = int(y)
             x, y = sorted([x,y])
             array += list(range(x, y+1))
             
          except ValueError:
            continue
        
      except Exception:
        continue
    
    array = sorted(set(array))
    array = [x for x in array if x >= self.minValue and x <= self.maxValue]
     
    return  array

  def _callback(self):
    
    
    if self.callback and self._isAltered:
      value = self.get()
      self.setText(self.valueToText(value))
      self.callback(value)
      self._isAltered = False
  
  def valueToText(self, values):
    
    values = sorted(values)
    groups = []
    group = []
    
    for x in values:
      if group:
        prev = group[-1]
        
        if x == prev+1:
          group.append(x)
          
        else:
          groups.append(group)
          group = [x]
      
      else:
        group = [x]
            
    if group:
      groups.append(group)  
    
    texts = []
    for group in groups:
      if len(group) == 1:
        texts.append('%d' % (group[0]))
      else:
        texts.append('%d-%d' % (group[0],group[-1]))
    
    return SEPARATOR.join(texts)


class FloatArrayEntry(FloatEntry):

  def __init__(self, parent, text='', callback=None, **kw):  
    
    FloatEntry.__init__(self, parent, text, callback, **kw)
  
  def textToValue(self, text):
  
    array = re.split(SPLIT_REG_EXP, text) or []
    return  [FloatEntry.textToValue(self, x) for x in array]
  
  def valueToText(self, values):
    
    texts = [FloatEntry.valueToText(self, x) for x in values]
    return SEPARATOR.join(texts)


class LabeledEntry(Frame):

  def __init__(self, parent, labelText, entryText='', callback=None, maxLength=32,  tipText=None, **kw):  
    
    Frame.__init__(self, parent, tipText=tipText, **kw)
  
    self.label = Label(self, labelText, tipText=tipText, grid=(0,0))
    self.entry = Entry(self, entryText, callback, maxLength,
                       tipText=tipText, grid=(0,1))
    
  def getLabel(self):

    return self.label.get()

  def setLabel(self, text):

    self.label.set(text)

  def getEntry(self):

    return self.entry.get()

  def setEntry(self, text):

    self.entry.set(text)

class LabeledIntEntry(LabeledEntry):

  def __init__(self, parent, labelText, entryText='', callback=None,
               minValue=-MAXINT, maxValue=MAXINT,  tipText=None, **kw):  
    
    Frame.__init__(self, parent, tipText=tipText, **kw)
  
    self.label = Label(self, labelText, tipText=tipText, grid=(0,0))
    self.entry = IntEntry(self, entryText, callback, minValue,
                          maxValue, tipText=tipText, grid=(0,1))


class LabeledFloatEntry(LabeledEntry):

  def __init__(self, parent, labelText, entryText='', callback=None,
               minValue=-MAXINT, maxValue=MAXINT, decimals=4, tipText=None, **kw):  
    
    Frame.__init__(self, parent, tipText=tipText, **kw)
  
    self.label = Label(self, labelText, tipText=tipText, grid=(0,0))
    self.entry = FloatEntry(self, entryText, callback, minValue,
                            maxValue, decimals, tipText=tipText, grid=(0,1))


if __name__ == '__main__':

  from gui.qtgui.Application import Application

  app = Application()

  window = QtGui.QWidget()
  
  def callback(value):
    print("Callback", value)
  
  Entry(window, 'Start Text', callback)
  
  ArrayEntry(window, ['A','C','D','C'], callback)
  
  IntEntry(window, 123, callback)
  
  IntRangesEntry(window,  [4, 5, 6, 7], callback)
  
  FloatEntry(window, 2.818, callback)
  
  e = FloatArrayEntry(window, [1,2,4], callback, decimals=2)
  e.set([1e12, -0.7e-5, 9.75])
  
  LabeledEntry(window, 'Text:', 'Initial val', callback)
  
  LabeledIntEntry(window, 'Int:', 0, callback)
  
  LabeledFloatEntry(window, 'Float:', 0.7295, callback, decimals=8)
  
  window.show()
  
  app.start()

