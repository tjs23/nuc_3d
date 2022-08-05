from PySide import QtCore, QtGui

from gui.qtgui.Base import Base

MAXINT = 2**31-1
INFINITY = float('Inf')

class _SpinBase(Base):
       
  def __init__(self, parent, value, minValue, maxValue,
               step=1, callback=None, suffix=None,
               prefix=None, listener=None, **kw):
    
    Base.__init__(self, parent, **kw)
    
    if suffix:
      self.setSuffix(suffix)
      
    if prefix:
      self.setPrefix(prefix)  
    
    self.callback = callback
    self.setRange(minValue, maxValue)   
    self.setStep(step)
    self.setValue(value)
    self.connect(self, QtCore.SIGNAL('editingFinished()'), self._callback)
      
    if listener:
      if isinstance(listener, (set, list, tuple)):
        for signal in listener:
          signal.connect(self.set)
        
      else:
        listener.connect(self.set)
  
  def _callback(self):

    if self.callback:
      self.callback(self.value())
  
  
  def _convertValue(self, value):
  
    return int(value)
    
    
  def get(self):
  
    return self.value()
  
  
  def set(self, value, doCallback=True):
    
    if not doCallback:
      callback = self.callback
      self.callback = None
      self.setValue(self._convertValue(value))
      self.callback = callback
      
    else:
      self.setValue(self._convertValue(value))
      
    return 
  
  def setStep(self, step):
  
    self.setSingleStep(self._convertValue(step))
  
  
  def setRange(self, minValue, maxValue):
    
    if minValue > maxValue:
      minValue, maxValue = maxValue, minValue
    
    value = self.get()
    
    minValue = self._convertValue(minValue)
    maxValue = self._convertValue(maxValue)
    
    if minValue <= value <= maxValue:
      callback = self.callback
      self.callback = None
      QtGui.QSpinBox.setRange(minValue, maxValue)
      self.callback = callback
    
    else:
      QtGui.QSpinBox.setRange(minValue, maxValue)
 


class IntSpinBox(QtGui.QSpinBox, _SpinBase):

  def __init__(self, parent, value=0.0, minValue=-MAXINT, maxValue=MAXINT, step=1,
               callback=None, suffix=None, prefix=None, multiplier=None, **kw):
    
    QtGui.QSpinBox.__init__(self, parent)
    _SpinBase.__init__(self, parent, value, minValue, maxValue, step,
                       callback, suffix, prefix, **kw)
  
    self.connect(self, QtCore.SIGNAL('valueChanged(int)'), self._callback)
    self.multiplier = multiplier

  def stepBy(self, steps):
    
    value = self.value()
    if self.multiplier and value:
    
      if steps > 0:
        value *= self.multiplier
      else:
        value /= self.multiplier
      
      self.setValue(value)
      
    else:
      QtGui.QSpinBox.stepBy(self, steps)
  
  def _convertValue(self, value):
  
    return int(value)


class FloatSpinBox(QtGui.QDoubleSpinBox, _SpinBase):

  def __init__(self, parent, value=0, minValue=-MAXINT, maxValue=MAXINT, step=1.0,
               callback=None, suffix=None, prefix=None, multiplier=None, decimals=None, **kw):
    
    QtGui.QDoubleSpinBox.__init__(self, parent)
    _SpinBase.__init__(self, parent, value, minValue, maxValue, step,
                       callback, suffix, prefix, **kw)

    self.connect(self, QtCore.SIGNAL('valueChanged(double)'), self._callback)
    self.multiplier = multiplier
    
    if decimals is not None:
      self.setDecimals(decimals)

  def stepBy(self, steps):
  
    value = self.value()
    if self.multiplier and value:
    
      if steps > 0:
        value *= self.multiplier
      else:
        value /= self.multiplier
      
      self.setValue(value)
      
    else:
      QtGui.QDoubleSpinBox.stepBy(self, steps)
  
  def _convertValue(self, value):
  
    return float(value)
   
  
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  
  app = Application()

  window = BasePopup()
  window.setSize(300, 120)
  
  def callback(*arg):
    print("Callback")

  IntSpinBox(window, 25, 0, 100, 10, callback)
  FloatSpinBox(window, 0.9, 0.1, 1.5, 0.1, callback)
  
  window.show()
  
  app.start()
