from PySide import QtCore, QtGui

from gui.qtgui.Base import Base

QPainter = QtGui.QPainter
Antialiasing = QPainter.Antialiasing

PEN = QtGui.QColor('#000000')
BRUSH = QtGui.QColor(255, 255, 255, 128)


class Slider(QtGui.QSlider, Base):

  def __init__(self, parent, startVal=0, endVal=100, value=None,
               direction='h', step=1, bigStep=None, callback=None,
               tracking=True, showNumber=True, tickInterval=None,
               tickPosition=None, listener=None, **kw):
    
    QtGui.QSlider.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.callback = callback

    if value is None:
      value = startVal
    
    if not bigStep:
      bigStep = step
    
    if tickInterval:
      if not tickPosition:
        if 'h' in direction.lower():
          tickPosition = self.TicksBelow
        else:
          tickPosition = self.TicksRight
    
      self.setTickInterval(tickInterval)
      self.setTickPosition(tickPosition)
    
    self.showNumber = showNumber
    self.setRange(startVal, endVal)
    self.setStep(step, bigStep)
    self.set(value)
    self.fontMetric = QtGui.QFontMetricsF(self.font())
    
    if 'h' in direction.lower():
      self.setOrientation(QtCore.Qt.Horizontal)
    else:
      self.setOrientation(QtCore.Qt.Vertical)

    # Callback continuously (True)
    # Or only at intervals (False)
    self.setTracking(tracking)
    
  
    if showNumber and not tracking:
      self.connect(self, QtCore.SIGNAL('sliderMoved(int)'), self._redraw)
    
    if showNumber:
      self.connect(self, QtCore.SIGNAL('sliderReleased()'), self.update)
    
    self.connect(self, QtCore.SIGNAL('valueChanged(int)'), self._callback)
    
    if listener:
      if isinstance(listener, (set, list, tuple)):
        for signal in listener:
          signal.connect(self.setValue)
        
      else:
        listener.connect(self.setValue)
      
  def _redraw(self, pos):
    
    self.update()
  
  def paintEvent(self, event):
    
    pos = self.pos()
    x = pos.x()
    y = pos.y()
    
    rect = QtCore.QRect(x, y, 12, self.height())
    
    event2 = QtGui.QPaintEvent(rect)
    retVal = QtGui.QSlider.paintEvent(self, event2)
    
    if self.showNumber and self.isSliderDown():
      painter = QPainter()
      painter.begin(self)
      painter.setRenderHint(Antialiasing)
 
      text = str(self.get())
      bbox = self.fontMetric.tightBoundingRect(text)
      
      h = bbox.height()
      w = bbox.width()
      y = self.height() / 2.0 +  h / 2.0
      x = self.width() / 2.0 - w / 2.0
      point = QtCore.QPointF(x, y)
      
      painter.setPen(PEN)
      painter.drawText(x, y, text)

      painter.setPen(BRUSH)
      painter.setBrush(BRUSH)
      painter.drawRect(x-2, y-h-2, w+4, h+4)
 
      painter.end()
     
    return retVal
    
  def setRange(self, startVal, endVal):
    
    startVal = int(startVal)
    endVal = int(endVal)
    
    assert startVal != endVal
    
    if startVal > endVal:
      self.setInvertedAppearance(True)
      startVal, endVal = endVal, startVal
    else:
      self.setInvertedAppearance(False)

    value = self.get()
    
    if startVal <= value <= endVal:
      callback = self.callback
      self.callback = None
      QtGui.QSlider.setRange(self, startVal, endVal)
      self.callback = callback
    
    else:
      QtGui.QSlider.setRange(self, startVal, endVal)
  
  def setStep(self, step, bigStep=None):
  
    self.setSingleStep(step)
    
    if bigStep:
      self.setPageStep(bigStep)
  
  def set(self, value, doCallback=True):
    
    if not doCallback:
      callback = self.callback
      self.callback = None
      self.setValue(int(value))
      self.callback = callback
      
    else:
      self.setValue(int(value))
  
  def get(self):
  
    return self.value()  
  
  def _callback(self, callback):
  
    if self.callback:
      self.callback(self.value())

  def disable(self):

    self.setDisabled(True)

  def enable(self):

    self.setEnabled(True)

  def setState(self, state):

    self.setEnabled(state)


class FloatSlider(Slider):

  def __init__(self, parent, startVal=0.0, endVal=1.0, value=None,
               direction='h', step=0.1, bigStep=None, callback=None,
               tracking=True, showNumber=True, tickInterval=None,
               tickPosition=None, listener=None, decimals=2, **kw):

    QtGui.QSlider.__init__(self, parent)
    Base.__init__(self, parent, **kw)

    assert 0 <= decimals < 5
    
    self.callback = callback
    self._fac = float(10**decimals)
    
    if value is None:
      value = startVal
    
    if not bigStep:
      bigStep = step
    
    if tickInterval:
      if not tickPosition:
        if 'h' in direction.lower():
          tickPosition = self.TicksBelow
        else:
          tickPosition = self.TicksRight
    
      self.setTickInterval(tickInterval)
      self.setTickPosition(tickPosition)
    
    self.showNumber = showNumber
    self.setRange(startVal, endVal)
    self.setStep(step, bigStep)
    self.set(value)
    self.fontMetric = QtGui.QFontMetricsF(self.font())
    
    if 'h' in direction.lower():
      self.setOrientation(QtCore.Qt.Horizontal)
    else:
      self.setOrientation(QtCore.Qt.Vertical)

    self.setTracking(tracking)
    
    self.callback = callback
  
    if showNumber and not tracking:
      self.connect(self, QtCore.SIGNAL('sliderMoved(int)'), self._redraw)
    
    if showNumber:
      self.connect(self, QtCore.SIGNAL('sliderReleased()'), self.update)
    
    self.connect(self, QtCore.SIGNAL('valueChanged(int)'), self._callback)
    
    if listener:
      listener.connect(self.set)
  
  def _callback(self, callback):
  
    if self.callback:
      self.callback(self.value()/self._fac)

  def setTickInterval(self, tickInterval):
  
     f = self._fac
     Slider.setTickInterval(self, tickInterval*f)
  
  def setStep(self, step, bigStep=None):
  
    f = self._fac
    
    if bigStep:
      bigStep *= f
      
    Slider.setStep(self, step*f, bigStep)

  def setRange(self, startVal, endVal):
    
    f = self._fac
    Slider.setRange(self, startVal*f, endVal*f)
  
  def set(self, value, doCallback=True):
    
    f = self._fac
    Slider.set(self, int(value*f), doCallback)
    
  def get(self):

    f = self._fac
    value = Slider.get(self)

    return value/f

if __name__ == '__main__':

  from gui.qtgui.Application import Application

  app = Application()

  window = QtGui.QWidget()
  
  s1 = Slider(window, 10, 1, tracking=False)
  s2 = Slider(window, 1, 10, 7, showNumber=False)
  s3 = Slider(window, -100, 100, -10, direction='v',
              step=10, bigStep=50)
  
  def endCallback(v):
    print("End Callback", v)
    
  def trackCallback1(v):
    print("Track Callback 1", v)
    
  def trackCallback2(v):
    print("Track Callback 2", v)
  
  s1.callback = endCallback
  s2.callback = trackCallback1
  s3.callback = trackCallback2
  
  window.show()
  
  app.start()

