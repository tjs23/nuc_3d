from PySide import QtCore, QtGui

from gui.qtgui.Base import Base, Icon
from gui.qtgui.Entry import FloatEntry
from gui.qtgui.Button import Button
from gui.qtgui.Label import Label  

QPainter = QtGui.QPainter
Antialiasing = QPainter.Antialiasing

PEN = QtGui.QColor('#000000')
BRUSH = QtGui.QColor(255, 255, 255, 128)

class Scrollbar(QtGui.QScrollBar, Base):

  def __init__(self, parent, startVal=0, endVal=100, loValue=None, hiValue=None,
               isVertical=False, callback=None, tracking=True,
               isInverted=False, **kw):
    
    QtGui.QScrollBar.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.setInvertedAppearance(isInverted)
    
    if startVal or endVal:
      if loValue is None:
        loValue = startVal
 
      if hiValue is None:
        hiValue = loValue + 1
 
      self.setRange(startVal, endVal)
      self.set(loValue, hiValue)

    self.isVertical = isVertical
    self.fontMetric = QtGui.QFontMetricsF(self.font())
    
    if isVertical:
      self.setOrientation(QtCore.Qt.Vertical)
      self.setSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
    else:
      self.setOrientation(QtCore.Qt.Horizontal)
      self.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)

    # Callback continuously (True)
    # Or only at intervals (False)
    self.setTracking(tracking)
    
    self.callback = callback
    
    self.valueChanged.connect(self._callback)
  
    #if showNumber and not tracking:
    #  self.sliderMoved.connect(self._sliderMoved)
    
    #if showNumber:
    #  self.sliderReleased.connect(self.update)
  
  def _callback(self):
  
    if self.callback and self.isVisible():
      lo, hi = self.get()
      self.callback(lo, hi)
      
  def _sliderMoved(self, pos):
    
    self.update()
    
  def setRange(self, startVal, endVal):
    
    startVal = int(startVal)
    endVal = int(endVal)
    
    assert startVal != endVal
    
    if startVal > endVal:
      startVal, endVal = endVal, startVal
    
    self.range = (startVal, endVal)    
    QtGui.QScrollBar.setRange(self, startVal, endVal)
  
  def set(self, lo, hi):
    
    if lo > hi:
      lo, hi = hi, lo
    
    startVal, endVal = self.range
    
    if hi > endVal:
      hi = endVal
 
    if lo < startVal:
      lo = startVal
     
    delta = hi-lo
    
    self.setValue(lo)
    self.setMaximum(endVal-delta)
    self.setPageStep(delta)
  
  def get(self):
  
    val = self.value()    
    return val, val+self.pageStep()

  def disable(self):

    self.setDisabled(True)

  def enable(self):

    self.setEnabled(True)

  def setState(self, state):

    self.setEnabled(state)

class FloatScrollbar(Scrollbar):

  def __init__(self, parent, startVal=-1.0, endVal=1.0,
               loValue=0.0, hiValue=0.1, decimals=2, **kw):
    
    Scrollbar.__init__(self, parent, None,
                       None, None, None, **kw)

    assert 0 <= decimals < 5
    
    self._fac = float(10**decimals)
    
    self.setRange(startVal, endVal)
    self.set(loValue, hiValue)
    
  def setRange(self, startVal, endVal):
    
    f = self._fac
    Scrollbar.setRange(self, startVal*f, endVal*f)
  
  def set(self, lo, hi):
    
    f = self._fac
    Scrollbar.set(self, lo*f, int(hi*f)+1)
    
  def get(self):

    f = self._fac
    lo, hi = Scrollbar.get(self)

    return lo/f, hi/f
    

class RegionScrollbar(QtGui.QWidget, Base):

  def __init__(self, parent, startVal=-1.0, endVal=1.0,
               loValue=0.0, hiValue=0.1, decimals=2,
               isVertical=False, callback=None,
               tracking=True, setMid=True, **kw):

    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.minSize = 10.0**-decimals
      
    self.expandButton = Button(self, '', self.expand, Icon('icons/range-expand.png'),
                               tipText='Expand range', grid=(0,0),
                               gridSpan=(2,1), stretch=(0,0))
    self.shrinkButton = Button(self, '', self.contract, Icon('icons/range-contract.png'),
                               tipText='Shrink range', grid=(0,1),
                               gridSpan=(2,1), stretch=(0,0))
                               
    self.shrinkButton.setSizePolicy(QtGui.QSizePolicy.Fixed,
                                    QtGui.QSizePolicy.Minimum)
                                                                      
    self.expandButton.setSizePolicy(QtGui.QSizePolicy.Fixed,
                                    QtGui.QSizePolicy.Minimum) 
                                                                       
    self.format = format = '%%.%df' % decimals
    
    self.startLabel = Label(self, format % startVal, grid=(1,2), hAlign=self.left)
    
    self.scrollbar = FloatScrollbar(self, startVal, endVal,
                                    loValue, hiValue, decimals,
                                    isVertical=isVertical, callback=callback,
                                    tracking=tracking, gridSpan=(1,3),
                                    grid=(0,2))

    self.endLabel = Label(self, format % endVal, grid=(1,4), hAlign=self.right)
    
    if setMid:
      self.entry = FloatEntry(self, 0.5*(loValue+hiValue), self.setMid, startVal,
                              endVal, decimals, grid=(1,3), stretch=(0,0))
      self.entry.setSizePolicy(QtGui.QSizePolicy.Fixed,
                               QtGui.QSizePolicy.Fixed)                                    
      self.scrollbar.sliderMoved.connect(self._sliderMoved)
  
  def _sliderMoved(self, val):
    
    a, b = self.scrollbar.get()
    mid = (a+b)/2.0
    self.entry.set(mid)
      
  def setMid(self, mid=None):
    
    if mid is None:
      mid = self.entry.get()
    
    a, b = self.scrollbar.get()
    delta = (b-a)/2.0
    minVal, maxVal = self.scrollbar.range
    
    if (mid-delta) < minVal:
      mid = minVal+delta
    
    elif (mid+delta) > maxVal:
      mid = maxVal-delta  
    
    self.scrollbar.set(mid-delta, mid+delta)
  
  def expand(self):
    
    a, b = self.scrollbar.get()
    delta = (b-a)/2.0
    mid = (a+b)/2.0
    
    delta = max(self.minSize, delta)
    delta *= 1.3
    self.scrollbar.set(mid-delta, mid+delta)
    
  def contract(self):

    a, b = self.scrollbar.get()
    delta = (b-a)/2.0
    mid = (a+b)/2.0
    
    delta /= 1.3
    delta = max(self.minSize, delta)
    self.scrollbar.set(mid-delta, mid+delta)
    
  def setRange(self, startVal, endVal):
   
    self.startLabel.set(self.format % startVal)
    self.endLabel.set(self.format % endVal)
    
    self.scrollbar.setRange(startVal, endVal)
  
  def set(self, lo, hi):
    
    self.scrollbar.set(lo, hi)
    
  def get(self):

    lo, hi = self.scrollbar.get(self)

    return lo, hi


    
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup

  app = Application()

  window = BasePopup()
  window.setSize(500,500)
  
  def endCallback(hi, lo):
    print("End Callback", hi, lo)
    
  def trackCallback(hi, lo):
    print("Track Callback", hi, lo) 
  
  s1 = Scrollbar(window, 10, 1, tracking=False,
                 callback=endCallback)
  s2 = Scrollbar(window, 1, 10, 3, 9,
                 callback=trackCallback)
  s3 = FloatScrollbar(window, -100.0, 100.0, 0.9, 17.3, isVertical=True,
                 callback=trackCallback)

  
  RegionScrollbar(window, 0.0, 20.0, 1.4, 5.7, callback=trackCallback)
  
  window.show()
  
  app.start()

