from PySide import QtGui, QtCore, QtSvg

from os import path

ICON_DIR = path.join('..','icons') 

class Align(object):

  top = 'top'
  left = 'left'
  right = 'right'
  center = 'center'
  bottom = 'bottom'
  
  _hAlignDict = {left:   QtCore.Qt.AlignLeft,
                center: QtCore.Qt.AlignHCenter,
                right:  QtCore.Qt.AlignRight}
 
  _vAlignDict = {top:    QtCore.Qt.AlignTop,
                center: QtCore.Qt.AlignVCenter,
                bottom: QtCore.Qt.AlignBottom}


class SizePolicy(object):

  fixed = 'fixed'      
  minimum = 'minimum'
  maximum = 'maximum'
  preferred = 'preferred'
  expanding = 'expanding'
  minimumExpanding = 'minimumExpanding'
  ignored = 'ignored'
                
  _policyDict = {fixed: QtGui.QSizePolicy.Fixed,
                minimum: QtGui.QSizePolicy.Minimum,
                maximum: QtGui.QSizePolicy.Maximum,
                preferred: QtGui.QSizePolicy.Preferred,
                expanding: QtGui.QSizePolicy.Expanding,
                minimumExpanding: QtGui.QSizePolicy.MinimumExpanding,
                ignored: QtGui.QSizePolicy.Ignored}              


class Base(Align, SizePolicy):
  
  def __init__(self, parent, tipText=None, grid=(None, None), gridSpan=(1, 1),
               stretch=(0,0), hAlign=None, vAlign=None, hPolicy=None,
               vPolicy=None, bgColor=None, isFloatWidget=False):
    
    from gui.qtgui.Layout import FlowLayout
    #Align.__init__(self)
    #SizePolicy.__init__(self)
    
    if tipText: # Tool tips
      self.setToolTip(tipText)
    
    if self.parent() and not isFloatWidget:
      self.getParent = self.parent
      # Setup gridding within parent
      layout = parent.layout()
 
      if grid and not layout:
        layout = QtGui.QGridLayout(parent)
        layout.setSpacing(2)
        layout.setContentsMargins(2,2,2,2)
        parent.setLayout( layout )
      
      elif not layout:
        layout = QtGui.QBoxLayout(QtGui.QBoxLayout.TopToBottom, parent)
        layout.setSpacing(2)
        layout.setContentsMargins(2,2,2,2)
        parent.setLayout( layout )        
      
      if isinstance(layout, QtGui.QGridLayout): 
        rowSpan, colSpan = gridSpan
        row, col = self._getRowCol(grid)
        align = self._getAlignment(hAlign, vAlign)
        
        rowStr, colStr = stretch
        layout.setRowStretch(row, rowStr)
        layout.setColumnStretch(col, colStr)
        layout.addWidget(self, row, col, rowSpan, colSpan, align)
       
      elif isinstance(layout, QtGui.QBoxLayout):
        if isinstance(stretch, (tuple, list)):
          wStretch = stretch[0]
        else:
          wStretch = stretch
        
        align = self._getAlignment(hAlign, vAlign)
        layout.addWidget(self, wStretch, align)

      elif isinstance(layout, FlowLayout):
        layout.addWidget(self)
        
    
    if hPolicy or vPolicy:
      policy = QtGui.QSizePolicy()
      
      if hPolicy:
        policy.setHorizontalPolicy(self._policyDict.get(hPolicy))
      
      if vPolicy:
        policy.setVerticalPolicy(self._policyDict.get(vPolicy))

      self.setSizePolicy(policy)               
 
    # Setup colour overrides (styles used primarily)
    if bgColor:
      self.setAutoFillBackground(True)
      rgb = QtGui.QColor(bgColor).getRgb()[:3]
      self.setStyleSheet("background-color: rgb(%d, %d, %d);" %  rgb)


  def _getAlignment(self, hAlign=None, vAlign=None):

    """
    Align options
      1 Qt::AlignLeft
      2 Qt::AlignRight
      4 Qt::AlignHCenter
      8 Qt::AlignJustify
     16 Qt::AlignAbsolute
     32 Qt::AlignTop
     64 Qt::AlignBottom
    128 Qt::AlignVCenter
    """
    
    alignment = self._hAlignDict.get(hAlign, 0) | self._vAlignDict.get(vAlign, 0)         
    
    return QtCore.Qt.Alignment(alignment)
  
  
  def _fixFileExtension(self, filePath, extensions):
    
    for i, ext in enumerate(extensions):
      ext = ext.lower()
      
      if ext[0] == '*':
        ext = ext[1:]
      if ext[0] != '.':
        ext = '.' + ext
        
      extensions[i] = ext  
    
    root, ext = path.splitext(filePath)
    ext = ext.lower()
    
    if (not ext) or (ext not in extensions):
      ext = extensions[0]
      filePath = root + ext
  
    return filePath
  
  def exportPdf(self, filePath=None):
    
    if not filePath:
      from gui.qtgui.FileSelect import selectSaveFile, FileType
      fileTypes = [FileType('PDF', ['*.pdf']),]
      filePath = selectSaveFile(self, 'Select PDF file name',
                                fileTypes=fileTypes)
 
    if not filePath:
      return
    
    filePath = self._fixFileExtension(filePath, ['.pdf',])
    device = QtGui.QPrinter()
    device.setOutputFileName(filePath)
    device.setOutputFormat(device.PdfFormat)

    painter = QtGui.QPainter()
    painter.begin(device)
    
    if isinstance(self, QtGui.QGraphicsView):
      self.render(painter)
    else:
      self.render(painter, QtCore.QPoint(0,0))
    
    painter.end() 
     
  def exportSvg(self, filePath=None):
    
    if not filePath:
      from gui.qtgui.FileSelect import selectSaveFile, FileType
      fileTypes = [FileType('SVG', ['*.svg']),]
      filePath = selectSaveFile(self, 'Select SVG file name',
                                fileTypes=fileTypes)
 
    if not filePath:
      return

    filePath = self._fixFileExtension(filePath, ['.svg',])
    device = QtSvg.QSvgGenerator()
    device.setFileName(filePath)
    
    if isinstance(self, QtGui.QGraphicsView):
      size = QtCore.QSize(210, 297) # A ratio
    else:
      size = self.size()
    
    device.setSize(size)

    painter = QtGui.QPainter()
    painter.begin(device)

    if isinstance(self, QtGui.QGraphicsView):
      self.render(painter)
    else:
      self.render(painter, QtCore.QPoint(0,0))
    
    painter.end() 
     
  def exportJpeg(self, filePath=None):
  
    exts = ['.jpg','.jpeg','.jpr']
    
    if not filePath:
      from gui.qtgui.FileSelect import selectSaveFile, FileType
      fileTypes = [FileType('JPEG', ),exts]
      filePath = selectSaveFile(self, caption='Select image file name',
                                directory=None, fileTypes=fileTypes)
 
                   
    if filePath:
      filePath = self._fixFileExtension(filePath, exts)
      widget = self
      pixmap = QtGui.QPixmap.grabWidget(widget, widget.rect())
      pixmap.save(filePath, 'JPEG')
    
    
  def exportPng(self, filePath=None):
  
    if not filePath:
      from gui.qtgui.FileSelect import selectSaveFile, FileType
      fileTypes = [FileType('PNG', ['*.png']),]
      filePath = selectSaveFile(self, caption='Select image file name',
                                directory=None, fileTypes=fileTypes)
                   
    if filePath:
      filePath = self._fixFileExtension(filePath, ['.png',])
      widget = self
      pixmap = QtGui.QPixmap.grabWidget(widget, widget.rect())
      pixmap.save(filePath, 'PNG')
  
  def _getRowCol(self, grid):
    
    layout = self.parent().layout()
        
    if grid:
      row, col = grid
      
      if row is None:
        row = layout.rowCount()
 
      if col is None:
        col = 0
  
    else:
      row = layout.rowCount()
      col = 0
     
    return row, col
  
  #def grid(self, row=None, col=0, rowSpan=None, colSpan=None, sticky=None):
  #  """ Maybe add Yucky Tkinter emulation """
    
  #  row, col, align = self._getGridData((row, col), sticky)
    
    
class Icon(QtGui.QIcon):

  def __init__(self, image=None, color=None, colorSize=22):
    
    assert image or color
    
    if color:
      image = QtGui.QPixmap(colorSize, colorSize)
      painter = QtGui.QPainter(image)
      
      if isinstance(color, QtGui.QColor):
        image.fill(color)
      
      elif isinstance(color, (tuple, list)):
        image.fill(color[0][:7])
        dx = colorSize/float(len(color))
        
        x = dx
        for i, c in enumerate(color[1:]):
          col = QtGui.QColor(c[:7])
          painter.setPen(col)
          painter.setBrush(col)
          painter.drawRect(x,0,x+dx,colorSize-1)
          x += dx
        
      else:  
        color = QtGui.QColor(color[:7])
        image.fill(color)
      
      painter.setPen(QtGui.QColor('#000000'))
      painter.setBrush(QtGui.QBrush())
      painter.drawRect(0,0,colorSize-1,colorSize-1)
      painter.end()
    
    elif not isinstance(image, QtGui.QIcon):
    
      if not path.exists(image):
        image = path.join(ICON_DIR, image)
      
    QtGui.QIcon.__init__(self, image)
    
    
  
