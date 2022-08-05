from PySide import QtCore, QtGui

from gui.qtgui.Base import Base
from gui.qtgui.FileSelect import selectSaveFile

class PrintDialog(QtGui.QPrintDialog):

  def __init__(self, parent=None, printer=None):
    
    QtGui.QPrintDialog.__init__(self, parent)
    
    self.setOption(QtGui.QAbstractPrintDialog.PrintToFile)
    
  def getPrinter(self):
  
    if self.exec_():
      return self.printer()

  def printWidget(self, widget):
  
    if self.exec_():
      printer = self.printer()
    
      if not printer:
        return

      painter = QtGui.QPainter()
      painter.begin(printer)
      
      widget.render(painter)
      
      painter.end()

  def printWidgetPdf(self, widget):
  
    if self.exec_():
      printer = self.printer()
      if not printer:
        return
      
      filePath = selectSaveFile(self, 'Select PDF file name',
                                fileTypes= ['*.pdf'])
      if not filePath:
        return
      
      printer.setOutputFileName(filePath)
      printer.setOutputFormat(printer.PdfFormat)

      painter = QtGui.QPainter()
      painter.begin(printer)
      
      widget.render(painter)
      
      painter.end()
      
  def printWidgetPixmap(self, widget):
 
    if self.exec_():
      printer = self.printer()
    
      if not printer:
        return
      
      painter = QtGui.QPainter()
      painter.begin(printer)
      #printer.newPage()
      
      pixmap =  QtGui.QPixmap.grabWidget(widget)
      painter.drawPixmap(0, 0, pixmap)
      painter.end()
  
def getFilePrinter(fileName=None, format='PDF', isHighRes=True,
                   isLandscape=False, isColor=True):
  
  if isHighRes:
    mode = QtGui.QPrinter.HighResolution
  else:
    mode = QtGui.QPrinter.ScreenResolution

  printer = QtGui.QPrinter(mode)
  
  if isColor:
    printer.setColorMode(QtGui.QPrinter.Color)
  else:
    printer.setColorMode(QtGui.QPrinter.GrayScale)
  
  if isLandscape:
    printer.setOrientation(QtGui.QPrinter.Landscape)
  else:
    printer.setOrientation(QtGui.QPrinter.Portrait)
  
  if format.upper() == 'PDF':
    printer.setOutputFormat(QtGui.QPrinter.PdfFormat)
    if not fileName:
      fileName = 'Output.pdf'
    
  else:
    printer.setOutputFormat(QtGui.QPrinter.PostScriptFormat)
    if not fileName:
      fileName = 'Output.ps'
    
  printer.setOutputFileName(fileName) 
  
  return printer

  
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.Graph import Graph, GraphAxis, GraphDataSet
  from gui.qtgui.Button import Button
  
  app = Application()
  
  window = BasePopup()
  window.setSize(600, 450)
  
  # Something to print
    
  numbers = [float(x/10.0) for x in range(1,21)]  
  squared = GraphDataSet([(x,x*x) for x in numbers], u'x\xb2', '#800000')
  inverse = GraphDataSet([(x,1.0/x) for x in numbers], u'x\u207b\xb9', '#008000',)
   
  xAxis = GraphAxis('X Value', labels=None, ticks=True)
  yAxis = GraphAxis('Y Value', labels=None, ticks=True, valRange=(0,12))
  
  graph = Graph(window, (xAxis, yAxis), (squared, inverse),
                title='Example Graph')  
  
  button = Button(window, '* Print This *')
  
 
  def printFunc():
  
    # Get the print options
    # Does nothing until called into action
 
    dialog = PrintDialog(window)
 
    # Do the printing for the widget
 
    dialog.printWidgetPixmap(graph)
  
  
  
  button.setCallback(printFunc)
  
  app.start()
