
def main():

  import sys
  
  from PySide import QtCore, QtGui
  
  from gui.NucMain import NucMain

  qtApp = QtGui.QApplication(('Nuc3D',))

  QtCore.QCoreApplication.setApplicationName("Nuc3D")

  filePath = None
  if len(sys.argv) > 1:
    fileNames = [fp for fp in sys.argv[1:] if fp.endswith('.nuc')]
    
    if fileNames:
      filePath = fileNames[0]
    

  main = NucMain(None, filePath)
  main.show()
  
  sys.exit(qtApp.exec_())


if __name__ == '__main__':
  
  main()
