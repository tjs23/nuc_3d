import sys, os

if __name__ == '__main__':
  sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from PySide import QtGui
from gui.ImagePanel import ImagePanel

if __name__ == '__main__':
  app = QtGui.QApplication(["Test"])
  main = QtGui.QMainWindow()
  widget = ImagePanel()
  main.setCentralWidget(widget)
  main.show()
  app.exec_()
