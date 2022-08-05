from PySide import QtCore, QtGui
from os import path

from gui.qtgui.Tree import FileSystemTreePanel
from gui.qtgui.MessageDialog import showError
from gui.qtgui.FileSelect import FileType
from gui.qtgui.Frame import Frame
from util.Io import checkRegularFile
from formats.Util import DATA_TRACK_FORMATS


class FileSystemPanel(Frame):

  def __init__(self, parent, openFunc):
  
    Frame.__init__(self, parent)
    
    fileTypes = [FileType('Nuc3D file', ['*.nuc']),
                 FileType('Nuc contacts', ['*.pfe','*.fend_pairs']),
                 FileType('BAM paired contacts', ['*.bam', '*.sam'])]

    for format in sorted(DATA_TRACK_FORMATS):
      fileTypes.append( FileType(format, DATA_TRACK_FORMATS[format]))
    
    fileTypes.append( FileType('All files', ['*',]))

    self.setObjectName('FileSystemPanel')
    
    self.openFunc = openFunc
    self.treePanel = FileSystemTreePanel(self, fileTypes, callback=self.openFile)
    self.treePanel.fileTypePulldown.setObjectName('FileTypePulldown')

    
  def openFile(self, filePath, haveModKey=False):
    
    replace = not haveModKey
  
    if path.isdir(filePath):
      return
 
    isOk, msg = checkRegularFile(filePath)
 
    if not isOk:
      showError('Could not load file', msg, parent=self)
      return
      
    self.openFunc(filePath, replace)
 
    
  def openDir(self, dirPath):
  
    self.treePanel.openDir(dirPath)

