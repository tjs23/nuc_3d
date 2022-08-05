from PySide import QtCore, QtGui

from gui.qtgui.Base import Base

def selectFile(parent, caption='Select file to open',
               directory=None, fileTypes=None, returnType=False):
  
  dialog = FileDialog(parent, caption,
                      directory=directory, doSave=False,
                      fileTypes=fileTypes,
                      returnType=returnType)
                      
  return dialog.getFile()


def selectSaveFile(parent, caption='Select save file',
                   directory=None, fileTypes=None,
                   default=None):

  dialog = FileDialog(parent, caption,
                      directory=directory, doSave=True,
                      default=default, fileTypes=fileTypes)
  
  return dialog.getFile()

def selectFiles(parent, caption='Select files to open',
                directory=None, fileTypes=None, returnType=False):

  dialog = FileDialog(parent, caption, multiFile=True,
                      directory=directory, doSave=False,
                      fileTypes=fileTypes, returnType=False)
  
  return dialog.getFiles()


def selectDirectory(parent, caption='Select directory', directory=None):
  
  dialog = DirectoryDialog(parent, caption, directory=directory)
  
  return dialog.getDirectory()



class FileType:
  
  def __init__(self, text='All files', extensions=None, showExtensions=True):
    
    self.text = text.strip()
    self.extensions = extensions or ['*',]
    self.showExtensions = showExtensions
    
    for i, ext in enumerate(self.extensions):
      ext = ext.strip()
      self.extensions[i] = ext

  def getFilterText(self):
    
    if self.showExtensions:
      return '%s (%s)' % (self.text, ' '.join(self.extensions))
    
    else:
      return '%s' % self.text

      

class FileDialog(QtGui.QFileDialog):

  def __init__(self, parent, caption='Select file',
               directory=None, doSave=False,
               showDetails=False, multiFile=False, suffix=None,
               default=None, fileTypes=None, sideBarUrls=None,
               returnType=False, warnOverwrite=True):
  
    QtGui.QFileDialog.__init__(self, parent, caption, directory)
    
    self.default = default
    self.returnType = returnType
    
    if directory:
      self.setDirectory(directory)
    
    if doSave:
      if warnOverwrite:
        self.setAcceptMode(QtGui.QFileDialog.AcceptSave)

      self.setFileMode(QtGui.QFileDialog.AnyFile)
    
    else:
      self.setAcceptMode(QtGui.QFileDialog.AcceptOpen)
      
      if multiFile:
        self.setFileMode(QtGui.QFileDialog.ExistingFiles)
      else:
        self.setFileMode(QtGui.QFileDialog.ExistingFile)
    
    if showDetails:
      self.setViewMode(QtGui.QFileDialog.Detail)
    else:
      self.setViewMode(QtGui.QFileDialog.List)
    
    if suffix:
      self.setDefaultSuffix(suffix)   
    
    if fileTypes:
      filters = []
      
      for fileType in fileTypes:
        if isinstance(fileType, FileType):
          fileType = fileType.getFilterText()
        
        filters.append(fileType)  
      
      self.setNameFilters(filters)
      
    if sideBarUrls:
      urls = []
      
      for url in sideBarUrls:
        if not isinstance(QtCore.QUrl, url):
          url = QtCore.QUrl(url)
          
        urls.append(url) 
      
      self.setSideBarUrls(urls)
    
      
  def getDirectory(self):

    if self.exec_():
      fileNames = self.selectedFiles()
      if fileNames: # A dir was clicked
        return fileNames[0]
      
      qDir = self.directory()
      if qDir: # Current dir
        return qDir.path()
  
  
  def getFile(self, default=None):
    
    
    if default or self.default:
      self.selectFile(default or self.default)
    
    if self.exec_():
      fileNames = self.selectedFiles()
      
      if self.returnType:
        return fileNames[0], self.selectedFilter()
      else:
        return fileNames[0]
    
    
  def getFiles(self):
    
    if self.exec_():
      fileNames = self.selectedFiles()
      
      if self.returnType:
        return fileNames, self.selectedFilter()
      else:
        return fileNames
  
  # setDirectory() exists in Qt
  
  def setFile(self, filePath):
  
    self.default = filePath
    # Comes out when self.selectFile()  
  
  
  
class DirectoryDialog(FileDialog):
  
  def __init__(self, parent, caption='Select directory',
               directory=None, doSave=False,
               showDetails=False, showFiles=False, *kw):
                 
    FileDialog.__init__(self, parent, caption, directory=directory, doSave=doSave,
                        showDetails=showDetails, *kw)
    
    self.setFileMode(QtGui.QFileDialog.Directory)
    
    if not showFiles:
      self.setOptions(QtGui.QFileDialog.ShowDirsOnly)

    
    
class FileWidget(QtGui.QWidget, Base):

  def __init__(self, parent, caption='Select file',
               directory=None, doSave=False,
               showDetails=False, multiFile=False, suffix=None,
               default=None, fileTypes=None, sideBarUrls=None,
               callback=None, **kw):
  
    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, parent, **kw)
    
    self.callback = callback
    self.fileBrowser = FileDialog(self, caption, directory, doSave,
                             showDetails, multiFile, suffix,
                             default, fileTypes, sideBarUrls)
    
    self.fileBrowser.setParent(self, QtCore.Qt.Widget)
    self.fileBrowser.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
    self.fileBrowser.accept = self._callback
    self.fileBrowser.done = self._callback
    self.move(0, 0)
    self.show()
  
  
  def _callback(self, accept=True):
    
    if accept and self.callback:
      self.callback()
    
    
  def getDirectory(self):

    fileNames = self.fileBrowser.selectedFiles()
    if fileNames: # A dir was clicked
      return fileNames[0]
    
    qDir = self.fileBrowser.directory()
    if qDir: # Current dir
      return qDir.path()
    
    
  def getFile(self, default=None):
        
    fileNames = self.fileBrowser.selectedFiles()
    return fileNames[0]    
    
    
  def getFiles(self):
    
    fileNames = self.fileBrowser.selectedFiles()
    return fileNames
  
  
  def setDirectory(self, directory):
  
    self.fileBrowser.setDirectory(directory)
  
  
  def setFile(self, filePath):
  
    self.fileBrowser.selectFile(filePath)
    
    
    
if __name__ == '__main__':

  import sys
  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.Button import Button

  argv = sys.argv
  app = Application(argv)
  popup = BasePopup(title='Test FileSelect Directory')
  

  def directoryCallback():

    directory = selectDirectory(popup, 'Select a directory')

    print(directory)

  def fileCallback():
   
    fileName = selectFile(popup, 'Select a file', fileTypes=fileTypes)

    print(fileName)
  
  
  button = Button(popup, text='Choose directory', callback=directoryCallback)
  button = Button(popup, text='Choose file', callback=fileCallback)
  popup.setSize(650, 550)
  
  
  def callback():
    fileTypes = [FileType('Python', ['*.py','*.pyw']),]
  
    dialog = FileDialog(popup, caption='Select file',
                        directory=None, doSave=False,
                        showDetails=True, multiFile=True, suffix=None,
                        default=None, fileTypes=fileTypes, sideBarUrls=None)
    
    print(dialog.getFiles())
    
  
  button = Button(popup, text='Choose multiple files', callback=callback)
  
  fileWidget = FileWidget(popup, gridSpan=(1,2))
  
  def widgetCallback():
    print(fileWidget.getFile())
  
  fileWidget.callback = widgetCallback
    
  app.start()
