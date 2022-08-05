from PySide import QtCore, QtGui

from gui.qtgui.Base import Base

class Text(QtGui.QTextEdit, Base):

  def __init__(self, parent, text='', callback=None, wrap=False, readOnly=False, **kw):
    
    QtGui.QTextEdit.__init__(self, text, parent)
    Base.__init__(self, parent, **kw)
    
    if wrap:
      self.setLineWrapMode(self.WidgetWidth)
    else:
      self.setLineWrapMode(self.NoWrap)
    
    self.setReadOnly(readOnly)
    self.callback = callback
    self._isUserModified = False
   
    self.connect(self, QtCore.SIGNAL('textChanged()'), self._textChanged)
  
  def _textChanged(self):
  
    self._isUserModified = True
  
  def leaveEvent(self, event):
  
    QtGui.QPlainTextEdit.leaveEvent(self, event)
    
    if self._isUserModified and self.callback:
      self.callback()
      self._isUserModified = False
    
  def clearText(self):
    
    self.clear()
    self._isUserModified = False

  def getText(self):

    return self.toPlainText()

  def getHtml(self):

    return self.toHtml()

  def setText(self, text):
    
    self.setText(text)
    self._isUserModified = False

  def setHtml(self, text):
    
    QtGui.QTextEdit.setHtml(self, text)
    self._isUserModified = False

  def setPlainText(self, text):
    
    QtGui.QTextEdit.setPlainText(self, text)
    self._isUserModified = False
    
  def addText(self, text):
  
    self.append(text)
    self._isUserModified = False

import sys
from code import InteractiveConsole
from io import StringIO

class HelpWrapper():

  def __call__(self, *args, **kwds):
      if args or kwds:
          return help(*args, **kwds)

  def __repr__(self):
      return help.__repr__()
  
  
class Console(QtGui.QPlainTextEdit, Base):
  
  def __init__(self, parent, message='', context=None, prompt='>>> ', closeFunc=None, **kw):
    
    major, minor, micro, rel, serial = sys.version_info
    startup = 'Python %d.%d.%d on %s' % (major, minor, micro, sys.platform)
    
    QtGui.QPlainTextEdit.__init__(self, startup, parent)
    Base.__init__(self, parent, **kw)

    self.setStyleSheet("background-color: rgb(64, 64, 64); color: rgb(255, 255, 0);")    
    
    self.setOverwriteMode(False)
    self.setLineWrapMode(self.WidgetWidth)
    self.setReadOnly(False)
    self.setUndoRedoEnabled(False)
    self.document().setDefaultFont(QtGui.QFont("monospace", 10, QtGui.QFont.Normal))
    
    self.prompt = prompt
    
    if context:
      self.setContext(context)
    else:
      self.shell = InteractiveConsole(locals=locals())
    
    self.message = message
    self.closeFunc = closeFunc
    self.history = []
    self.historyPos = 0
    self._pos = 0
    self.promptLen = len(prompt)
    self.inBlock = False
    self.buffer = []
        
    self._startMessage()
   
  def _startMessage(self):
  
    if self.message:
      self.write(self.message)
      self.flush()
  
    self._prompt()
  
  def setContext(self, context):
  
    if context:
      if isinstance(context, dict):
        context = dict(context)
      else:
        context = dict([(n, getattr(context, n)) for n in dir(context)])
    else:
      context = locals()
    
    context['help'] = HelpWrapper()
    context['clear'] = self.clearConsole
    
    self.shell = InteractiveConsole(locals=context)  
  
  def _prompt(self):
  
    if self.inBlock:
      prompt = '.' * (self.promptLen-1)
      prompt+= ' '
    else:
      prompt = self.prompt
    
    self.appendPlainText(prompt)
    self.moveCursor(QtGui.QTextCursor.End)
    self._pos = self.textCursor().position()
  
  def flush(self):
  
    text = ''.join(self.buffer).rstrip()
    
    if text:
      self.appendPlainText(text)
      
    self.buffer = []
  
  def write(self, text):
    
    if text:
      self.buffer.append(text)
                
  def _issueLine(self):
  
    text = self.getLastLine()
    
    if text:
      self.appendHistory(text)
    
    stdout = sys.stdout
    stderr = sys.stderr
    
    sys.stdout = self
    sys.stderr = self

    self.inBlock = self.shell.push(text)
       
    sys.stdout = stdout
    sys.stderr = stderr
    
    self.flush()
    
    QtCore.QCoreApplication.processEvents()
      
    self._prompt()    
    
  def keyPressEvent(self, event):
    
    key = event.key()
    
    pos = self.textCursor().position()
    
    if pos < self._pos:
      self.moveCursor(QtGui.QTextCursor.End)
      
    if key in (QtCore.Qt.Key_Return, QtCore.Qt.Key_Enter):
      self._issueLine()
      return
    
    elif key == QtCore.Qt.Key_Home:
      self.setLinePos(0)
      return
      
    elif key == QtCore.Qt.Key_PageUp:
      return
        
    elif key in (QtCore.Qt.Key_Left, QtCore.Qt.Key_Backspace):
      if self.getLinePos() == 0:
        return
          
    elif key == QtCore.Qt.Key_Up:
      self.prevHistory()
      return
        
    elif key == QtCore.Qt.Key_Down:
      self.nextHistory()
      return
        
    elif key in (QtCore.Qt.Key_D, QtCore.Qt.Key_Z) and event.modifiers() == QtCore.Qt.ControlModifier:
      self.clearConsole()
      self._startMessage()
      
      if self.closeFunc:
        self.closeFunc()
       
    elif key == QtCore.Qt.Key_L and event.modifiers() == QtCore.Qt.ControlModifier:
      self.clearConsole()
      self._prompt()
     
    QtGui.QPlainTextEdit.keyPressEvent(self, event)
    
  def getLastLine(self):
  
    document = self.document()
    last = document.lineCount() - 1
    text = document.findBlockByLineNumber(last).text()
    text = text.rstrip()
    
    return u'%s' % (text[self.promptLen:])
    
  def setLastLine(self, text):

    self.moveCursor(QtGui.QTextCursor.End)
    self.moveCursor(QtGui.QTextCursor.StartOfLine, QtGui.QTextCursor.KeepAnchor)
    for i in range(self.promptLen):
      self.moveCursor(QtGui.QTextCursor.Right, QtGui.QTextCursor.KeepAnchor)
      
    self.textCursor().removeSelectedText()
    self.textCursor().insertText(text)
    self.moveCursor(QtGui.QTextCursor.End)
    
  def appendHistory(self, text):
  
     if text:
       if self.history:
         if self.history[-1] != text:
           self.history.append(text)
       
       else:
         self.history.append(text)
       
     self.historyPos = len(self.history)

  def prevHistory(self):
  
    if self.history:
      self.historyPos -= 1
      self.historyPos = max(0, self.historyPos)
      
      self.setLastLine( self.history[self.historyPos] )

  def nextHistory(self):
  
    if self.history:
      n = len(self.history)
      self.historyPos += 1
      self.historyPos = min(n-1, self.historyPos)
      
      if self.historyPos < n:
        self.setLastLine( self.history[self.historyPos] )

  def getLinePos(self):
    
    return self.textCursor().columnNumber() - self.promptLen
    
  def setLinePos(self, p):
    
    self.moveCursor(QtGui.QTextCursor.StartOfLine)
    
    for i in range(self.promptLen + p):
      self.moveCursor(QtGui.QTextCursor.Right)
 
  def clearConsole(self):
    
    self.clear()
    self.shell.resetbuffer()
    

if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.Menu import Menu
  
  import string
  
  def callback():
    print("Editing done")
  
  app = Application()

  window = QtGui.QWidget()
  
  text = Text(window, 'Initial text', callback)
  
  text.addText('Pugh,\nPugh,\nBarney McGrew,\nCuthbert,\nDibble,\nGrubb.\n')

  text.addText('\n\n%s\n' % (string.letters))
  

  text2 = Text(window, 'Initial text', callback, wrap=True)
  text2.setText('With wrapping\n%s\n' % (string.letters))
  
  console = Console(window, '[Python Console Widget]', app)
  
  window.show()
  
  app.start()

