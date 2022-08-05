from PySide import QtCore, QtGui

from gui.qtgui.Base import Base, Icon
from gui.qtgui.BasePopup import BasePopup
from gui.qtgui.Button import Button
from gui.qtgui.WebBrowser import WebBrowser

class ButtonList(QtGui.QWidget, Base):

  def __init__(self, parent=None, texts=None, callbacks=None, icons=None,
               tipTexts=None, direction='h', commands=None, images=None, **kw):

    QtGui.QWidget.__init__(self, parent)
    Base.__init__(self, parent,  **kw)

    if commands:
      print("qtgui.ButtonList.commands is deprecated; use .callbacks")
      callbacks = commands

    if images:
      print("qtgui.ButtonList.images is deprecated; use .icons")
      icons = images

    if texts is None:
      texts = []

    if callbacks is None:
      callbacks = []

    assert len(texts) == len(callbacks)

    direction = direction.lower()
    self.direction = direction
    
    if tipTexts is None:
      tipTexts = []
      
    if icons is None:
      icons = []

    while len(tipTexts) < len(texts):
      tipTexts.append(None)  

    while len(icons) < len(texts):
      icons.append(None)  

    self.buttons = []
    self.addButtons(texts, callbacks, icons, tipTexts)
  
  def addButtons(self, texts, callbacks, icons=None, tipTexts=None):
  
    if tipTexts is None:
      tipTexts = []
      
    if icons is None:
      icons = []

    while len(tipTexts) < len(texts):
      tipTexts.append(None)  

    while len(icons) < len(texts):
      icons.append(None)  
    
    j = len(self.buttons)
    for i, text in enumerate(texts):
      if 'h' in self.direction:
        grid = (0,i+j)
      else:
        grid = (i+j,0)
        
      button = Button(self, text, callbacks[i], icons[i],
                      tipText=tipTexts[i], grid=grid)
      button.setMinimumWidth(20)
      self.buttons.append(button)
    
      
class UtilityButtonList(ButtonList):

  def __init__(self, parent,
               webBrowser=None, helpUrl=None, helpMsg=None,
               doClone=True, doHelp=True, doClose=True,
               cloneText=None, helpText=None, closeText=None,
               cloneCmd=None, helpCmd=None, closeCmd=None,
               cloneTip='Duplicate winddow', helpTip='Show help', closeTip='Close window',
               *args, **kw):
    
    ButtonList.__init__(self, parent, *args, **kw)
    
    self.helpUrl = helpUrl
    self.helpMsg = helpMsg
    
    self.popup = parent.window()
    if not isinstance(self.popup, BasePopup):
      self.popup = None
   
    if self.popup and not webBrowser:
      webBrowser = WebBrowser(self.popup)
 
    self.webBrowser = webBrowser

    _callbacks = []
    _texts    = []
    _icons   = []
    _tipTexts = []
    
    _doActions = [(doClone, cloneCmd, self.clonePopup, cloneText, 'icons/window-duplicate.png', cloneTip),
                  (doHelp,   helpCmd, self.launchHelp, helpText,  'icons/system-help.png',       helpTip),
                  (doClose, closeCmd, self.closePopup, closeText, 'icons/window-close.png',     closeTip),]
    
    for doAction, userCmd, defaultCmd, text, image, tipText in _doActions:
      if doAction:
        _callbacks.append(userCmd or defaultCmd)
        _tipTexts.append(tipText)
        _texts.append(text)
        
        if image:
          icon = Icon(image)
        else:
          icon = None  
        
        _icons.append(icon)       

    self.addButtons(_texts, _callbacks, _icons, _tipTexts)

  def clonePopup(self):
     
    if self.popup:
      try:
        newPopup = self.popup.__class__(self.popup.parent)
        x,y, w, h = self.getGeometry()
        newPopup.setGeometry(x+25, y+25, w, h)
      
      except:
        pass
  
  def launchHelp(self):
  
    if self.helpUrl and self.webBrowser:
      self.webBrowser.open(self.helpUrl)
      
    elif self.popup:
      from gui.qtgui.WebView import WebViewPopup
      WebViewPopup(self.helpMsg or 'http://www.ccpn.ac.uk/documentation')  
  
  def closePopup(self):
  
    if self.popup:
      self.popup.close()
      
    else:
      self.destroy()  
  

if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup

  def callback(text):
    print('callback', text)

  texts = ['abc', 'def', 'ghi']
  callbacks = [lambda t=text: callback(t) for text in texts]
  icons = [None, None, 'icons/applications-system.png']

  app = Application()
  popup = BasePopup(title='Test ButtonList')
  popup.setSize(200,60)
  buttons = ButtonList(parent=popup, texts=texts, callbacks=callbacks, icons=icons)
  utils = UtilityButtonList(parent=popup, texts=texts, callbacks=callbacks)
  
  app.start()

