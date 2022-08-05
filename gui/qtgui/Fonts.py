from PySide import QtCore, QtGui

from gui.qtgui.Base import Base, Icon
from gui.qtgui.Menu import Menu
from gui.qtgui.Button import Button

COMMON = ('Arial',
          'Times New Roman',
          'Monospace',
          'Georgia',
          'Courier 10 Pitch',
          'Verdana')

# Lots of font menu lookups currently segfault ! ! ! !
# - PySide v0.4.1 under Ubuntu 10.10
#
# Due to:
# "Fail to add dynamic slot to QObject
#  PySide support at most 50 dynamic slots.#
#
# Should be fixed by v1.0

LATIN = QtGui.QFontDatabase.Latin

# Original gui.FontList had modes for
# OpenGl, Tk, PostScript ect

# Font tree?

# Consider QFontComboBox -> "FontPanel"

class FontDialog(QtGui.QFontDialog):

  def __init__(self, parent=None, **kw):
    
    QtGui.QFontDialog.__init__(self, parent)
    self.fontDb = QtGui.QFontDatabase()
    
  def setFont(self, font):
    # color can be name, #hex, (r,g,b) or (r,g,b,a)

    self.setCurrentFont(font)
  
  def getFont(self):
  
    font, wasNotCanceled = QtGui.QFontDialog.getFont(self)
    
    if wasNotCanceled:
      return font
  
  def getFontName(self):
    
    font = self.getFont()
    
    if font:
      data = (font.family(), self.fontDb.styleString(font), font.pointSize())
      return '%s %s %d' % data


class FontPulldown(Button):

  def __init__(self, parent, text='Select...', minSize=8, maxSize=16, 
               writingSystem=LATIN, families=COMMON,
               showBold=False, showItalic=False, **kw):
    
    Button.__init__(self, parent, text, **kw)
    
    self.menuCallback = self.callback
    self.setCallback(None)
    self.fontDb = QtGui.QFontDatabase()
    
    fontMenu = FontMenu(self, text, minSize, maxSize, writingSystem,
                        families, showBold, showItalic, self._menuCallback)
    
    self.setMenu(fontMenu)
  
  def _menuCallback(self, font):
     
    text = '%s %s %d' % (font.family(),
                         self.fontDb.styleString(font),
                         font.pointSize())
    
    self.setText(text)
    
    if self.menuCallback:
      self.menuCallback(font)
    

class FontMenu(Menu):

  
  def __init__(self, parent=None, text='', minSize=8, maxSize=16, 
               writingSystem=LATIN, families=COMMON, default=None,
               showBold=False, showItalic=False, callback=None, **kw):

    Menu.__init__(self, parent, text, **kw)
    
    if self.icon().isNull():
      self.setIcon(Icon('icons/preferences-desktop-font.png'))
    
    self.fontDb = QtGui.QFontDatabase()
    self.fontCallback = callback
    self.callback = self._fontCallback
    
    self.default = default
    
    allFamilies = self.fontDb.families(writingSystem)
    
    families = set(families) & set(allFamilies)
    families = list(families)
    families.sort()

    for family in families:
      familyMenu = Menu(self, family)
      styles = self.fontDb.styles(family)

      filterStyles = []
      for style in styles:
        if not showItalic:
          if 'Italic' in style:
            continue
          if 'Oblique' in style:
            continue
            
        if not showBold:
          if 'Bold' in style:
            continue
                        
        sizes = self.fontDb.pointSizes(family, style)
        sizes = [x for x in sizes if minSize <= x <= maxSize]
        
        if not sizes:
          continue
          
        filterStyles.append((style, sizes))
      
      if not filterStyles:
        continue
        
         
      for style, sizes in filterStyles:  
        if len(filterStyles) > 1:
          subMenu = Menu(familyMenu, style)

        else:
          subMenu = familyMenu
          
        for size in sizes:
          text = '%d' % (size)
          data = (family, style, size)
          
          subMenu.addItem(text, self._fontCallback, object=data)
          
    if default:
      self.addItem('Default', self._fontCallback, object=default)

  def _fontCallback(self, obj):
    
    if self.fontCallback and obj:
      
      if obj == self.default:
        font = obj
        
      else:
        family, style, size = obj
        font = self.fontDb.font(family, style, size)
 
      if font:
        self.fontCallback(font)

    
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.MainWindow import MainWindow
  def callback(font):
    print(font)

  def openDialog():
    dialog = FontDialog()
    print(dialog.getFontName()) # opens the dialog
    print(dialog.selectedFont())

  app = Application()
  
  mainWindow = MainWindow(title='Font Test')
  menuBar = mainWindow.menuBar()
  
  fontMenu = FontMenu(menuBar, callback=callback)        
  
  FontPulldown(mainWindow.mainFrame, callback=callback)
  
  Button(mainWindow.mainFrame, 'Open Font Dialog', callback=openDialog)
  
  app.start()

