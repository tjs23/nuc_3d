
from PySide import QtCore, QtGui, QtWebKit

from gui.qtgui.BasePopup import BasePopup
from gui.qtgui.Base import Base

# # # #  T D  # # # # 
# Back
# Forward
# Reload
# Find
# Show plain text
#
# Remove HelpPopup

class WebViewPanel(QtWebKit.QWebView, Base):

  def __init__(self, parent, **kw):
  
    QtWebKit.QWebView.__init__(self, parent=parent)
    Base.__init__(self, parent, **kw)
    
class WebViewPopup(BasePopup):

  def __init__(self, url=None, **kw):
  
    BasePopup.__init__(self, title='Web View', **kw)

    self.webViewPanel = WebViewPanel(self)
    
    if url:
      self.setUrl(url)
   
  def setUrl(self, urlText):
  
    qUrl = QtCore.QUrl(urlText)
    self.webViewPanel.setUrl(qUrl)

if __name__ == '__main__':
  
  from gui.qtgui.Application import Application
  
  app = Application()

  popup = WebViewPopup('http://www.ccpn.ac.uk')
  
  app.start()
