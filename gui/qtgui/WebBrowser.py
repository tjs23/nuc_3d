import webbrowser as wb
from gui.qtgui.PulldownList import PulldownList
from gui.qtgui.WebView import WebViewPopup

browserNames = ['firefox','netscape','mozilla','konqueror','kfm','mosaic',
                'grail','w3m','windows-default','internet-config']

class WebBrowser:

  def __init__(self, parent, name=None, url=None):
  
    names = getBrowserList()
    if names and (not name):
      name = names[0]
    
    self.name   = name
    
    if url:
      self.open(url)
  
  def open(self, url):
    
    try:
      browser = wb.get(self.name)
      browser.open(url)
      
    except:
      WebViewPopup(url)
  
class WebBrowserPulldown(PulldownList):

  def __init__(self, parent, browser=None, **kw):

    self.browserList = getBrowserList()

    if not browser:
      browser = getDefaultBrowser()

    if self.browserList:
      if (not browser) or (browser not in self.browserList):
        browser = self.browserList[0]

    self.browser = browser

    PulldownList.__init__(self, parent, **kw)

    if self.browserList:
      self.setup(self.browserList, self.browserList, self.browserList.index(self.browser))
    
    self.callback = self.setWebBrowser

  def setWebBrowser(self, name):

    if name != self.browser:
      self.browser = name

  def destroy(self):

    pass

def getBrowserList():

  browsers = []
  default  = getDefaultBrowser()
  if default:
    browsers = [default,]
  
  for name in browserNames:
    if name == default:
      continue
  
    try:
      wb.get(name)
      browsers.append(name)
    except:
      
      try:
        if wb._iscommand(name):
          wb.register(name, None, wb.Netscape(name))
          wb.get(name)
          browsers.append(name)
      except:
        continue

  return browsers

def getDefaultBrowser():

  try:
    br = wb.get()
  except:
    return
  
  if not hasattr(br, 'name'):
    # Max OS X
    return
  
  try:
    wb.get(br.name)
  except:
    wb.register(br.name, None, br)
  
  return br.name
