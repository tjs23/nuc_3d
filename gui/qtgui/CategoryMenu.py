from PySide import QtCore, QtGui

from gui.qtgui.Button import Button
from gui.qtgui.Menu import Menu

class CategoryMenu(Button):

  def __init__(self, parent, texts=None, objects=None, categories=None,
               icons=None, callback=None, index=0, **kw):
    
    Button.__init__(self, parent, **kw)
    
    self.callback = callback
    
    menu = Menu(self, callback=self._callback)
    self.setMenu(menu)   
    self.menu = menu
    
    self.index = None
    self.selected = None
    self.object = None   
    
    self.texts = []
    self.objects = []
    self.categories = []
     
    self.setStyleSheet("text-align: left; padding: 4px;")
    self.setData(texts, objects, categories, index, icons)
    
  def _callback(self, obj):
  
    self.object = obj
    self.index = list(self.objects).index(obj)
    self.selected = self.texts[self.index]
    self.setText(self.selected)
    
    if self.callback:
      self.callback(obj)
  
  def setCurrentIndex(self, index):
   
    self.index = index
    if self.texts:
      self.selected = self.texts[self.index]
      self.object = self.objects[self.index]
      self.setText(self.selected)
 
    
  def get(self):
  
    return self.currentObject()

  def currentObject(self):
    
    return self.object

  def currentData(self):

    return (self.selected, self.object)

  def select(self, item):

    index = None
    
    if item in self.texts:
      index = list(self.texts).index(item)

    elif item in self.objects:
      index = list(self.objects).index(item)

    if index is not None:
      self.setCurrentIndex(index)

  def set(self, item):

    self.select(item)
    
  def setData(self, texts=None, objects=None, categories=None,
              index=None, icons=None):
    
    texts = texts or []
    objects = objects or []
    categories = categories or []
    
    n = len(texts)
    
    if objects:
      msg = 'len(texts) = %d, len(objects) = %d'
      assert n == len(objects), msg % (n, len(objects))
      
    else:
      objects = texts[:]
    
    if icons:
      while len(icons) < n:
        icons.append(None)
        
    else:
      icons = [None] * n
    
    if categories:
      while len(categories) < n:
        categories.append(None)
        
    else:
      categories = [None] * n
    
    self.menu.clear()
    
    uniqCats = set()
    for cat in categories:
      if isinstance(cat, (set, list, tuple)):
        uniqCats.update(cat)
      else:
        uniqCats.add(cat)
    
    uniqCats = list(uniqCats)
    uniqCats.sort()
    
    subMenuDict = {None:self.menu}
    for i, text in enumerate(texts):
      cats = categories[i]
      
      if not isinstance(cats, (set, list, tuple)):
        cats = [cats]
      
      for cat in cats:
        if cat is None:
          subMenuDict[cat].addItem(text, callback=None,
                                   object=objects[i],
                                   icon=icons[i])
    
    for cat in uniqCats:
      if cat is not None:
        subMenuDict[cat] = Menu(self.menu, text=cat)
    
    # Add cascade after  
    for i, text in enumerate(texts):
      cats = categories[i]
      
      if not isinstance(cats, (set, list, tuple)):
        cats = [cats]
    
      for cat in cats:
        if cat is None:
          continue
        
        subMenuDict[cat].addItem(text, callback=None,
                                 object=objects[i],
                                 icon=icons[i])
    
    self.texts = texts
    self.objects = objects
    self.categories = categories
    self.icons = icons
    
    if index is not None:  
      self.setCurrentIndex(index)

    
if __name__ == '__main__':

  from gui.qtgui.Application import Application
  from gui.qtgui.Menu import Menu
  from gui.qtgui.Base import SizePolicy

  app = Application()

  window = QtGui.QWidget()

  def clickObj(obj):
    print("Selected", obj)
  
  texts = ['Jan','Feb','Mar',
           'Apr','May','Jun',
           'Jul','Aug','Sep',
           'Oct','Nov','Dec']
  cats = ['31',('28','29'),'31',
          '30','31','30',
          '31','31','30',
          '31','30','31',]         
  b = CategoryMenu(window, texts=texts, categories=cats,
                   callback=clickObj,
                   tipText='A two-level cascase',
                   grid=(0, 4), vPolicy=SizePolicy.expanding)
  
  window.show()
  
  app.start()

