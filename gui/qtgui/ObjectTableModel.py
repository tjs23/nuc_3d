from gui.qtgui.TableModel import TableModel

class ObjectTableModel(TableModel):

  # should be overridden in subclass
  header = None  # tuple or list of column headings
  getters = None  # gets data from object and puts in table

  # can be overridden in subclass
  setters = None  # gets data from edited table and puts in object
  colsUsed = None  # the columns in the header that are used

  def __init__(self, objects=None, header=None, colsUsed=None):

    if not objects:
      objects = []

    TableModel.__init__(self)

    self.objects = objects

    if header:
      self.header = header

    if colsUsed:
      self.colsUsed = colsUsed
    elif not self.colsUsed:
      self.colsUsed = range(len(self.header))

  ############################################################
  # implementation
  ############################################################

  def numberRows(self):
    return len(self.objects)

  def numberCols(self):
    return len(self.colsUsed)

  def headerForCol(self, col):
    actualCol = self.colsUsed[col]
    return self.header[actualCol]

  def dataForCell(self, row, col):
    actualCol = self.colsUsed[col]
    getter = self.getters[actualCol]
    object = self.objects[row]
    return getter(self, object)

  def setDataForCell(self, row, col, value):
    if self.setters:
      actualCol = self.colsUsed[col]
      setter = self.setters[actualCol]
      if setter:
        object = self.objects[row]
        return setter(self, object, value)

    return False

  def sortRows(self, col, isDescending=False):
    actualCol = self.colsUsed[col]
    getter = self.getters[actualCol]
    keyFunc = lambda object: getter(self, object)
    objects = sorted(self.objects, key=keyFunc)
    if isDescending:
      objects.reverse()
    self.objects = objects

  def isEditableCell(self, row, col):
    actualCol = self.colsUsed[col]
    return self.setters and self.setters[actualCol]

if __name__ == '__main__':

  import sys

  from gui.qtgui.Application import Application
  from gui.qtgui.BasePopup import BasePopup
  from gui.qtgui.Table import Table

  class Fruit:
    def __init__(self, name, weight):
      self.name = name
      self.weight = weight

  class FruitTableModel(ObjectTableModel):

    header = ('fruit', 'weight')

    def getName(self, fruit):
      return fruit.name

    def setName(self, fruit, name):
      fruit.name = name
      return True

    def getWeight(self, fruit):
      return fruit.weight

    def setWeight(self, fruit, weight):
      try:
        fruit.weight = int(weight)
        return True
      except:
        return False

    getters = (getName, getWeight)
    setters = (setName, setWeight)

  objects = (Fruit('apple', 50), Fruit('pear', 200), Fruit('orange', 300))
  model = FruitTableModel(objects)

  app = Application(sys.argv)
  popup = BasePopup(title='Test Table')
  table = Table(parent=popup, model=model)
  table.setMinimumSize(400, 300)
  app.start()


