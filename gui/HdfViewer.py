import sys, os
from gui.qtgui.Application import Application
from gui.qtgui.BasePopup import BasePopup
from gui.qtgui.FileSelect import selectFile, FileType
from gui.qtgui.Button import Button
from gui.qtgui.Table import ObjectTable, Column

from gui.HdfTree import HdfTreePanel
from h5py import File, Group, Dataset, AttributeManager
from numpy import cumprod, array, ndarray

dirname =  os.path.dirname
   
class HdfViewer(BasePopup):

  def __init__(self, parent):
  
    BasePopup.__init__(self, parent, title='HDF file viewer')

    self.tree = HdfTreePanel(self, callback=self.selectItem, grid=(0,0))
    self.fileDir = None
    self.data = None
    self.obj = None
    self.transpose = False
    
    but = Button(self, 'Select HDF file', callback=self.selectHdfFile, grid=(1,0))
    
    self.defaultCols = [Column('Index 0', self.getVal),
                        Column('Index 1', self.getVal),
                        Column('Value',   self.getVal),]
    self.table = ObjectTable(self, self.defaultCols, [], grid=(0,1))

    self.transpButton = Button(self, 'Transpose', callback=self.toggleTranspose, grid=(1,1))


  def toggleTranspose(self):
  
    self.transpose = not self.transpose
  
    if self.obj is not None:
      self.selectItem(self.obj)
  
  
  def getIndex(self, indices, col=0):
  
    return int(indices[col])
    
    
  def getVal(self, indices):
    
    if self.data is not None:
      return  str(self.data[indices])


  def getArrayVal(self, indices, col=0):
    
    if self.data is not None:
      row = self.data[indices[:-1]]
      return  str(row[col])


  def singular(self, data):
    
    return str(data)

    
  def selectItem(self, obj):
    
    self.obj = obj
    
    if hasattr(obj, 'dtype'):
      self.data = array(obj)
      
      columns = []
      shape   = self.data.shape
      dtype   = self.data.dtype
      
      if not shape:
        columns.append(Column('Value', self.singular))
        rows = [self.data]
        self.transpButton.disable()
      
      elif len(shape) == 1:
        columns.append(Column('Index', self.getIndex))
        columns.append(Column('Value', self.getVal))
        rows = [(i,) for i in range(shape[0])]
        self.transpButton.disable()
      
      else:
        self.transpButton.enable()
        if self.transpose:
          self.data = self.data.T
          shape = self.data.shape
      
        rows = []
        rowsAppend = rows.append
        ndim = self.data.ndim
        dims = range(ndim)
        
        if shape[-1] < 9:
          step = shape[-1]
          for i, dim in enumerate(shape[:-1]):
            dim = ndim-i-1 if self.transpose else i
            column = Column('Index %d' % dim, lambda x, c=i:self.getIndex(x, c))
            columns.append(column)

          dim = 0 if self.transpose else ndim - 1
          for i in range(shape[-1]):
            column = Column('Value %d[%d]' % (dim,i), lambda x, c=i:self.getArrayVal(x, c))
            columns.append(column)
          
        else:
          step = 1
          for i in dims:
            dim = ndim-i-1 if self.transpose else i
            column = Column('Index %d' % dim, lambda x, c=i:self.getIndex(x, c))
            columns.append(column)
 
          columns.append(Column('Value', self.getVal))
          
        strides = [0] * ndim

        n = 1
        for i in range(ndim-1, -1, -1):
          strides[i] = n
          n *= shape[i]
        
        for i in range(0, n, step):
          indices = [0] * ndim
          index = i
 
          for j in dims:
            indices[j], index = divmod(index, strides[j])
 
          rowsAppend(tuple(indices))
        
      self.table.setObjects([])
      self.table.setColumns(columns)
      self.table.setObjects(rows)
      self.table.resizeColumnToContents(0)
      
    
    else:
      self.transpButton.disable()
      self.data = None
      self.table.setObjects([])
      self.table.setColumns(self.defaultCols)
      
  def selectHdfFile(self):
  
    fileTypes = [FileType('Nuc files',  ['*.nuc']),
                 FileType('HDF5 files', ['*.hdf', '*.h5', '*.hdf5', '*.he5']),
                 FileType('All files',  ['*.*']),]
                 
    msg = 'Select HDF file'
    
    if not self.fileDir:
      self.fileDir = dirname(dirname(dirname(__file__)))
    
    filePath = selectFile(popup, msg, self.fileDir, fileTypes)
    
    if filePath:
      self.fileDir = os.path.dirname(filePath)
      self.tree.setFile(filePath)
      

if __name__ == '__main__':
   
  argv = sys.argv
  app = Application(argv[0])
  
  popup = HdfViewer(None)
  
  if len(argv) > 1:
    popup.tree.setFile(argv[1])
  
  app.start()
  
