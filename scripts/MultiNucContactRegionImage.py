import sys

from math import ceil
from os.path import abspath, dirname, join, splitext
sys.path.append(abspath(dirname(dirname(__file__))))

from numpy import array, dstack, uint8, zeros, log
from glob import glob
from NucApi import Nucleus
from util.Image import pixmapToImage


def getContactRegionPixmap(nuc, groupName, chromo, start, end, binSize=int(5e4)):
  
  contacts = nuc.getContactMatrix(chromo, chromo, binSize, groupName, (start, end), (start, end))

  p1 = int(start) / binSize
  p2 = 1 + int(end) /binSize
  
  matrix = array(contacts, float)
  
  if nuc.areContactsSingleCell(groupName):
    matrix[matrix.nonzero()] = 1.0
  else:
    matrix = log(1.0 + matrix)
  
  maxVal = matrix.max()
  if maxVal:
    matrix /= maxVal 
  
  r = g = array(matrix)
  b = zeros(r.shape)
  
  return dstack([r, g, b])
  

def getDirContactRegionImages(popFile, dirPath, chromo='1', start=3e6, end=4e6, format='PNG'):

  filePaths = [popFile] + glob(join(dirName, '*10x_100kb.nuc'))
  pixmaps = []
  
  for i, filePath in enumerate(filePaths):
    print(i, filePath)
    
    fileRoot =  splitext(filePath)[0]
    nuc = Nucleus(filePath)
    
    groupName, group = nuc.getContactGroups()[0]
    
    #print('Make image')
    pixmap = getContactRegionPixmap(nuc, groupName, chromo, start, end)
    pixmaps.append(pixmap)
    
    if i == 0:
      mean = zeros(pixmap.shape, float)
    else:
      mean += pixmap
    
  mean /= mean.max()
  pixmaps.append(mean ** 0.5)
  
  nImg = len(pixmaps)
  
  n = len(pixmaps[0])
  
  nCols = 5
  nRows = int(ceil(nImg/float(nCols)))
  
  x = n * nCols + nCols + 1
  y = n * nRows + nRows + 1
  
  
  pixmap = zeros((y,x,3), uint8) + [32, 32, 32]
  
  for i, p in enumerate(pixmaps):
    row = int(i / nCols)
    col = i % nCols 
    
    a = 1 + row + row * n
    b = 1 + col + col * n
    
    pixmap[a:a+n,b:b+n,:] = array(255*p, uint8)
  
  image = pixmapToImage(pixmap)
   
  sc = 2
  h, w = image.size
  image = image.resize((h * sc, w * sc))
  
  fileRoot = 'DomainContacts'
  
  imgFileName = '%s_Chr%s:%.2f-%.2fMb.%s' % (fileRoot, chromo, start/1e6, end/1e6, format.lower())
  
  print('Save image "%s"' % imgFileName)
  image.save(imgFileName, format)
  

if __name__ == '__main__':
  
  popFile = '/home/tjs23/nucleus/SLX-7671_hapsort_pop_mm9.nuc'
  dirName = '/home/tjs23/nucleus/Nextera_01/best_100kb/'

  getDirContactRegionImages(popFile, dirName, '1', 72e6, 76e6)
