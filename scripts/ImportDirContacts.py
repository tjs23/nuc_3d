import sys, gc
from os.path import abspath, dirname, join, splitext
sys.path.append(abspath(dirname(dirname(__file__))))

from glob import glob
from NucApi import Nucleus


def importDirContacts(dirPath):
  
  filePaths = glob(join(dirName, '*N706_S502_uniques.bam'))
  #filePaths += glob(join(dirName, '*.fend_pairs'))
  
  for filePath in filePaths:
    
    nucFile, fExt = splitext(filePath)
    nucFile += '.nuc'
    
    if fExt == '.bam':
      format = 'sam'
    else:
      format = 'text'
    
    nuc = Nucleus(nucFile)
    nuc.importContacts(filePath, format)
    nuc.save()
    
    nuc = None
    gc.collect()

if __name__ == '__main__':
  
  dirName = '/data/hi-c/sc2/'
  #dirName = '/home/tjs23/Desktop/Nextera_02'
  
  importDirContacts(dirName)
