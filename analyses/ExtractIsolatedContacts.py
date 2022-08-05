import sys, os
from os.path import dirname, join

nucPath = dirname(dirname(__file__))
sys.path.append(nucPath)

from NucApi import Nucleus

def extractIsolatedContacts(bamFileName, threshold=int(1e6), groupName='singleCell'):

  print('Processing "%s"...' % (bamFileName,))
  
  outFileName = bamFileName[:-4]+'_iso.txt'
  
  # Temporary .nuc file
  nucFileName = bamFileName+'_temp.nuc'
  if os.path.exists(nucFileName):
    os.remove(nucFileName)
  
  # Make nucleus object and import contacts
  nuc = Nucleus(nucFileName)
  nuc.importContactFile(bamFileName, format='sam', groupName=groupName)
  
  numContacts = nuc.getNumContacts(groupName)
  if numContacts > 5e4:
    outFileObj = open(outFileName, 'w')
    
    # Get dicts contating all contacts and indices of thoseisolated 
    contactDict = nuc.getContacts(groupName)
    indicesDict = nuc.getIsolatedContacts(groupName, threshold=threshold)
    
    for key in indicesDict:
      chrA, chrB = key
      indices = indicesDict[key]
      contacts = contactDict[key].T # Contacts are shape ((posA, posB, nObs), nContacts), so transpose
      isolated = contacts[indices]
      
      for posA, posB, numObs in isolated:
        line = '%s\t%d\t%s\t%d\n' % (chrA, posA, chrB, posB)
        outFileObj.write(line)
    
    outFileObj.close()
  
  os.remove(nucFileName)
  print(' ...done')


if __name__ == '__main__':

  dirName = join(nucPath, 'data', 'SiCUP', 'sample_1527', 'merged_processing')
  filePaths = [join(dirName, fn) for fn in os.listdir(dirName) if fn.endswith('.bam')]

  for filePath in filePaths:
    extractIsolatedContacts(filePath)
