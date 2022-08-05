from formats.Util import getFileObj
from numpy import array, uint32

def exportContacts(filePath, contactDict):
  
  # extend this to support ambiguity
    
  fileObj = getFileObj(filePath, 'w')
  
  write = fileObj.write
  write('#%s\t%s\t%s\t%s\n' % ('chr_A', 'pos_A', 'chr_B', 'pos_B'))
  #write('#%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('chr_A', 'start_A', 'end_A', 'chr_B', 'start_B', 'end_b', 'num_obs'))

  chromo_pairs = sorted(contactDict.keys())

  for chromoPair in contactDict:
    chrA, chrB = chromoPair
    
    for row in contactDict[chromoPair].T:
      posA, posB, numObs = row[:3]
      line = 'chr%s\t%d\tchr%s\t%d\n' % (chrA, posA, chrB, posB)
      #line = 'chr%s\t%d\t%d\tchr%s\t%d\t%d\t%d\n' % (chrA, posA, posA+49999, chrB, posB, posB+49999, numObs)
      write(line)
  
  fileObj.close() 


def importContacts(filePath, binSize=None):
  
  # TBD: Fast text file reader in Cython
  
  fileObj = getFileObj(filePath, 'r')
  
  contactDict = {}

  if binSize:
    binDict = {}
    
    for line in fileObj:
      data = line.split()
      nItems = len(data)
 
      if not nItems:
        continue
 
      elif data[0] == '#':
        continue

      elif nItems != 5:
        raise Exception('Cannot import contacts from file %s' % filePath)
 
      else:
        chrA, posA, chrB, posB, numObs = data
        posA = binSize * int(posA//binSize)
        posB = binSize * int(posB//binSize)
 
        key = tuple(sorted([(chrA, posA), (chrB, posB)]))
 
        if key in binDict:
          binDict[key] += 1
        else:
          binDict[key] = 1
  
    for key in binDict:
      pairA, pairB = key
      chrA, posA = pairA
      chrB, posB = pairB
      numObs = binDict[key]
      chromoPair = (chrA, chrB)

      if chromoPair in contactDict:
        contactDict[chromoPair].append( (posA, posB, numObs) )
      else:
        contactDict[chromoPair] = [(posA, posB, numObs)]
        
  else:
   
    for line in fileObj:
      data = line.split()
      nItems = len(data)
 
      if not nItems:
        continue
 
      elif line[0] == '#':
        continue

      elif nItems < 4:
        raise Exception('Cannot import contacts from file %s' % filePath)
 
      else:
        if nItems == 4:
          chrA, posA, chrB, posB = data
          numObs = 1
        else:
          chrA, posA, chrB, posB, numObs = data[:5]
          numObs = int(numObs)
          
        posA = int(posA)
        posB = int(posB)
 
 
        if chrA > chrB:
          chrA, chrB = chrB, chrA
          posA, posB = posB, posA
 
        chromoPair = (chrA, chrB)
 
        if chromoPair in contactDict:
          contactDict[chromoPair].append((posA, posB, numObs))
        else:
          contactDict[chromoPair] = [(posA, posB, numObs)]
     
  fileObj.close()
  
  for chromoPair in contactDict:
    contactDict[chromoPair] = array(contactDict[chromoPair], uint32)
  
  return contactDict

 
 
        
          
             

