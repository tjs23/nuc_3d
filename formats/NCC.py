from formats.Util import getFileObj
from numpy import array, uint32

def exportContacts(filePath, contactDict):
  
  # extend this to support ambiguity
    
  fileObj = getFileObj(filePath, 'w')
  
  write = fileObj.write
  
  i = 1
  for chromoPair in contactDict:
    chrA, chrB = chromoPair
    
    if len(contactDict[chromoPair][0]) == 4: # Has ambiguity data
      for posA, posB, numObs, ambig in contactDict[chromoPair]:
        line = '%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n' % (chrA, posA, posA+1, 0, 0, '+', chrB, posB, posB+1, 0, 0, '+', ambig, i, 0)
        i += 1
        write(line)
   
    else:
      for posA, posB, numObs in contactDict[chromoPair]:
        line = '%s\t%d\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n' % (chrA, posA, posA+1, 0, 0, '+', chrB, posB, posB+1, 0, 0, '+', i, i, 0)
        i += 1
        write(line)
 
  
  fileObj.close() 


def importContacts(filePath, binSize=None):
  
  # TBD: Fast text file reader in Cython
  
  from collections import defaultdict
  
  fileObj = getFileObj(filePath, 'r')
  
  contactDict = {}

  if binSize and binSize > 1:
    binDict = {}
    
    for line in fileObj:
      chrA, f_start_a, f_end_a, start_a, end_a, strand_a, chrB, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_code, pair_id, swap_pair = line.split()
     
      if strand_a == '+':
        posA = int(f_start_a)
      else:
        posA = int(f_end_a)
      
      if strand_b == '+':
        posB = int(f_start_b)        
      else:
        posB = int(f_end_b)
           
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
    ambig = defaultdict(int)

    ambig_group = 0
    for line in fileObj:
      row = line.split()
      ambig_code = row[12]
           
      if '.' in ambig_code:
        count, selected = ambig_code.split('.')
        if count != '0':
          ambig_group += 1
       
        if selected == '0':
          continue
    
      else:
        ambig_group = int(ambig_code)
    
      ambig[ambig_group] += 1
         
    ambig_group = 0
    fileObj.seek(0)
    
    for line in fileObj:
      chrA, f_start_a, f_end_a, start_a, end_a, strand_a, chrB, f_start_b, f_end_b, start_b, end_b, strand_b, ambig_code, pair_id, swap_pair = line.split()
      
      if strand_a == '+':
        posA = int(f_start_a)
      else:
        posA = int(f_end_a)
      
      if strand_b == '+':
        posB = int(f_start_b)       
      else:
        posB = int(f_end_b)
      
      if '.' in ambig_code:
        count, selected = ambig_code.split('.')
        if count != '0':
          ambig_group += 1
       
        if selected == '0':
          continue
        
      else:
        ambig_group = int(ambig_code)
      
      if ambig[ambig_group] > 1:
        continue
       
      if chrA > chrB:
        chrA, chrB = chrB, chrA
        posA, posB = posB, posA
 
      chromoPair = (chrA, chrB)
 
      if chromoPair in contactDict:
        contactDict[chromoPair].append((posA, posB, 1, ambig_group))
      else:
        contactDict[chromoPair] = [(posA, posB, 1, ambig_group)]
     
  fileObj.close()
  
  for chromoPair in contactDict:
    contactDict[chromoPair] = array(contactDict[chromoPair], uint32)
  
  return contactDict

 
 
        
          
             

