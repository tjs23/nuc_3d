from formats.Util import getFileObj
from numpy import array, uint32, float32, empty, uint8

def exportInteractions(filePath, dataDict, name, separator='\t'):

  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write
  
  heads = ('#chrA', 'startA', 'endA', 'chrB', 'startB', 'endB', 'origValue','value','label')
  template = separator.join(['%s','%d','%d','%s','%d','%d','%.6f','%.6f','%s\n'])
  line = separator.join(heads)
  write(line + '\n')
  
  for chromoPair in dataDict:
    chrA, chrB = chromoPair
    regionData, valueData, annotations = dataDict[chromoPair]
 
    for i, region in enumerate(regionData):
      a1, a2, b1, b2 = region
        
      if annotations is not None:
        annotation = annotations[i]
      else:
        annotation = ''
 
      origValue, value = valueData[i]
 
      line = template % (chrA, a1, a2, chrB, b1, b2, origValue, value, annotation)
      write(line)

  fileObj.close()


def importInteractions(filePath, separator='\t'):

  fileObj = getFileObj(filePath, 'r')
  
  dataDict = {}
  
  line = fileObj.readline()

  if line.startswith('#'):
    if ':' in line:
      name, heads = line[1:].split(':')
    else:
      name = None
    
  else:
    name = None
    fileObj.seek(0)
  
  for line in fileObj:
    chrA, a1, a2, chrB, b1, b2, origVal, val, label = line[:-1].split(separator)
    
    a1 = int(a1)
    a2 = int(a2)
    b1 = int(b2)
    b2 = int(b2)
    origVal = float(origVal)
    val = float(val)
    chromoPair = (chrA, chrB)
    
    if chromoPair in dataDict:
      chromoData = dataDict[chromoPair]
    else:
      chromoData = [[], [], []]
      dataDict[chromoPair] = chromoData
 
    chromoData[0].append([a1, a2, b1, b2])
    chromoData[1].append([origVal, val])
    chromoData[2].append(label)
 
  for chromoPair in dataDict:
    regionData, valueData, annotations = dataDict[chromoPair]
 
    dataDict[chromoPair][0] = array(regionData, uint32)
    dataDict[chromoPair][1] = array(valueData, float32)

    aSet = set(annotations)
    if len(aSet) == 1:
      if aSet == set(['']):
        dataDict[chromoPair][2] = None
      elif aSet == set(['-']):
        dataDict[chromoPair][2] = None
      elif aSet == set(['.']):
        dataDict[chromoPair][2] = None
 
  fileObj.close()
  
  return dataDict, name  
  
  
def exportDataTrack(filePath, dataDict, name, stranded, separator='\t'):

  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write
  
  heads = ('#%s:chr' % name,'start','end','origValue','value','label')
  template = separator.join(['%s','%d','%d','%.6f','%.6f','%s\n'])
  line = separator.join(heads)
  write(line + '\n')
  
  for chromo in dataDict:
    regionData, valueData, annotations = dataDict[chromo]
 
    for i, region in enumerate(regionData):
      start, end = region
      
      if not stranded:
        if start > end:
          start, end = end, start
        
      if annotations is not None:
        annotation = annotations[i]
      else:
        annotation = ''
 
      origValue, value = valueData[i]
 
      line = template % (chromo, start, end, origValue, value, annotation)
      write(line)

  fileObj.close()


def importDataTrack(filePath, binSize=None):
    
  fileObj = getFileObj(filePath, 'r')
  
  dataDict = {}
  
  line = fileObj.readline()
  
  if line.startswith('#'):
    name, heads = line[1:].split(':')
  
  else:
    name = None
    fileObj.seek(0)
  
  if binSize:
    for line in fileObj:
      chromo, start, end, origVal, val, label = line.split()
      start = int(start)
      end = int(end)
      
      if start > end:
        start, end = end, start
      
      origVal = float(origVal)
      val = float(val)
 
      if chromo in dataDict:
        chromoData = dataDict[chromo]
      else:
        chromoData = {}
        dataDict[chromo] = {}
      
      keyA = int(start/binSize)
      keyB = int(end/binSize)
      
      if keyA == keyB:
        if keyA in chromoData:
          chromoData[keyA][0] += origVal
          chromoData[keyA][1] += val
          
          if label:
            prev = chromoData[keyA][2]
            
            if prev and label != prev:
              chromoData[keyA][2] = ''
            else:
              chromoData[keyA][2] = label
        
        else:
          chromoData[keyA] = [origVal, val, label]

      else:
        delta = float(end - start)

        for key in range(keyA, keyB+1):
          a = max(start, key*binSize)
          b = min(end, (key+1)*binSize)
          f = (b-a)/delta

          if key in chromoData:
            chromoData[key][0] += f * origVal
            chromoData[key][1] += f * val
            
            if label:
              prev = chromoData[key][2]
              
              if prev and label != prev:
                chromoData[key][2] = ''
              else:
                chromoData[key][2] = label
            
          else:
            chromoData[key] = [f * origVal, f * val, label]
     
    for chromo in dataDict:
      chromoData = dataDict[chromo]
      keys = sorted(chromoData.keys())
      values = [chromoData[k][:2] for k in keys]
      annotations = [chromoData[k][2] for k in keys]
      
      starts = array(keys) * binSize
      ends = starts + (binSize-1)
      
      regionData = array([starts, ends], uint32).T
      valueData = array(values, float32)
      
      aSet = set(annotations)
      if len(aSet) == 1:
        if aSet == set(['']):
          annotations = None
        elif aSet == set(['-']):
          annotations = None
        elif aSet == set(['.']):
          annotations = None
      
      dataDict[chromo] = regionData, valueData, annotations
      
  else:
    for line in fileObj:
      chromo, start, end, origVal, val, label = line.split()
      start = int(start)
      end = int(end)
      origVal = float(origVal)
      val = float(val)
 
      if chromo in dataDict:
        chromoData = dataDict[chromo]
      else:
        chromoData = [[], [], []]
        dataDict[chromo] = chromoData
 
      chromoData[0].append([start, end])
      chromoData[1].append([origVal, val])
      chromoData[2].append(label)
 
    for chromo in dataDict:
      regionData, valueData, annotations = dataDict[chromo]
 
      dataDict[chromo][0] = array(regionData, uint32)
      dataDict[chromo][1] = array(valueData, float32)

      aSet = set(annotations)
      if len(aSet) == 1:
        if aSet == set(['']):
          dataDict[chromo][2] = None
        elif aSet == set(['-']):
          dataDict[chromo][2] = None
        elif aSet == set(['.']):
          dataDict[chromo][2] = None
  
  fileObj.close()

  return dataDict, name


def exportCoords(filePath, posDict, coordsDict):
  
  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write

  for chromo in posDict:
    chromoCoords = coordsDict[chromo]
    chromoPos = posDict[chromo]
    
    nModels = len(chromoCoords)
    nCoords = len(chromoPos)
    
    line = '%s\t%d\t%d\n' % (chromo, nCoords, nModels)
    write(line)
    
    for j in range(nCoords):
      data = chromoCoords[:,j].ravel().tolist()
      data = '\t'.join('%.8f' % d for d in  data)
      
      line = '%d\t%s\n' % (chromoPos[j], data)
      write(line)

  fileObj.close()
  

def importCoords(filePath):
    
  fileObj = getFileObj(filePath, 'r')
  
  posDict = {}
  coordsDict = {}  
  chromo = None
  
  for line in fileObj:
    
    data = line.split()
    nItems = len(data)
    
    if not nItems:
      continue
    
    elif data[0] == '#':
      continue
    
    elif nItems == 3:
      chromo, nCoords, nModels = data
      nCoords = int(nCoords)
      nModels = int(nModels)
      
      chromoPos = []
      chromoCoords = empty((nModels, nCoords, 3), float32)
      
      coordsDict[chromo] = chromoCoords
      posDict[chromo] = chromoPos
      
      nCoords = int(nCoords)
      nModels = int(nModels)
      check = (nModels * 3) + 1
      i = 0
      
    elif not chromo:
      raise Exception('Missing chromosome record in file %s' % filePath)
     
    elif nItems != check:
      msg = 'Data size in file %s does not match Position + Models * Positions * 3'
      raise Exception(msg % filePath)
    
    else:
      chromoPos.append(int(data[0]))
      
      coord = [float(x) for x in data[1:]]
      coord = array(coord).reshape(nModels, 3)
      chromoCoords[:,i] = coord
      i += 1
  
  fileObj.close()
  
  return posDict, coordsDict


def exportContactMatrix(filePath, matrixDict, startDict, binSize):
  
  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write
 
  for chromoPair in matrixDict:
    chrA, chrB = chromoPair
    matrix = matrixDict[chromoPair]
    startA, startB = startDict[chromoPair]
    n, m = matrix.shape
   
    colHeads = [str(startB+i*binSize) for i in xrange(m)]
    line = '%s|%s\t%s\n' % (chrA, chrB, '\t'.join(colHeads))
    write(line)
    
    for i in range(n):
      items = [str(matrix[i,j]) for j in xrange(m)]
      line = '%d\t%s\n' % (startA+i*binSize, '\t'.join(items))
      write(line)
 
  fileObj.close()


def importContactMatrix(filePath):
    
  fileObj = getFileObj(filePath, 'r')
  
  matrix = []
  matrixDict = {}
  startDict = {}
  chromoPair = None
  startA = None
  startB = None
  
  for line in fileObj:
    data = line.split()
    nItems = len(data)
    
    if not nItems:
      continue
    
    elif data[0] == '#':
      continue
    
    elif "|" in data[0]:
      if matrix and chromoPair:
        matrixDict[chromoPair] = array(matrix, int)
        startDict[chromoPair] = (startA, startB)
      
      chrA, chrB = data[0].split('|')
      chromoPair = (chrA, chrB)
      binSize = int(data[2]) - int(data[1])
      startB = int(data[1])
      startA = None
      m = len(data)-1
      matrix = []
      
    elif not chromoPair:
      raise Exception('Cannot read contact matrix file %s' %filePath)
    
    else:
      if startA is None:
        startA = int(data[0])
      
      vals = [int(x) for x in data[1:]]
      
      if len(vals) != m:
        raise Exception('Inconsistent number of columns in contact matrix file %s' %filePath)
      
      matrix.append(vals) 
      
  if matrix and chromoPair:
    matrixDict[chromoPair] = array(matrix, int)

  fileObj.close()
  
  return matrixDict, startDict, binSize
