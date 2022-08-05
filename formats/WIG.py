import re
from formats.Util import getFileObj
from numpy import array, uint32, float32

PARAM_PATT = re.compile('\S*\s*(\S+)=(.+?)(\s+\S+=.+|\n|$)')

def exportDataTrack(filePath, dataDict, name, stranded):

  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write
  template = 'variableStep chrom=chr%s span=%d\n%d %.3f\n'
        
  for chromo in dataDict:
    regionData, valueData, strands, annotations = dataDict[chromo]
    
    for i, region in enumerate(regionData):
     start, end = sorted(region)
     origValue, value = valueData[i]
     delta = end-start
     line = template % (chromo, delta, start, origValue)
 
     write(line)

  fileObj.close()


def _getParamDict(line):
  
  paramDict = {}
  matchObj = PARAM_PATT.match(line)
  
  while matchObj:
    key = matchObj.group(1)
    val = matchObj.group(2)
    line = matchObj.group(3)
    
    if val[0] in '\'"':
      val = val[1:-1]
    
    paramDict[key] = val
    matchObj = PARAM_PATT.match(line)
  
  return paramDict
  

def importDataTrack(filePath, binSize=None):
  
  fileObj = getFileObj(filePath)
  isFixed = False
  name = None
  regionDict = {}
  valueDict = {}
  span = 1
  step = 1
  
  for line in fileObj:
  
    if line.startswith('track'):
      paramDict = _getParamDict(line)
      name = paramDict.get('name')
      
    elif line.startswith('variableStep'):
      isFixed = False
      paramDict = _getParamDict(line)
      chromo = paramDict['chrom']
      
      if chromo not in regionDict:
        if binSize:
          regionDict[chromo] = {}
         
        else:
          regionDict[chromo] = []
          valueDict[chromo] = []
      
      span = int(paramDict.get('span', 1))
      
    elif line.startswith('fixedStep'):
      isFixed = True
      paramDict = _getParamDict(line)
      chromo = paramDict['chrom']

      if chromo not in regionDict:
        if binSize:
          regionDict[chromo] = {}
         
        else:
          regionDict[chromo] = []
          valueDict[chromo] = []
     
      pos = int(paramDict['start'])
      step = int(paramDict['step'])
      span = int(paramDict.get('span', 1))
    
    else:
      data = line.split()
       
      if isFixed:
        if len(data) != 1:
          continue
 
        val = float(data[0])
        
        if binSize:
          posA = pos
          posB = pos + span
          keyA = posA / binSize
          keyB = posB / binSize
          
          if keyA == keyB:
            if keyA in regionDict[chromo]:
              regionDict[chromo][keyA] += val
            else:
              regionDict[chromo][keyA] = val
 
          else:
            delta = float(posB - posA)
 
            for key in range(keyA, keyB+1):
              a = max(posA, key*binSize)
              b = min(posB, (key+1)*binSize)
              f = val * (b-a)/delta
 
              if key in regionDict[chromo]:
                regionDict[chromo][key] += val
              else:
                regionDict[chromo][key] = val
        
        else:
          regionDict[chromo].append( [pos, pos+span] )
          valueDict[chromo].append( [val, val])
 
        pos += step
 
      else:
        if len(data) != 2:
          continue
 
        pos, val = data
        pos = int(pos)
        val = float(val)
 
        if binSize:
          posA = pos
          posB = pos + span
          keyA = int(posA / binSize)
          keyB = int(posB / binSize)
 
          if keyA == keyB:
            if keyA in regionDict[chromo]:
              regionDict[chromo][keyA] += val
            else:
              regionDict[chromo][keyA] = val
 
          else:
            delta = float(posB - posA)
 
            for key in range(keyA, keyB+1):
              a = max(posA, key*binSize)
              b = min(posB, (key+1)*binSize)
              f = val * (b-a)/delta
 
              if key in regionDict[chromo]:
                regionDict[chromo][key] += f
              else:
                regionDict[chromo][key] = f
        
        else:
          regionDict[chromo].append( [pos, pos+span] )
          valueDict[chromo].append( [val, val])
  
  chromos = list(regionDict.keys())
  dataDict = {}
  if binSize:
    for chromo in chromos:
      positions = []
      values = []
 
      for key in regionDict[chromo]:
        pos = key * binSize
        val = regionDict[chromo][key]
        positions.append([pos, pos+binSize])
        values.append([val, val])
 
      dataDict[chromo] = [array(positions, uint32),
                          array(values, float32),
                          None]
      del regionDict[chromo]
  
  else:
    for chromo in chromos:
      dataDict[chromo] = [array(regionDict[chromo], uint32),
                          array(valueDict[chromo], float32),
                          None]
      del regionDict[chromo]
      del valueDict[chromo]
  

  fileObj.close()

  return dataDict, name


