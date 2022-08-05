from formats.Util import getFileObj
from numpy import array, float32, uint32

def exportDataTrack(filePath, dataDict, name, stranded):

  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write
  write('track type=bedGraph name="%s"\n' % name)
  template = 'chr%s\t%d\t%d\t%.5f\n' # chr, start, end, value

  for chromo in dataDict: 
    regionData, valueData, annotations = dataDict[chromo]
   
    for i, region in enumerate(regionData):
      start, end = sorted(region)
      origValue, value = valueData[i] 
      line = template % (chromo, start, end, origValue)
      write(line)

  fileObj.close()


def importDataTrack(filePath, binSize=None):
  
  name = None  
  fileObj = getFileObj(filePath, 'r')
  dataDict = {}
  
  if binSize:
    for line in fileObj:
      chromo, start, end, val = line.split()
      start = int(start)
      end = int(end)
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
          chromoData[keyA] += val
         
        else:
          chromoData[keyA] = val

      else:
        delta = float(end - start)

        for key in range(keyA, keyB+1):
          a = max(start, key*binSize)
          b = min(end, (key+1)*binSize)
          f = (b-a)/delta

          if key in chromoData:
            chromoData[key] += f * val
            
          else:
            chromoData[key] = f * val
     
    for chromo in dataDict:
      chromoData = dataDict[chromo]
      keys = sorted(chromoData.keys())
      values = [chromoData[k] for k in keys]
      
      starts = array(keys) * binSize
      ends = starts + (binSize-1)
      
      regionData = array([starts, ends], uint32).T
      valueData = array([values, values], float32).T
            
      dataDict[chromo] = regionData, valueData, None
  
  else:
    for line in fileObj:
      chromo, start, end, val = line.split()
      start = int(start)
      end = int(end)
      val = float(val)
 
      if chromo in dataDict:
        chromoData = dataDict[chromo]
      else:
        chromoData = [[], [], None]
        dataDict[chromo] = chromoData
 
      chromoData[0].append([start, end])
      chromoData[1].append([val, val])
 
    for chromo in dataDict:
      regionData, valueData, annotations = dataDict[chromo]
 
      dataDict[chromo][0] = array(regionData, uint32)
      dataDict[chromo][1] = array(valueData, float32)
  
  fileObj.close()

  return dataDict, name


