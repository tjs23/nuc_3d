from formats.Util import getFileObj
from numpy import array, uint32, float32

def exportDataTrack(filePath, dataDict, name, stranded):

  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write
  template = 'chr%s\t%d\t%d\t%s\t%d\t%s\t%.5f\t-1\t-1\n' # pValue, qValue

            
  for chromo in dataDict:
    regionData, valueData, annotations = dataDict[chromo]
    
    for i, region in enumerate(regionData):
      start, end = region
       
      if stranded:
        if start > end:
          start, end = end, start
          strand = '-'
        else:
          strand = '+'
      
      else:
        if start > end:
          start, end = end, start
          
        strand = '.'
 
      if annotations is not None:
        label = annotations[i]
      
      else:
        label = '.'
 
      origVal, value = valueData[i]
      score = int(value * 1000)
      line = template % (chromo, start, end, label, score, strand, origVal) # , pValue, qValue)
      
      write(line)

  fileObj.close()


def importDataTrack(filePath, binSize=None):
  
  name = None  
  fileObj = getFileObj(filePath, 'r')
  
  dataDict = {}
  
  if binSize:
    for line in fileObj:
      chromo, start, end, label, val, strand, origVal = line.split()
      start = int(start)
      end = int(end)
      val = float(val)
      origVal = float(origVal)
 
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
      chromo, start, end, label, score, strand, origVal = line.split()
      start = int(start)
      end = int(end)
      score = float(score)
      origVal = float(origVal)
 
      if strand == '-':
        if start < end:
          start, end = end, start
 
      if chromo in dataDict:
        chromoData = dataDict[chromo]
      else:
        chromoData = [[], [], []]
        dataDict[chromo] = chromoData
 
      chromoData[0].append([start, end])
      chromoData[1].append([origVal, score])
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


