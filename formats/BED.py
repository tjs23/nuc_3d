from formats.Util import getFileObj
from numpy import array, uint32, float32, ones

def exportDataTrack(filePath, dataDict, name, stranded):

  fileObj = getFileObj(filePath, 'w')
  write = fileObj.write

  template = 'chr%s\t%d\t%d\t%s\t%d\t%s\n' # chr, start, end, label, score, strand
  
  vMax = 0
  for chromo in dataDict:
    vMax = max(vMax, dataDict[chromo][1][0].max())
  
  if vMax <= 1.0:
    scale = 1000
  else:
    scale = 1  
    
  for chromo in dataDict:
    regionData, valueData, strands, annotations = dataDict[chromo]
    
    for i, region in enumerate(regionData):
      start, end = region
      
      if start > end:
        strand = '-'
        start, end = end, start
      
      else:
        strand = '+'
 
      if annotations is not None:
        label = annotations[i]
      
      else:
        label = '%d' % i
 
      origValue, value = valueData[i]
      score = int(value * scale)
  
      line = template % (chromo, start, end, label, score, strand)
      write(line)

  fileObj.close()
  

def importDataTrack(filePath, binSize=None):
  
  name = None  
  fileObj = getFileObj(filePath, 'r')

  line = fileObj.readline()
  fileObj.seek(0)
  nFields = len(line.split())
  haveAnno = nFields > 3
  haveVal = nFields > 4
  haveStrand = nFields > 5
  
  dataDict = {}
  
  max_value = 0
  
  if binSize:
  
    for line in fileObj:
      if line[0] == '#':
        continue
    
      data = line.split()
      chromo = data[0]
      start = int(data[1])
      end = int(data[2])
      
      if haveAnno:
        label = data[3]
      else:
        label = ''
      
      if haveVal:
        val = float(data[4])
      else:
        val = 1.0
       
      if start > end:
        start, end = end, start

      if chromo in dataDict:
        chromoData = dataDict[chromo]
      else:
        chromoData = {}
        dataDict[chromo] = {}
      
      keyA = int(start/binSize)
      keyB = int(end/binSize)
      
      if keyA == keyB:
        if keyA in chromoData:
           chromoData[keyA][0] += val
           
           if label:
             prev = chromoData[keyA][1]
             
             if prev and label != prev:
               chromoData[keyA][1] = ''
             else:
               chromoData[keyA][1] = label
         
        else:
          chromoData[keyA] = [val, label]

      else:
        delta = float(end - start)

        for key in range(keyA, keyB+1):
          a = max(start, key*binSize)
          b = min(end, (key+1)*binSize)
          f = (b-a)/delta

          if key in chromoData:
            chromoData[key][0] += f * val
            
            if label:
              prev = chromoData[key][1]
              
              if prev and label != prev:
                chromoData[key][1] = ''
              else:
                chromoData[key][1] = label
            
          else:
            chromoData[key] = [f * val, label]
     
    for chromo in dataDict:
      chromoData = dataDict[chromo]
      keys = sorted(chromoData.keys())
      values = [chromoData[k][0] for k in keys]
      
      if not values:
        continue
      
      starts = array(keys) * binSize
      ends = starts + (binSize-1)
      
      regionData = array([starts, ends], uint32).T
      valueData = array([values, values], float32).T
      max_value = max(max_value, max(values))
      
      if haveAnno:
        annotations = [chromoData[k][1] for k in keys]
        aSet = set(annotations)
        
        if len(aSet) == 1:
          if aSet == set(['']):
            annotations = None
          elif aSet == set(['-']):
            annotations = None
          elif aSet == set(['.']):
            annotations = None
      
      else:
        annotations = None
            
      dataDict[chromo] = [regionData, valueData, annotations]
 
  else:
    
    prev = None
    
    for line in fileObj:
      if line[0] == '#':
        continue
        
      data = line.split()
      chromo = data[0]
      start = int(data[1])
      end = int(data[2])
 
      if haveVal:
        score = float(data[4])
      else:
        score = 1.0
 
      if chromo in dataDict:
        chromoData = dataDict[chromo]
      else:
        prev = None
        chromoData = [[], [], []]
        dataDict[chromo] = chromoData
 
      if haveStrand:
        strand = data[5]
        if strand == '-':
          if start < end:
            start, end = end, start
 
        elif strand == '+':
          if start > end:
            start, end = end, start
      
      #if end - start < 1e4:
      #  continue

      #if abs(end - start) > 10e6:
      #  continue
      
      #if prev and start < prev:
      #  start = prev + 1
      
      chromoData[0].append([start, end])
      chromoData[1].append([score, score])
      max_value = max(max_value, score)
      prev = end
      
      if haveAnno:
        chromoData[2].append(data[2])
 
    for chromo in dataDict:
      regionData, valueData, annotations = dataDict[chromo]
 
      dataDict[chromo][0] = array(regionData, uint32)
      dataDict[chromo][1] = array(valueData, float32)
      
      if annotations:
        aSet = set(annotations)
        if len(aSet) == 1:
          if aSet == set(['']):
            dataDict[chromo][2] = None
          elif aSet == set(['-']):
            dataDict[chromo][2] = None
          elif aSet == set(['.']):
            dataDict[chromo][2] = None
 
      else:
        dataDict[chromo][2] = None
  
  if max_value == 0:
    for chromo in dataDict:
      dataDict[chromo][1] = ones(dataDict[chromo][1].shape, float32)
   
  elif max_value > 100:
    for chromo in dataDict:
      dataDict[chromo][1] = dataDict[chromo][1]/100
  
  elif max_value > 10:
    for chromo in dataDict:
      dataDict[chromo][1] = dataDict[chromo][1]/10
  
    
  fileObj.close()

  return dataDict, name

