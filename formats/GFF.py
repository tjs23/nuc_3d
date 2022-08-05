from formats.Util import getFileObj
from numpy import array, uint8, float32, uint32, loadtxt
from collections import defaultdict

def getFeatures(filePath):
  
  #from time import time
  
  #t0 =  time()
  
  #data = loadtxt(filePath, 'S32', comments='#', delimiter='\t', converters=None, skiprows=0, usecols=(2,), unpack=False).tolist()

  features = defaultdict(int)
  
  for line in getFileObj(filePath):
    if line[0] == '#':
      continue
    
    features[line.split('\t')[2]] += 1
  
  #print "Time taken: %.3f" % (time()-t0)
    
  return features
  

def importDataTrack(filePath, binSize=None, feature=None, onlyChromosomes=True):
  # Should work with GFF and GTF
  # No binning
  
  from util.Genome import getChromosomeNames
    
  dataDict = {}
  name = None

  fileObj = getFileObj(filePath)
  
  sep1 = ';' # v3
  sep2 = '='
  
  for line in fileObj:
    if line[0] == '#':
      continue
    data = line[:-1].split('\t')
    
    if len(data) > 8:
      attribs = data[8]
      
      if ';' in attribs:
        if '; ' in attribs: 
          sep1 = '; ' # v2
 
        if '=' not in attribs.split(';')[0]:
          sep2 = None
      
      else:
        sep1 = None
      
      break
      
  fileObj.seek(0)
  
  for line in fileObj:
    if line[0] == '#':
      continue
    
    data = line[:-1].split('\t')
    n = len(data)
    
    if n < 8:
      continue
      
    if feature and (data[2] != feature):
      continue
      
    chromo, source, feat, start, end, score, strand, frame = data[:8]
    
    if n > 8 and sep1:
      attribs = data[8].split(sep1)
      attribDict = dict([a.split(sep2) for a in attribs])
      label = attribDict.get('Name', feat)
      
    else:
      label = feat
        
    if chromo in dataDict:
      chromoData = dataDict[chromo]
    else:
      chromoData = [[], [], []]
      dataDict[chromo] = chromoData
    
    if score == '.':
      score = 1.0
      val = 1.0
    else:
      score = float(score)
      val = score/1000.0
    
    if strand == '-':
      start, end = end, start
    
    chromoData[0].append([start, end])
    chromoData[1].append([score, val])
    chromoData[2].append(label)
  
  chromoIds = list(dataDict)
  nameDict = getChromosomeNames(chromoIds)
  
  for chromo in chromoIds:
    chrName = nameDict[chromo]
    
    if chrName.startswith('NW_'):
      del dataDict[chromo]
      continue
  
    if chrName.startswith('NT_'):
      del dataDict[chromo]
      continue
   
    regionData, valueData, annotations = dataDict[chromo]
    dataDict[chromo][0] = array(regionData, uint32)
    dataDict[chromo][1] = array(valueData, float32)
    
  fileObj.close()

  return dataDict, name


