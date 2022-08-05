from formats.Util import binSparseContacts, binTrackData
from numpy import load, savez, array, uint32, float32, uint8


def exporInteractions(filePath, dataDict, name):
    
  kwArgs = {}
  
  for chromoPair in dataDict:
    chromoA, chromoB = chromoPair
    regionData, valueData, strands, annotations = dataDict[chromoPair]
    
    rKey = 'interactions/regions/%s/%s/%s' % (name, chromoA, chromoB)
    vKey = 'interactions/values/%s/%s/%s' % (name, chromoA, chromoB)
    
    kwArgs[rKey] = array(regionData, uint32)
    kwArgs[vKey] = array(valueData, float32)  
      
    if annotations is not None:
      aKey = 'interactions/annotations/%s/%s/%s' % (name, chromoA, chromoB)
      kwArgs[aKey] = array(annotations)
    
  savez(filePath, **kwArgs)
   

def importInteractions(filePath):

  inDataDict = load(filePath)

  dataDict = {}
  
  for key in inDataDict:
    if not key.startswith('interactions/'):
      msg = 'File %s does not contain Nuc interactions data'
      raise Exception(msg % filePath)
   
    null, typ, name, chromoA, chromoB = key.split('/')      
    chromoPair = tuple([chromoA, chromoB])
    
    rKey = 'interactions/regions/%s/%s/%s' % (name, chromoA, chromoB)
    vKey = 'interactions/values/%s/%s/%s' % (name, chromoA, chromoB)
    aKey = 'interactions/annotations/%s/%s/%s' % (name, chromoA, chromoB)
    
    regionData = array(inDataDict[rKey], uint32)
    valueData = array(inDataDict[vKey], float32) 

    if aKey in dataDict:
      annotations = inDataDict[aKey].tolist()
    else:
      annotations = None
 
    dataDict[chromoPair] = regionData, valueData, annotations

  return dataDict, name  
  
  
def exportDataTrack(filePath, dataDict, name, stranded=True):

  kwArgs = {}
  
  for chromo in dataDict:
    regionData, valueData, strands, annotations = dataDict[chromo]
    
    rKey = 'dataTrack/regions/%s/%s' % (name, chromo)
    vKey = 'dataTrack/values/%s/%s' % (name, chromo)
    
    kwArgs[rKey] = array(regionData, uint32)
    kwArgs[vKey] = array(valueData, float32)

    if annotations is not None:
      aKey = 'dataTrack/annotations/%s/%s' % (name, chromo)
      kwArgs[aKey] = array(annotations)
    
  savez(filePath, **kwArgs)


def importDataTrack(filePath, binSize=None):

  inDataDict = load(filePath)
  dataDict = {}
    
  chromos = set()
  for key in inDataDict:
    if not key.startswith('dataTrack/'):
      msg = 'File %s does not contain a Nuc data track'
      raise Exception(msg % filePath)

    null, typ, name, chromo = key.split('/')
    chromos.add(chromo)

  for chromo in sorted(chromos):
    rKey = 'dataTrack/regions/%s/%s' % (name, chromo)
    vKey = 'dataTrack/values/%s/%s' % (name, chromo)
    aKey = 'dataTrack/annotations/%s/%s' % (name, chromo)
    
    regionData = array(inDataDict[rKey], uint32)
    valueData = array(inDataDict[vKey], float32)
  

    if binSize:
      regionData, valueData = binTrackData(regionData, valueData, binSize)
      annotations = None 
    elif aKey in dataDict:
      annotations = inDataDict[aKey].tolist()
    else:
      annotations = None
 
    dataDict[chromo] = regionData, valueData, annotations

  return dataDict, name


def exportCoords(filePath, posDict, coordsDict):
  
  kwArgs = {}
  
  for chromo in posDict:
    pKey = 'positions/%s' % chromo
    cKey = 'coords/%s' % chromo
    
    kwArgs[pKey] = array(posDict[chromo], uint32)
    kwArgs[cKey] = array(coordsDict[chromo], float32)
    
  savez(filePath, **kwArgs)
   

def importCoords(filePath):

  dataDict = load(filePath)
  
  chromosomes = [cKey[7:] for cKey in dataDict if cKey[:7] == 'coords/']
  
  if not chromosomes:
    msg = 'File %s does not contain Nuc structure data'
    raise Exception(msg % filePath)
  
  posDict = {}
  coordsDict = {}
  
  for chromo in chromosomes:
    pKey = 'positions/%s' % chromo
    cKey = 'coords/%s' % chromo
    
    posDict[chromo] = array(dataDict[pKey], uint32)
    coordsDict[chromo] = array(dataDict[cKey], float32)
    
  return posDict, coordsDict


def exportContacts(filePath, contactDict):
    
  kwArgs = {}
  
  for chromoPair in contactDict:
    
    key = 'contacts/%s/%s' % chromoPair
    
    kwArgs[key] = array(contactDict[chromoPair], uint32)
    
  savez(filePath, **kwArgs)
   

def importContacts(filePath, binSize=None):

  dataDict = load(filePath)

  contactDict = {}
  
  for key in dataDict:
    if not key.startswith('contacts/'):
      msg = 'File %s does not contain Nuc contact data'
      raise Exception(msg % filePath)
      
    chromoPair = tuple(key.split('/')[1:])
    data = array(dataDict[key], uint32)
      
    if binSize and binSize > 1:
      data = binSparseContacts(data, binSize)
   
    if data is None:
      continue
      
    contactDict[chromoPair] = data

  return contactDict
  
  
def exportContactMatrix(filePath, matrixDict, startDict, binSize):
    
  kwArgs = {}
  
  for chromoPair in matrixDict:
    chrA, chrB = chromoPair
    startA, startB = startDict[chromoPair]
    key = 'contactMatrix/%s/%s/%d/%d/%d' % (chrA, chrB, startA, startB, binSize)
    
    kwArgs[key] = array(matrixDict[chromoPair], uint32)
    
  savez(filePath, **kwArgs)
   

def importContactMatrix(filePath):

  dataDict = load(filePath)
  
  startDict = {}
  matrixDict = {}
  
  for key in dataDict:
    if not key.startswith('contactMatrix/'):
      msg = 'File %s does not contain Nuc contact matrix data'
      raise Exception(msg % filePath)
    
    head, chrA, chrB, startA, startB, binSize = key.split('/')
      
    chromoPair = (chrA, chrB)
    matrixDict[chromoPair] = array(dataDict[key], uint32)
    startDict[chromoPair] = (int(startA), int(startB))
    
  binSize = int(binSize)
  
  return matrixDict, startDict, binSize
