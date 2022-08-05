import json
from numpy import array, uint32, float32
from formats.Util import binSparseContacts, binTrackData
from formats.Util import getFileObj

def exportInteractions(filePath, dataDict, name):

  subDict = {}
  dataDict = {'interactions': subDict, 'trackName': name}
  
  for chromoPair in sorted(dataDict):
    regionData, valueData, annotations = dataDict[chromoPair]
    key = '/'.join(chromoPair)
    
    pairDict = {}
    subDict[key] = pairDict
    
    pairDict['regions'] = array(regionData, int).tolist()
    pairDict['values'] = array(valueData, float).tolist()
    
    if annotations is not None:
      pairDict['annotations'] = annotations.tolist()
 
  fileObj = getFileObj(filePath, 'w')
  
  json.dump(dataDict, fileObj)
  
  fileObj.close()
  

def importInteractions(filePath):

  fileObj = getFileObj(filePath, 'r')
  
  inDataDict = json.load(fileObj)
    
  if 'interactions' not in inDataDict:
    msg = 'File %s does not contain a Nuc interactions set'
    raise Exception(msg % filePath)
  
  dataDict = {}
  subDict = inDataDict['interactions']
  name = inDataDict['trackName']
  
  for chromoKey in subDict:
    chromoPair = tuple(chromoKey.split('/'))
    regionData = array(subDict[chromoKey]['regions'], uint32)
    valueData = array(subDict[chromoKey]['values'], float32)
    
    if 'annotations' in subDict:
      annotations = subDict[chromoKey]['annotations']
    else:
      annotations = None 
    
    dataDict[chromoPair] = regionData, valueData, annotations

  fileObj.close()

  return dataDict, name


def exportDataTrack(filePath, dataDict, name, stranded):

  subDict = {}
  dataDict = {'dataTrack': subDict, 'trackName': name}
  
  for chromo in sorted(dataDict):
    regionData, valueData, annotations = dataDict[chromo]
    chromoDict = {}
    subDict[chromo] = chromoDict
    
    chromoDict['regions'] = array(regionData, int).tolist()
    chromoDict['values'] = array(valueData, float).tolist()
    
    if annotations is not None:
      chromoDict['annotations'] = annotations.tolist()
 
  fileObj = getFileObj(filePath, 'w')
  
  json.dump(dataDict, fileObj)
  
  fileObj.close()


def importDataTrack(filePath, binSize=None):

  fileObj = getFileObj(filePath, 'r')
  
  inDataDict = json.load(fileObj)
    
  if 'dataTrack' not in inDataDict:
    msg = 'File %s does not contain a Nuc data track'
    raise Exception(msg % filePath)
  
  dataDict = {}
  subDict = inDataDict['dataTrack']
  name = inDataDict['trackName']
  
  for chromo in subDict:
    regionData = array(subDict['regions'], uint32)
    valueData = array(subDict['values'], float32)
    
    
    if binSize:
      regionData, valueData = binTrackData(regionData, valueData, binSize)
      annotations = None 
    elif 'annotations' in subDict:
      annotations = subDict['annotations']
    else:
      annotations = None 
    
    dataDict[chromo] = regionData, valueData, annotations

  fileObj.close()

  return dataDict, name


def exportCoords(filePath, posDict, coordsDict):

  subDict = {}
  dataDict = {'structure': subDict}
  
  for chromo in posDict:
    subDict[chromo] = {'positions':array(posDict[chromo], int).tolist(),
                        'coords':array(coordsDict[chromo], float).tolist()}
  
  fileObj = getFileObj(filePath, 'w')
  
  json.dump(dataDict, fileObj)
  
  fileObj.close()
  

def importCoords(filePath):

  fileObj = getFileObj(filePath, 'r')
  
  dataDict = json.load(fileObj)
  
  if 'structure' not in dataDict:
    msg = 'File %s does not contain Nuc structure data'
    raise Exception(msg % filePath)
  
  subDict = dataDict['structure']
  posDict = {}
  coordsDict = {}

  for chromo in subDict:
      
    posDict[chromo] = array(subDict[chromo]['positions'], uint32)
    coordsDict[chromo] = array(subDict[chromo]['coords'], float32)

  fileObj.close()
  
  return posDict, coordsDict


def exportContacts(filePath, contactDict):
    
  fileObj = getFileObj(filePath, 'w')
  
  subDict = {}
  dataDict = {'contacts': subDict}
  
  for chromoPair in contactDict:
    key = '/'.join(chromoPair)
    subDict[key] = array(contactDict[chromoPair], int).tolist()
  
  json.dump(dataDict, fileObj)
  
  fileObj.close() 


def importContacts(filePath, binSize=None):

  fileObj = getFileObj(filePath, 'r')
  
  dataDict = json.load(fileObj)
  
  if 'contacts' not in dataDict:
    msg = 'File %s does not contain Nuc contact data'
    raise Exception(msg % filePath)
  
  subDict = dataDict['contacts']  
  contactDict = {}
    
  for key in subDict:
    data = array(subDict[key], uint32)
   
    if binSize and binSize > 1:
      data = binSparseContacts(data, binSize)
      
    if data is None:
      continue
  
    chromoPair = tuple(key.split('/'))
    contactDict[chromoPair] = data
  
  fileObj.close()
  
  return contactDict


def exportContactMatrix(filePath, matrixDict, startDict, binSize):
    
  fileObj = getFileObj(filePath, 'w')
  
  subDict = {}
  dataDict = {'contactMatrix': subDict, 'binSize': binSize}
  
  for chromoPair in matrixDict:
    key = '/'.join(chromoPair)
    data = array(matrixDict[chromoPair], int).tolist()
    seqStarts = startDict[chromoPair],
    subDict[key] = {'data':data, 'seqStarts': seqStarts}
  
  json.dump(dataDict, fileObj)
  
  fileObj.close() 


def importContactMatrix(filePath):

  fileObj = getFileObj(filePath, 'r')
  
  dataDict = json.load(fileObj)
  
  if 'contactMatrix' not in dataDict:
    msg = 'File %s does not contain Nuc contact matrix data'
    raise Exception(msg % filePath)
  
  binSize = dataDict['binSize']
  subDict = dataDict['contactMatrix']
  startDict = {}
  matrixDict = {}
  
  for key in subDict:
    chromoPair = tuple(key.split('/'))
    data = matrixDict[chromoPair]['data']
    seqStarts = matrixDict[chromoPair]['seqStarts']
    
    matrixDict[chromoPair] = array(data, uint32)
    startDict[chromoPair] = tuple(seqStarts)
  
  fileObj.close()
  
  return matrixDict, startDict, binSize
