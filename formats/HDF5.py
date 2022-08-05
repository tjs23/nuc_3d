from formats.Util import binSparseContacts, binTrackData

from numpy import array, uint32, uint8, float32, int32, string_
from h5py import File, Group, special_dtype

VL_str = special_dtype(vlen=str)

def exportInteractions(filePath, dataDict, name):

  root = File(filePath, mode='w')
  
  group = root.create_group('interactions')
  group.attrs['trackName'] = string_(name)
 
  for chromoPair in dataDict:
    chrA, chrB = chromoPair
    
    if chrA in group:
      chromoGroup = group[chrA]
    else:
      chromoGroup = group.create_group(chrA)
    
    if chrB in chromoGroup:
      chromoGroupB = chromoGroup[chrB]
    else:
      chromoGroupB = chromoGroup.create_group(chrB)
      
    regionData, valueData, annotations = dataDict[chromoPair]
    
    chromoGroupB.create_dataset('regions', dtype=uint32, data=regionData)
    chromoGroupB.create_dataset('values', dtype=float32, data=valueData)
    
    if annotations is not None:
      chromoGroupB.create_dataset('annotations', dtype=VL_str, data=annotations)

  root.close()


def importInteractions(filePath):

  root = File(filePath, mode='r')
  
  if 'interactions' not in root:
    msg = 'File %s does not contain Nuc interactions'
    raise Exception(msg % filePath)

  group = root['interactions']
  name = str(group.attrs['trackName'])
  dataDict = {}

  for chrA in group:
    chromoGroup = group[chrA]
 
    for chrB in chromoGroup:
      regionData = array(chromoGroup[chrB]['regions'], uint32)
      valueData = array(chromoGroup[chrB]['values'], float32)
      
      if 'annotations' in chromoGroup[chrB]:
        annotations = array(chromoGroup[chrB]['annotations']).tolist()
      else:
        annotations = None
 
      dataDict[(chrA, chrB)] = (regionData, valueData, annotations)
 
  return dataDict, name


def exportDataTrack(filePath, dataDict, name, stranded):

  root = File(filePath, mode='w')
  
  group = root.create_group('dataTrack')
  group.attrs['trackName'] = string_(name)

  for chromo in sorted(dataDict):
    regionData, valueData, annotations = dataDict[chromo]
    chromoGroup = group.create_group(chromo)
    
    chromoGroup.create_dataset('regions', dtype=uint32, data=regionData)
    chromoGroup.create_dataset('values', dtype=float32, data=valueData)
    
    if annotations is not None:
      chromoGroup.create_dataset('annotations', dtype=VL_str, data=annotations)

  root.close()


def importDataTrack(filePath, binSize=None):

  root = File(filePath, mode='r')
  
  if 'dataTrack' not in root:
    msg = 'File %s does not contain Nuc data tracks'
    raise Exception(msg % filePath)

  group = root['dataTrack']
  name = str(group.attrs['trackName'])
  dataDict = {}
  
  for chromo in group:
    chromoGroup = group[chromo]    
    regionData = array(chromoGroup['regions'], uint32)
    valueData = array(chromoGroup['values'], float32)
    
  
    if binSize:
      regionData, valueData = binTrackData(regionData, valueData, binSize)
      annotations = None 
    elif 'annotations' in chromoGroup:
      annotations = array(chromoGroup['annotations']).tolist()
    else:
      annotations = None 
    
    dataDict[chromo] = regionData, valueData, annotations

  return dataDict, name
  

def exportCoords(filePath, posDict, coordsDict):

  root = File(filePath, mode='w')
  
  group = root.create_group('structure')
  
  for chromo in posDict:
    chromoGroup = group.create_group(chromo)
    chromoGroup.create_dataset('positions', dtype=uint32, data=posDict[chromo])
    chromoGroup.create_dataset('coords', dtype=float32, data=coordsDict[chromo])
  
  root.close()
  

def importCoords(filePath):

  root = File(filePath, mode='r')
  
  if 'structure' not in root:
    msg = 'File %s does not contain Nuc structure data'
    raise Exception(msg % filePath)
  
  group = root['structure']  
  posDict = {}
  coordsDict = {}
  
  for chromo in group:
    chromoGroup = group[chromo]

    posDict[chromo] = array(chromoGroup['positions'], uint32)
    coordsDict[chromo] = array(chromoGroup['coords'], float32)
    
  return posDict, coordsDict


def exportContacts(filePath, contactDict):
    
  root = File(filePath, mode='w')
  
  group = root.create_group('contacts')
  
  for chromoPair in contactDict:
    chrA, chrB = chromoPair
    
    if chrA in group:
      chromoGroup = group[chrA]
    else:
      chromoGroup = group.create_group(chrA)
      
    chromoGroup.create_dataset(chrB, dtype=uint32, data=contactDict[chromoPair])
  
  root.flush()
  root.close()


def importContacts(filePath, binSize=None):
 

  root = File(filePath, mode='r')

  if 'contacts' not in root:
    msg = 'File %s does not contain Nuc contact data'
    raise Exception(msg % filePath)
 
  group = root['contacts']  
  contactDict = {}

  for chrA in group:
    chromoGroup = group[chrA]
 
    for chrB in chromoGroup:
      data = array(chromoGroup[chrB], uint32)
      
      if binSize and binSize > 1:
        data = binSparseContacts(data, binSize)
      
      if data is None:
        continue
       
      contactDict[(chrA, chrB)] = data
 
  return contactDict


def exportContactMatrix(filePath, matrixDict, startDict, binSize):
    
  root = File(filePath, mode='w')
  
  group = root.create_group('contactMatrix')
  group.attrs['binSize'] = binSize

  for chromoPair in matrixDict:
    chrA, chrB = chromoPair
    
    if chrA in group:
      chromoGroup = group[chrA]
    else:
      chromoGroup = group.create_group(chrA)
      
    dataSet = chromoGroup.create_dataset(chrB, dtype=float32, data=matrixDict[chromoPair])
    dataSet.attrs['seqStarts'] = startDict[chromoPair]
    
  root.flush()
  root.close()


def importContactMatrix(filePath):

  root = File(filePath, mode='r')
  
  if 'contactMatrix' not in root:
    msg = 'File %s does not contain Nuc contact matrix data'
    raise Exception(msg % filePath)
  
  group = root['contactMatrix']
  binSize = group.attrs['binSize'] 
  startDict = {}
  matrixDict = {}
  
  for chrA in group:
    chromoGroup = group[chrA]
    
    for chrB in chromoGroup:
      chromoPair = (chrA, chrB)
      dataSet = chromoGroup[chrB]
      matrixDict[chromoPair] = array(dataSet, float32)
      startDict[chromoPair] = tuple(dataSet.attrs['seqStarts'])
  
  return matrixDict, startDict, binSize

