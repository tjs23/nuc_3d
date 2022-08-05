
def readVariableWigFile(filePath):

  fileObj = open(filePath)
  
  null = fileObj.readline()
  
  dataDict = {}
  chromo = None
  maxCount = 0
  
  for line in fileObj:
    line = line.strip()
    
    if line.startswith('variableStep'):
      key, chromo = line.split()[:2]
      chromo = chromo.split('=chr')[1]
      if chromo not in dataDict:
        dataDict[chromo] = []
    
    else:
      pos, count = line.split()
      count = int(count)
      
      if count:
        dataDict[chromo].append( (int(pos), count) )
        
        if count > maxCount:
          maxCount = count
  
  fileObj.close()
  
  for chromo in dataDict:
    dataDict[chromo].sort()
  
  return dataDict


def getMapabilityDict(mapFile, binSize, chromosome):

  fileObj = open(mapFile)
  null = fileObj.readline()

  mapability = {}
  for line in fileObj:
    chromo, pos, value = line.strip().split()
 
    if chromo == chromosome:
      key = int(pos) // binSize
 
      if key not in mapability:
        mapability[key] = []
 
      mapability[key].append( 1.0 - float(value) )

  for key in mapability:
    values = mapability[key]
    mapability[key] = sum(values)/float(len(values))

  return mapability

def readEnsembl(fileName='/home/tjs23/chromoVista/data/ensembl/EnsembleMm38.txt', chromosomes=None):

  fileObj = open(fileName, 'r')
  
  geneDict = {}
  
  if not chromosomes:
    chromosomes = [str(x) for x in range(19)] + ['X','Y']
  
  for chro in chromosomes:
    geneDict[chro] = {}
  
  for line in fileObj:
    
    data = line.strip().split('\t')
    
    if len(data) != 9:
      continue
    
    gId, tIdm, chro, st, en, desc, name, goname, godom = data
  
    if chro in chromosomes:
      st = int(st)
      en = int(en)
      pos = (st+en)/2
      geneDict[chro][pos] = (name, goname)

  return geneDict

def regionsToGenes(regions, chromosomes=None):

  geneDict = readEnsembl(chromosomes=chromosomes)
  genes = set()
  
  for chromo, posA, posB in regions:
    for posG in geneDict.get(chromo, []):
      if posA < posG < posB:
        info = geneDict[chromo][posG]
        genes.add((chromo, posG, info))
  
  genes = list(genes)
  genes.sort()
  
  return genes
