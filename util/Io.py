from os import path, stat, access, R_OK

def pathExists(filePath):

  try:
    return path.exists(filePath)
  
  except TypeError as err:
    filePath = filePath.replace('\0', '')
    return path.exists(filePath)
  
  except UnicodeEncodeError:
    return path.exists(filePath.encode('utf-8'))


def pathStat(filePath):

  try:
    return stat(filePath)
  
  except TypeError as err:
    filePath = filePath.replace('\0', '')
    return stat(filePath)

  except UnicodeEncodeError:
    return path.exists(filePath.encode('utf-8'))
  
  
def checkRegularFile(fileName):

  msg = ''
  
  if not pathExists(fileName):
    msg = 'File "%s" does not exist' % fileName
    return False, msg
  
  if not path.isfile(fileName):
    msg = 'Location "%s" is not a regular file' % fileName
    return False, msg
  
  if pathStat(fileName).st_size == 0:
    msg = 'File "%s" is of zero size '% fileName
    return False, msg
    
  if not access(fileName, R_OK):
    msg = 'File "%s" is not readable' % fileName
    return False, msg
  
  return True, msg


def guessTextFileFormat(filePath, separator, numCheckLines=None):

  if filePath.endswith('.gz'):
    import gzip
    fileObj = gzip.open(filePath, 'rU')
  
  else:
    fileObj = open(filePath, 'rU')
  
  if numCheckLines:
    lines = fileObj.readlines(numCheckLines)
  else:
    lines = fileObj.readlines()
  
  typeList = []
  convDict = {int:int, float:float, str:None}
  
  for line in lines:
    items = line.strip().split(separator)
    types = []
    
    for item in items:
      try:
        value = int(item)
      
      except ValueError:
        try:
          value = float(item)
        
        except ValueError:
          value = item
      
      types.append(type(value))
    
    typeList.append(types)  

  if typeList[0] == typeList[1]:
    skipFirst = False
    first = typeList[0]
  else:
    skipFirst = True    
    first = typeList[1]
  
  converters = [convDict[t] for t in first]
  
  for types in typeList[1:]:
    if types != first:
      for i, t1 in enumerate(first):
        if t1 is not types[i]:
          converters[i] = None
  
  if converters[:3] == (None, int, int):
    numDubious = 0
    
    for line in lines[1:]:
      items = line.strip().split(separator)
      a = int(items[1])
      b = int(items[2])
      
      if b < 0.5*a:
        numDubious += 1
    
    if numDubious > 0.8 * numCheckLines:
      haveEnds = False
    else:
      haveEnds = True  
    
  else:
    haveEnds = False
  
  if None not in converters:
    converters[0] = None
  
  return converters, skipFirst, haveEnds
  
  
def readListFile(fileName, converters, separator=None, skipFirst=True):
  
  dataList = []
  
  if fileName.endswith('.gz'):
    import gzip
    fileObj = gzip.open(fileName, 'rU')
    
  else:
    fileObj = open(fileName, 'rU')
  
  if skipFirst:
    header = fileObj.readline()     # Extract first line

  if separator:
    for line in fileObj:           # Loop through remaining lines
      line = line.rstrip()
      data = line.split(separator)

      for index, datum in enumerate(data):
        convertFunc = converters[index]
 
        if convertFunc:
          data[index] = convertFunc(datum)

      dataList.append(data)
  
  else:
    for line in fileObj:           # Loop through remaining lines
      data = line.split()

      for index, datum in enumerate(data):
        convertFunc = converters[index]
 
        if convertFunc:
          data[index] = convertFunc(datum)

      dataList.append(data)
  
   
  return dataList

WHITESPACE_AND_NULL =  set(['\x00', '\t', '\n', '\r', '\x0b', '\x0c'])

def isFileBinary(filePath):

  fileObj = open(filePath, 'rb')
  firstData = fileObj.read(1024)
  fileObj.close()
  
  testData = set([c for c in firstData]) - WHITESPACE_AND_NULL
  if min([ord(c) for c in testData]) < 32:
    # probably binary
    return True
    
  else:
    return False  
 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # #                       # # # # # # # # # # # # # # #
# # # # # # # # # # # # #   S A M   T O O L S   # # # # # # # # # # # # # # #
# # # # # # # # # # # # #                       # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def readBamSamFile(filePath):
  
  # TBD: Needs more checks to make sure it is a properly processed and paired BAM file
  
  from pysam import Samfile
  
  if isFileBinary(filePath):
    samFile = Samfile(filePath, 'rb')
    
  else:
    samFile = Samfile(filePath, 'r')
    
  chrDict = {}
  dataList = []
  dataListAppend = dataList.append
  for alignedRead in samFile:
    posA = alignedRead.pos
    posB = alignedRead.mpos
    chrIdA = alignedRead.tid
    chrIdB = alignedRead.mrnm
    
    if chrIdA in chrDict:
      chrA = chrDict[chrIdA]
    else:
      chrA = samFile.getrname(chrIdA)
      chrDict[chrIdA] = chrA
   
    if chrIdB in chrDict:
      chrB = chrDict[chrIdB]
    else:
      chrB = samFile.getrname(chrIdB)
      chrDict[chrIdB] = chrB
    
    dataListAppend((chrA, posA, chrB, posB))
  
  return dataList 
  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # #                       # # # # # # # # # # # # # # #
# # # # # # # # # # # # #  C S M  L E G A C Y   # # # # # # # # # # # # # # #
# # # # # # # # # # # # #                       # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

class AnnealProtocol:

  def __init__(self):
  
    self.steps = []

  def getCoolingScheme(self):
  
    data = []
 
    for i, step in enumerate(self.steps):
      datum = [i+1,
               step.startTemp,
               step.endTemp,
               step.coolSteps,
               step.mdSteps,
               step.mdTau,
               step.repulseScale]
 
      data.append(datum)
 
    return data


class AnnealStep:

  def __init__(self, protocol, startTemp, endTemp,
               coolSteps, mdSteps, mdTau, repulseScale):

    self.startTemp    = startTemp
    self.endTemp      = endTemp
    self.coolSteps    = coolSteps
    self.mdSteps      = mdSteps
    self.mdTau        = mdTau
    self.repulseScale = repulseScale

    self.protocol = protocol
    protocol.steps.append(self)

  def delete(self):
  
    self.protocol.steps.remove(self)
    del self

  def getIndex(self):
  
    return self.protocol.steps.index(self)
    
  def setIndex(self, i):
  
    steps = self.protocol.steps
    i = max(0, min(i, len(steps)-1))
    
    steps.remove(self)
    steps.insert(i,self)


class Sample:

  def __init__(self, name):
  
    self.name = name
    self.interactions = []
    self.chromoNodes = []
    self.chromoNodeDict = {}
    self.chromosomes = []
    self.annealProtocol = None
    self.dataLayers = []
    self.dataLayerDir = None

  def delete(self):
  
    for interaction in self.interactions:
      interaction.delete()
      
    for tagCluster in self.chromoNodes:
      tagCluster.delete()
          
    del self  
  
  def loadDataLayers(self):
  
    if self.dataLayerDir and os.path.exists(self.dataLayerDir):
      for fileName in os.listdir(self.dataLayerDir):
        if fileName.endswith('.dat'):
          layer = DataLayer(self, fileName)
          filePath = os.path.join(self.dataLayerDir, fileName)
          layer.readDatFile(filePath)
          
  def getCisInteractionLoci(self, refPoint=None, refWidth=1e6, chromosome='X'):

    loci = set()
 
    for interaction in self.interactions:
      if not interaction.numObs:
        continue
 
      for pair in interaction.nodePairs:
        nodeA, nodeB = pair.nodes
        chromoA = nodeA.chromosome
        chromoB = nodeB.chromosome
 
        if chromoA[-1] == 'B':
          continue
 
        if chromoB[-1] == 'B':
          continue
 
        if chromoA != chromoB:
          continue
 
        if chromoA != chromosome:
          continue
 
        if refPoint:
          if abs(nodeA.locus-refPoint) < refWidth:
            loci.add(nodeB.locus)
 
            if abs(nodeB.locus-refPoint) < refWidth:
              loci.add(nodeA.locus)
 
          elif abs(nodeB.locus-refPoint) < refWidth:
            loci.add(nodeA.locus)
 
          else:
            continue
 
        else:
          loci.add(nodeA.locus)
          loci.add(nodeB.locus)
 
    loci = list(loci)
    loci.sort()
 
    return loci
  
  def getChromosomes(self):
  
    names = self.chromoNodeDict.keys()
    names.sort()
    
    return names
  
  def getNumModels(self):
  
    m = 0
    
    for node in self.chromoNodes:
      coords = node.coords
         
      if coords:
        m = len(coords)
        break
    
    return m
  
class DataLayer:

  def __init__(self, sample, name, color='#FF0000',
               shape='star', scale=1.0, threshold=0.5,
               showText=True, details=None, isShown=False,
               isPersistent=True):
  
    self.sample = sample
    self.name = name
    self.color = color
    self.shape = shape
    self.scale = 1.0
    self.threshold = threshold
    self.showText = showText
    self.details = details
    self.isShown = isShown
    self.isPersistent = isPersistent
    
    self.dataDict = {}
    
    sample.dataLayers.append(self)
  
  def readDatFile(self, fileName):
  
    fileObj = open(fileName, 'rU')
    dataDict = self.dataDict
    
    
    for line in fileObj:
      
      if line.startswith('#DataLayer'):
        data = line.strip().split('|')
        name, color, shape, scale, threshold, showText, details = data[1:]
        
        scale = float(scale)
        threshold = float(threshold)
        showText = bool(showText)

        self.name = name
        self.color = color
        self.shape = shape or None
        self.scale = scale
        self.threshold = threshold
        self.showText = showText
        self.details = details
        continue

      data = line.strip().split()
      
      if len(data) > 3:
        chromo, pos, value = data[:3]
        anno = ' '.join(data[3:])
      else:
        chromo, pos, value = data
        anno = None
      
      pos = int(pos)
      value = float(value)
      
      key = (chromo, pos)
      dataDict[key] = (value, anno)
  
  def getCoords(self):
   
    sample = self.sample
    models = range(sample.getNumModels())
    threshold = self.threshold
    layerDict = self.dataDict
    chromoData = {}
    
    # Group data positions by chromosome
    for key in layerDict:
      chromo, pos = key
      val, ann = layerDict[key]
      
      if val < threshold:
        continue
      
      if chromo not in chromoData:
        chromoData[chromo] = []
      
      chromoData[chromo].append(pos) 
    
    for chromo in chromoData:
      chromoData[chromo].sort()
    
    # Get data coordinates in structure
    
    dataCoords = []
    nodeDict = sample.chromoNodeDict
    for chromo in nodeDict:
      if chromo not in chromoData:
        continue
      
      loci = [l for l in nodeDict[chromo] if nodeDict[chromo][l].coords]
      loci.sort()
 
      i = 0
      n = len(chromoData[chromo])
      pos = chromoData[chromo][0]
 
      for j, locus in enumerate(loci[:-1]):
        node = nodeDict[chromo][locus]
        locusB = loci[j+1]
        nodeB = nodeDict[chromo][locusB]
        
        while (pos < locusB) and (i<n):
          delta = float(locusB-locus)
          p = (pos-locus) / delta
          q = 1.0-p
          
          coords = []
          for m in models:
            x1, y1, z1 = node.coords[m]
            x2, y2, z2 = nodeB.coords[m]
            
            x = p*x1 + q*x2
            y = p*y1 + q*y2
            z = p*z1 + q*z2
            
            coords.append((x, y, z))
            
          dataCoords.append(coords)
          
          pos = chromoData[chromo][i]
          i += 1
    
    return dataCoords 
  
  def delete(self):
  
    self.sample.dataLayers.remove(self)
    
    del self
  
class ChromoNode:

  def __init__(self, sample, locus, chromosome, isRestrained=True, isFixed=False):
  
    self.sample = sample
    self.locus  = locus
    self.chromosome = chromosome
    self.nodePairs = {}
    self.coords = []
    self.density = 0.0
    self.depths = []
    self.mapability = 1.0
    self.isRestrained = isRestrained
    self.isFixed = isFixed
    self.name = '%s:%s' % (chromosome,locus)
    
    if not sample.chromoNodeDict.has_key(chromosome):
      sample.chromoNodeDict[chromosome] = {}
      sample.chromosomes.append(chromosome)
      sample.chromosomes.sort()
    
    if sample.chromoNodeDict[chromosome].get(locus):
      data = (chromosome,locus,sample.name)
      print('Duplicate ChromoNode %s:%s for Sample %s' % data)
      #raise Exception('Duplicate ChromoNode %s:%s for Sample %s' % data)
    
    sample.chromoNodeDict[chromosome][locus] = self
    sample.chromoNodes.append(self)  

    self.seqId = len(sample.chromoNodeDict[chromosome])

  def delete(self):
    
    for nodePair in self.nodePairs.keys():
      nodePair.delete()
    
    sample = self.sample
      
    sample.chromoNodes.remove(self)
    del sample.chromoNodeDict[self.chromosome][self.locus] 
    
    if not sample.chromoNodeDict[self.chromosome]:
      del sample.chromoNodeDict[self.chromosome]
      sample.chromosomes.remove(self.chromosome)

  def __str__(self):
  
    return '<ChromoNode %s %s:%s>' % (self.sample.name, self.chromosome, self.locus)
      
      
class TagInteraction: # PET restraint observation 

  def __init__(self, sample, numObs=None, targetDist=None, upperDist=None, lowerDist=None):
  
    self.sample = sample
    self.numObs = numObs
    self.targetDist = targetDist
    self.upperDist = upperDist
    self.lowerDist = lowerDist
    self.nodePairs = []
        
    self.sample.interactions.append(self)
  
  def deleteNodePairs(self):
  
    for nodePair in self.nodePairs:
       nodePair.delete()
    self.nodePairs = []
    
    
  def delete(self):
  
    self.sample.interactions.remove(self)
    
    for nodePair in self.nodePairs:
       nodePair.delete()
    
    del self

  def __str__(self):
    
    data = (self.numObs, self.targetDist, self.upperDist, self.lowerDist)
    txt = '<Interact %d T:%.2f U:%.2f L:%.2f>\n' % data
    for nodePair in self.nodePairs:
      txt += '    ' + str(nodePair) + '\n'
    
    return txt
      
class NodePair: # Ambiguous restraint item

  def __init__(self, interaction, nodeA, nodeB):
    
    assert nodeA.sample is nodeB.sample
    
    self.nodes = nodeA, nodeB
    self.interaction = interaction
    
    if nodeA.chromosome == nodeB.chromosome:
      self.distance = abs(nodeA.locus-nodeB.locus)
    else:
      self.distance = None

    interaction.nodePairs.append(self)
    nodeA.nodePairs[self] = True
    nodeB.nodePairs[self] = True
  
  def delete(self):
  
    self.interaction.nodePairs.remove(self)
    nodeA, nodeB = self.nodes
    del nodeA.nodePairs[self]
    del nodeB.nodePairs[self]
    del self

  def __str__(self):
  
    nodeA, nodeB = self.nodes
    return '<Pair> %s - %s' % (str(nodeA), str(nodeB))
    
    
DefaultAnnealProtocol = [[  5000, 3000,  2, 1000, 0.001, 0],
                         [  3000,  300,100, 1200, 0.001, 0],
                         [   300,   25,100, 3400, 0.001, 0.1],
                         [    25,   25,  1, 22500, 0.001,0.1],
                         [    10,   10,  1, 22500, 0.001,0.1],
                         [  0.01, 0.01,  1, 22500,0.0005,0.1]]

TRUE = 'True'

def loadCsmDataSet(fileName):
  
  file = open(fileName)
  
  """
  try:
    sample = cPickle.load(file)
  except EOFError:
    showWarning('Warning','Data load failed', parent=self)
    sample = None
  """
  
  sample = None
  proto = None
  dataLayer = None
  
  line = file.readline()
  while line:
    a = line.split()
    key = a[0]
    
    if key == 'Sample':
      name = line[len(key)+1:].strip()
      sample = Sample(name)
    
    if sample:
      getNode = sample.chromoNodeDict.get
    
      if key == 'AnnealProtocol':
        proto = AnnealProtocol()
        sample.annealProtocol = proto
    
      elif (key == 'AnnealStep') and proto:
        step = AnnealStep(proto, float(a[1]), float(a[2]),
                          int(float(a[3])), int(float(a[4])),
                          float(a[5]), float(a[6]))

      elif key == 'DataLayerDir':
        dataLayerDir = ' '.join(a[1:])
        sample.dataLayerDir = dataLayerDir

      elif key == 'DataLayer':
        name, color, shape, scale, threshold, showText = a[1:7]
        details = ' '.join(a[7:-1])
        isShown = a[-1]
        if details == 'None':
          details = None
          
        dataLayer = DataLayer(sample, name, color, shape,
                              float(scale), float(threshold),
                              showText==TRUE, details, isShown==TRUE)
                             
      elif (key == 'LayerDatum') and dataLayer:
        chromo, pos, value, anno = a[1:]
        key = (chromo, int(pos))
        
        if anno == 'None':
          anno = None
        
        dataLayer.dataDict[key] = (float(value), anno)
  
      elif key == 'ChromoNode':
        # sample locus chromosome x1:y1:z1,x2:y2:z2,x3:y3:z3 den map d1:d2:d3
        
        if a[5] == 'None':
          coords = []
        
        else:
          coords = [[float(c) for c in xyz.split(':')] for xyz in a[5].split(',')]
        
        if a[4] == 'True':
          isFixed = True
        else:
          isFixed = False

        if a[3] == 'True':
          isRestrained = True
        else:
          isRestrained = False
        
        if len(a) > 6:
          density = float(a[6])
        else:
          density = 0.0
        
        if len(a) > 7:
          mapability = float(a[7])
        else:
          mapability = 1.0
        
        if len(a) > 8:
           if a[8] == 'None':
             depths = []
           else:
             depths = [float(d) for d in a[8].split(',')]
        else:
          depths = []
         
        tag = ChromoNode(sample, int(a[1]), a[2], isRestrained, isFixed)
        tag.coords = coords
        tag.density = density
        tag.mapability = mapability
        tag.depths = depths

      elif key == 'TagInteraction':
      
        if a[1] == 'None':
          obs = None
        else:
          obs = int(float(a[1]))
      
        interation = TagInteraction(sample, obs)
       
        if a[2] == 'None':
          t = None
        else:
          t = float(a[2])
        
        if a[3] == 'None':
          u = None
        else:
          u = float(a[3])
        
        if a[4] == 'None':
          l = None
        else:
          l = float(a[4])
        
        interation.targetDist = t
        interation.upperDist = u
        interation.lowerDist = l

      elif key == 'NodePair':
        
        locA = int(a[1])
        chrA = getNode(a[2])
        locB = int(a[3])
        chrB = getNode(a[4])
        
        if chrA and chrB:
          nodeA = chrA.get(locA)
          nodeB = chrB.get(locB)
  
          if nodeA and nodeB:
            pair = NodePair(interation, nodeA, nodeB)
            
  
    line = file.readline()
  
  
  return sample
            
def saveCsmDataSet(sample, fileName):

  fileObj = open(fileName, 'w')
  #cPickle.dump(sample, fileObj)

  write = fileObj.write

  write('Sample %s\n' % sample.name)
  
  if sample.annealProtocol:
    write('AnnealProtocol\n')
    
    for step in sample.annealProtocol.steps:
      data = (step.startTemp,
              step.endTemp,
              step.coolSteps,
              step.mdSteps,
              step.mdTau,
              step.repulseScale)
      
      write('AnnealStep %f %f %d %d %f %f\n' % data)
  
  if sample.dataLayerDir:
    write('DataLayerDir %s\n' % sample.dataLayerDir)
      
  for dataLayer in sample.dataLayers:
    if not dataLayer.isPersistent:
      continue
  
    name = dataLayer.name.replace(' ', '_')
    data = (name, dataLayer.color,
            dataLayer.shape, dataLayer.scale,
            dataLayer.threshold, str(dataLayer.showText),
            dataLayer.details or 'None', dataLayer.isShown)
    write('DataLayer %s %s %s %f %f %s %s %s\n' % data)
    
    for key in dataLayer.dataDict:
      chromo, pos = key
      value, anno = dataLayer.dataDict[key]
      
      data = (chromo, pos, value, anno or 'None')
      write('LayerDatum %s %d %f %s\n' % data)
   
  for node in sample.chromoNodes:
  
    if node.coords:
      coords = ','.join(['%f:%f:%f' % tuple(x) for x in node.coords])
    else:
      coords = 'None'
    
    if node.depths:
      depths = ','.join(['%f' % x for x in node.depths])
    else:
      depths = 'None'
    
    data = (node.locus, node.chromosome, node.isRestrained,
            node.isFixed, coords, node.density,
            node.mapability, depths)
    write('ChromoNode %d %s %s %s %s %f %f %s\n' % data)
  
  for interaction in sample.interactions:
    
    if interaction.numObs is None:
      obs = 'None'
    else:
      obs = '%d' % interaction.numObs

    if interaction.targetDist is None:
      targetDist = 'None'
    else:
      targetDist = '%f' % interaction.targetDist

    if interaction.upperDist is None:
      upperDist = 'None'
    else:
      upperDist = '%f' % interaction.upperDist

    if interaction.lowerDist is None:
      lowerDist = 'None'
    else:
      lowerDist = '%f' % interaction.lowerDist

    data = (obs, targetDist, upperDist, lowerDist)
            
    write('TagInteraction %s %s %s %s\n' % data)

    for nodePair in interaction.nodePairs:
    
      nodeA, nodeB = nodePair.nodes
      
      data = (int(nodeA.locus), nodeA.chromosome,
              int(nodeB.locus), nodeB.chromosome)
              
      write('NodePair %s %s %s %s\n' % data)
            
  
  fileObj.close()

  
def loadCsmInteractions(fileName, sampleName, minSeparation=1):
  
  def _readColumnsFile(fileName):

    data = []
 
    fileObj = open(fileName)
    null = fileObj.readline()

    for line in fileObj:
      if line[0] not in '!*#':
        datum  = line.split()
 
        if len(datum) > 3:
          data.append(datum)
 
    fileObj.close()
 
    return data

  data = _readColumnsFile(fileName)
  
  sample = Sample(sampleName)

  nodeDict = {}
  
  
  if len(data[0]) == 4:
  
    for chrA, locA, chrB, locB in  data:
      chrA = chrA.upper()
      chrB = chrB.upper()
 
      bpA = int(locA)
      bpB = int(locB)

      if chrA == chrB:
        if abs(bpA-bpB) < minSeparation:
          continue

      # TBD ambig Chromo
      
      if (chrA, bpA) == (chrB, bpB):
        continue


      if chrA in ('X','Y'):
        if chrB in ('X','Y'):
          loci1 =[(chrA, bpA),]
          loci2 =[(chrB, bpB),]
          
        else:  
          loci1 =[(chrA, bpA),]
          loci2 =[(chrB+'A', bpB), (chrB+'B', bpB),]
      
      elif chrB in ('X','Y'):
        loci1 =[(chrA+'A', bpA), (chrA+'B', bpA)]
        loci2 =[(chrB, bpB),]

      else:
        loci1 =[(chrA+'A', bpA), (chrA+'B', bpA)]
        loci2 =[(chrB+'A', bpB), (chrB+'B', bpB)]

      pairs = []
      for loc1 in loci1:
        key1 = '%s:%d' % loc1
        node1 = nodeDict.get(key1)
        
        if node1 is None:
          chromo1, posn1 = loc1
          node1 = ChromoNode(sample, posn1, chromo1)
          nodeDict[key1] = node1
      
        for loc2 in loci2:
          key2 = '%s:%d' % loc2
          node2 = nodeDict.get(key2)

          if node2 is None:
            chromo2, posn2 = loc2
            node2 = ChromoNode(sample, posn2, chromo2)
            nodeDict[key2] = node2
          
          pairs.append((node1, node2))


      ti = TagInteraction(sample, 10)
      
      for node1, node2 in pairs:
        NodePair(ti, node1, node2)
  
  elif len(data[0]) == 5:
  
    for chrA, locA, chrB, locB, numObs in  data:
 
      chrA = chrA.upper()
      chrB = chrB.upper()
      
      if 'X' not in (chrA, chrB):
        continue
      
      #if chrA not in ('X','3','4'):
      #  continue
 
      #if chrB not in ('X','3','4'):
      #  continue
       
      numObs = int(numObs)
 
      bpA = int(locA)
      bpB = int(locB)

      if chrA == chrB:
        if abs(bpA-bpB) < minSeparation:
          continue

      # TBD ambig Chromo
      
      if (chrA, bpA) == (chrB, bpB):
        continue


      if chrA in ('X','Y'):
        if chrB in ('X','Y'):
          loci1 =[(chrA, bpA),]
          loci2 =[(chrB, bpB),]
          
        else:  
          loci1 =[(chrA, bpA),]
          loci2 =[(chrB+'A', bpB), (chrB+'B', bpB),]
      
      elif chrB in ('X','Y'):
        loci1 =[(chrA+'A', bpB), (chrA+'B', bpB)]
        loci2 =[(chrB, bpA),]

      else:
        loci1 =[(chrA+'A', bpA), (chrA+'B', bpA)]
        loci2 =[(chrB+'A', bpB), (chrB+'B', bpB)]


      pairs = []
      for loc1 in loci1:
        key1 = '%s:%d' % loc1
        node1 = nodeDict.get(key1)
        
        if node1 is None:
          chromo1, posn1 = loc1
          node1 = ChromoNode(sample, posn1, chromo1)
          nodeDict[key1] = node1
      
        for loc2 in loci2:
          key2 = '%s:%d' % loc2
          node2 = nodeDict.get(key2)

          if node2 is None:
            chromo2, posn2 = loc2
            node2 = ChromoNode(sample, posn2, chromo2)
            nodeDict[key2] = node2
          
          pairs.append((node1, node2))


      if numObs > 0:
        ti = TagInteraction(sample, numObs)
        
        for node1, node2 in pairs:
          NodePair(ti, node1, node2)
      
        
  elif len(data[0]) == 7:
    for fendA, chrA, locA, fendB, chrB, locB, numObs in data:
 
      chrA = chrA.upper()
      chrB = chrB.upper()
 
      #if chrA not in ('X','3','4'):
      #  continue
 
      #if chrB not in ('X','3','4'):
      #  continue
       
      numObs = int(numObs)
 
      bpA = int(locA)
      bpB = int(locB)

      if chrA == chrB:
        if abs(bpA-bpB) < minSeparation:
          continue

      # TBD ambig Chromo
      
      if (chrA, bpA) == (chrB, bpB):
        continue


      if chrA in ('X','Y'):
        if chrB in ('X','Y'):
          loci1 =[(chrA, bpA),]
          loci2 =[(chrB, bpB),]
          
        else:  
          loci1 =[(chrA, bpA),]
          loci2 =[(chrB+'A', bpB), (chrB+'B', bpB),]
      
      elif chrB in ('X','Y'):
        loci1 =[(chrA+'A', bpB), (chrA+'B', bpB)]
        loci2 =[(chrB, bpA),]

      else:
        loci1 =[(chrA+'A', bpA), (chrA+'B', bpA)]
        loci2 =[(chrB+'A', bpB), (chrB+'B', bpB)]


      pairs = []
      for loc1 in loci1:
        key1 = '%s:%d' % loc1
        node1 = nodeDict.get(key1)
        
        if node1 is None:
          chromo1, posn1 = loc1
          node1 = ChromoNode(sample, posn1, chromo1)
          nodeDict[key1] = node1
      
        for loc2 in loci2:
          key2 = '%s:%d' % loc2
          node2 = nodeDict.get(key2)

          if node2 is None:
            chromo2, posn2 = loc2
            node2 = ChromoNode(sample, posn2, chromo2)
            nodeDict[key2] = node2
          
          pairs.append((node1, node2))


      if numObs > 0:
        ti = TagInteraction(sample, numObs)
        
        for node1, node2 in pairs:
          NodePair(ti, node1, node2)
          
  protocol = AnnealProtocol()
  sample.annealProtocol = protocol
    
  for startT, endT, coolSteps, mdSteps, mdTau, rep in DefaultAnnealProtocol:
    step = AnnealStep(protocol, startT, endT, coolSteps, mdSteps, mdTau, rep)

  return sample



