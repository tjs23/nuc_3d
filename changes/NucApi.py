import os, sys, shutil, time, colorsys, re
from datetime import datetime

from h5py import File, Group, special_dtype, h5o
from numpy import array, float32, uint32, int32, int16, uint8, ones, zeros, vstack, arange
from numpy import dstack, abs, sqrt, dot, append, random, empty, log, log2, outer, concatenate
from numpy import power, eye, argsort, argwhere, cumsum, sqrt, hstack, exp, clip, ndarray, string_

from random import seed, shuffle, randint, uniform
from math import sin, cos, atan, ceil

from solve.SimAnneal import runDynamics
from util.NaturalSort import naturalKey

# TBD Command line
#
# nuc3d f1.hdf5 merge f2.hdf5 f3.hdf5
# nuc3d f1.hdf5 align (f2.hdf5, ...)
# nuc3d f1.hdf5 anneal (tMax, tMin, nSteps)
# nuc3d f1.hdf5 set [backbone] value
# nuc3d f1.hdf5 list [id|info|chr|point|contact|coord|data|transform|viol|domain] (value)
# nuc3d f1.hdf5 num [chr|model|point|contact|data|viol] 
# nuc3d f1.hdf5 matrix [contact|distance]
# nuc3d f1.hdf5 image [model|contact|distance] (format)
# nuc3d f1.hdf5 movie (style) (format) (outfile)
# nuc3d f1.hdf5 export [pdb|wig|bed] (value|model|chr) (outFile)
# nuc3d f1.hdf5 copy f2.hdf5 # New ID
# nuc3d f1.hdf5 valid
# nuc3d f1.hdf5 print*
# nuc3d f1.hdf5 show*
# nuc3d f1.hdf5 text*
# nuc3d f1.hdf5 test

VL_str = special_dtype(vlen=str)
PROGRAM = 'Nuc3D'
FILE_EXT = '.nuc'
VERSION = 0.1

PI = 3.14159265358979323846
TAU = 2.0 * PI

def load(fileName):
  
  root = Nucleus(fileName)
  
  return root

def readContacts(fileName, outFileName=None):
  
  if not fileName:
    outFileName = os.path.splitext(fileName)[0] + FILE_EXT
  
  nuc = Nucleus(outFileName)
  nuc.importContactFile(fileName, binSize=int(5e5))
  
  nuc.setRestraints(binned=False) # Default value for all chromosomes
  #nuc.setMapability() # Default source, based on genome build
  nuc.setSpiralCoords()
  
  return nuc

def merge(inFileNameA, inFileNameB, outFileName):

  pass

def modelAlign(fileNameA, fileName):

  pass

DERIVED  = 'derived'
EXTERNAL = 'external'
INNATE   = 'innate'

COLOR_MODES = ['Seq. Position','Chromsome ID','Density',
               'Model number','Data track', 'Faint', 'RMSD']

DISPLAY_MODES = ['Ball and Stick', 'Line', 'Tube', 'Surface']

DATA_TRACK_SYMBOLS = ('Circle','Star','Square','Triangle','Hexagon')

DATA_TRACK_PEAK_TYPES = ('Histogram', 'Region strip', 'Contact box')

DATA_TRACK_SYMBOL_PARAMS = ((TAU/10.0, 10), (TAU/2.5, 5), (TAU/4.0, 4),
                             (TAU/3.0, 3), (TAU/6.0, 6)) # (dAngle, nPoints)

RESTRAINT_COLOR_MODES = ['Cis/Trans', 'Distance', 'Seq. Separation', 'Chromosome ID']
  
STRUC_CALC_DEFAULTS = {'numModels':1, 'bboneReg':0, 'bboneSpace':100, 'restrBinned':0,
                       'powerLaw':-0.33, 'seqUnitScale':10, 'restrDist':1.0,
                       'restrErr':0.2, 'tempMax':5000, 'tempMin':10, 'tempSteps':256,
                       'dynSteps':256, 'hierProtocol':1, 'hierStart':4,
                       'hierSteps':1, 'startStruc':0, 'randRad':100.0, 'randSeed':7}

DATA_TRACK_CACHE_BIN = 100000


"""
/ - id
  sample/ - name
    protocol/
    tissue/
    organism/
  
  display/
    colorSchemes/
      
  contacts/
    original/
      <groupName>/ - isSingleCell, filePath, binSize
        chrA/
          chrB/ - posA, posB, numObs 
          
    working/
      <groupName>/
        <chrA>/
          <chrB>/ - posA, posB, numObs, model
    
  structures/
    calculation/
  
    chromosomes/
      <chrName>/ - display, color
        backbone  - 0/1
        positions - seqPos
 
    restraints/
      <chrA>/
        <chrB>/ idxA, idxB, upper, lower, weight, ..
  
    transforms/
      global - x, y, z
      models - model, x, y, z
      
    coords/
      <chr>/ model, x, y, z

  dataTracks/
    derived/ # From contacts/structures
      <code>/ - display, options
        <chr>/ - indicesMb
          regions - first, last
          values  - origVal, normVal
          annotations (optional)
          models (optional)

    innate/  # From sequence/molecule
      <code>/
        <chr>/ ...
    
    external/
      <code>/
        <chr>/ ...
    
  images/
    foci/
      <code> - x, y, z, shape, intensity

    grid/ - origin, gridSizes
      <code> - value [i,j,k]
   
    coords/
      <code> - x, y, z, value 

"""

def backwardCompatibility(root):
  # backward compatibility after model changes
  
  if 'tissueRef' in root.attrs:
    root.attrs['experimentRef'] = string_(root.attrs['tissueRef'])
    del root.attrs['tissueRef']
    
  if 'dataTracks' not in root:
    root.create_group('dataTracks')
    
    if 'genomeData' in root:
      for group in root['genomeData']:
        root['dataTracks'].create_group(group)
        
        for code in root['genomeData'][group]:
          root['dataTracks'][group].create_group(code)
          
          if 'options' not in root['genomeData'][group][code].attrs:
            root['dataTracks'][group][code].attrs['stranded'] = 0
            root['dataTracks'][group][code].attrs['options']  = (0, 0, 1, 0, 0, 0, 0, 0) # shown, labels, symbol, trackType
            root['dataTracks'][group][code].attrs['display']  = (0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0) # r, g, b, scale, threshold

          else:
            root['dataTracks'][group][code].attrs['stranded'] = root['genomeData'][group][code].attrs['stranded']
            root['dataTracks'][group][code].attrs['options']  = root['genomeData'][group][code].attrs['options']
            root['dataTracks'][group][code].attrs['display']  = root['genomeData'][group][code].attrs['display']
 
          
          for chromo in list(root['genomeData'][group][code]):
            root['dataTracks'][group][code].create_group(chromo)
            
            prevGrp = root['genomeData'][group][code][chromo]
            newGrp  = root['dataTracks'][group][code][chromo]
            
            data = array(prevGrp['regions'], uint32)
            newGrp.create_dataset('regions', dtype=uint32, data=data)
             
            data = array(prevGrp['values'], float32)
            newGrp.create_dataset('values', dtype=float32, data=data)
             
            if 'annotations' in prevGrp:
              newGrp.create_dataset('annotations', dtype=VL_str, data=prevGrp['annotations'])
 
            if 'models' in prevGrp:
              data = array(prevGrp['models'], uint32)
              newGrp.create_dataset('models', dtype=uint32, data=data)
            
            del root['genomeData'][group][code][chromo]
                      
      del root['genomeData']     
  
  
  for group in root['dataTracks']:
    for code in root['dataTracks'][group]:
      for chromo in root['dataTracks'][group][code]:
        attrs = root['dataTracks'][group][code][chromo].attrs
        
        if 'indicesMb' not in attrs:
          attrs['indicesMb'] = []
        
  
  if 'display' not in root:
    root.create_group('display')

  if 'original' not in root['contacts']:
    root['contacts'].create_group('original')
    
  if 'working' not in root['contacts']:
    root['contacts'].create_group('working')
  
  if 'voxels' in root:
    del root['voxels']
    group = root.create_group('images')
    group.create_group('foci')
    group.create_group('grid')
    group.create_group('coords')
    
  original = root['contacts']['original']
  groupNames = [a for a in root['contacts'] if a not in ('original', 'working')]
  
  for groupName in groupNames:
    
    oldGroup = root['contacts'][groupName]
    newGroup = original.create_group(groupName)
    
    if groupName == 'singleCell':
      newGroup.attrs['isSingleCell'] = 1
    else:
      newGroup.attrs['isSingleCell'] = 0        
    
    if 'cis' in oldGroup: # was cis/trans split to singleCell/population contacts
      for chrA in oldGroup['cis']:
        data = array(oldGroup['cis'][chrA])
        subGroup = newGroup.create_group(chrA)
        subGroup.create_dataset(chrA, dtype=uint32, data=data,
                                compression='gzip')
        
      for chrA in oldGroup['trans']:
        if chrA in newGroup:
          subGroup = newGroup[chrA]
        else:
          subGroup = newGroup.create_group(chrA)
      
        for chrB in oldGroup['trans'][chrA]:
          data = array(oldGroup['trans'][chrA][chrB])
          subGroup.create_dataset(chrB, dtype=uint32, data=data,
                                  compression='gzip')
      
      del oldGroup['cis']
      del oldGroup['trans']
    
    else: # was contacts stored as 'chrA chrB' keys in singleCell/poopulation contacts
      keys = [k for k in oldGroup]
      
      for key in keys:
        if not len(oldGroup[key]):
          del oldGroup[key] 
          continue

        chromos = sorted(key.split())
        data = array(oldGroup[key]).T
        
        chrA = chromos[0]
        chrB = chromos[-1]
        
        if chrA in newGroup:
          subGroup = newGroup[chrA]
        else:
          subGroup = newGroup.create_group(chrA)
        
        subGroup.create_dataset(chrB, dtype=uint32, data=data,
                                compression='gzip')
      
        del oldGroup[key] 
    
    del root['contacts'][groupName]
  
  group = root['structures']['restraints']
  if 'cis' in group:
    
    for chrA in group['cis']:
      data = array(group['cis'][chrA])
      subGroup = group.create_group(chrA)
      subGroup.create_dataset(chrA, dtype=float32, data=data)
    
    for chrA in group['trans']:
      if chrA in group:
        subGroup = group[chrA]
      else:
        subGroup = group.create_group(chrA)
    
      for chrB in group['trans'][chrA]:
        data = array(group['trans'][chrA][chrB])
        subGroup.create_dataset(chrB, dtype=float32, data=data)
    
    del group['cis']
    del group['trans']
    
  else:
    keys = [x for x in group if ' ' not in x] # single chromos first
    for key in keys:
      item = group[key]
      
      if not isinstance(item, Group): # Was stored directly with 'chrA chrB' keys
        chromos = key.split()
      
        if len(group[key]):
          chrA = chromos[0]
          chrB = chromos[-1]
          data = array(group[key]).T
          del group[key]

          if chrA in group:
            subGroup = group[chrA]
          else:
            subGroup = group.create_group(chrA)
          
          subGroup.create_dataset(chrB, dtype=float32, data=data)
        
        else:
          del group[key]         

    keys = [x for x in group if ' ' in x]
    for key in keys:
      item = group[key]
      
      if not isinstance(item, Group): # Was stored directly with 'chrA chrB' keys
        chromos = key.split()
      
        if len(group[key]):
          chrA = chromos[0]
          chrB = chromos[-1]
          data = array(group[key]).T
          del group[key]

          if chrA in group:
            subGroup = group[chrA]
          else:
            subGroup = group.create_group(chrA)
          
          subGroup.create_dataset(chrB, dtype=float32, data=data)
        
        else:
          del group[key]         
      
  for chromo in root['structures']['chromosomes']:
    if 'active' in root['structures']['chromosomes'][chromo]:
      del root['structures']['chromosomes'][chromo]['active']
      
    if 'colours' in root['structures']['chromosomes'][chromo]:
      del root['structures']['chromosomes'][chromo]['colours']
      
    if 'annotations' in root['structures']['chromosomes'][chromo]:
      del root['structures']['chromosomes'][chromo]['annotations']
  
  for section in root['contacts']:
    for groupName in root['contacts'][section]:
      attrs = root['contacts'][section][groupName].attrs
      
      if 'isSingleCell' not in attrs:
        attrs['isSingleCell'] = 1
    
      if 'binSize' not in attrs:
        if attrs['isSingleCell']:
          attrs['binSize'] = 1
        
        else:
          group = root['contacts'][section][groupName]
          binSize = 1
          for chrA in group:
            for chrB in group[chrA]:
              if chrA == chrB:
                dataset = array(group[chrA][chrB])
                pos = sorted(set(dataset[0].tolist()))
                
                if len(pos) > 1:
                  pos = array(pos)
                  deltas = pos[1:]-pos[:-1]
                  binSize = deltas.min() 
                  break
                  
            else:
              continue
            break        
          
          attrs['binSize'] = binSize
          

def usingTimer(function):

  def timer(*args, **kw):
    start = time.time()
    output = function(*args, **kw)
    end = time.time()
    print('Function "%s" took %.4f seconds' % (str(function), end-start))
    
    return output
  
  return timer


class AnnealStage:

  def __init__(self, schedule, domainDict, bboneSep, tempStart, tempEnd,
               tempSteps, dynSteps, timeStep, useTrans):

    self.schedule = schedule
    self.domainDict = domainDict  # per chromosome, list of regions (start, end)
    self.bboneSep = bboneSep      # If not binned set < 2 for no backbone spacers
    self.tempStart = tempStart
    self.tempEnd = tempEnd
    self.tempSteps = tempSteps
    self.dynSteps = dynSteps
    self.timeStep = timeStep
    self.useTrans = useTrans
    
    schedule.annealStages.append(self)
    schedule.numAnnealStages += 1
  
  
  def delete(self):
  
    self.schedule.numAnnealStages -= 1
    self.schedule.annealStages.remove(self)
  
  
  def getTempSchedule(self, scheme='exp'):
    
    # construct temperature sub-schedule
    temps = []
    n = self.tempSteps

    if scheme == 'exp':
      decay = log(self.tempStart/self.tempEnd)

      for i in range(n):
        frac = i/float(n)
        temp = self.tempStart * exp(-decay*frac)
        temps.append( temp )
    
    else: # linear
      
      for i in range(n):
        frac = (i-n)/float(n)
        temp = (1.0-frac) * self.tempStart + frac * self.tempEnd
        temps.append( temp )
 
    return temps
   
   
  def getRepulsionSchedule(self, scheme='sigmoid'):

    # construct trepulsion sub-schedule
    repScales = []
    n = self.tempSteps
    
    if scheme == 'sigmoid':
      adj = 1.0 / atan(10.0)

      for i in range(n):
        frac = i/float(n)
        repScale = 0.5 + adj * atan(frac*20.0-10) / PI
        repScales.append(repScale)
    
    else: # linear
      
      for i in range(n):
        frac = 0.99 * (i-n)/float(n)
        repScales.append(0.01 + frac)
        
    return repScales
     
        
class StrucCalcParams:
           
  def __init__(self, distPowerLaw=-0.33, distLower=0.8, distUpper=1.2,
               bboneLower=0.1,  bboneUpper=1.1, seqScale=int(1e4),         
               randSeed=None, randStart=True, randWalk=False, randRad=100.0,
               minNumObs=1, maxPopDist=5.0):
    
    self.distPowerLaw = distPowerLaw
    self.distLower = distLower
    self.distUpper = distUpper
    self.bboneLower = bboneLower
    self.bboneUpper = bboneUpper
    self.seqScale = seqScale
    self.randSeed = randSeed
    self.randStart = randStart
    self.randWalk = randWalk
    self.randRad = randRad
    self.minNumObs = minNumObs
    self.maxPopDist = maxPopDist
    
    self.numAnnealStages = 0
    self.annealStages = []
    self._istep = 0
  
  def __iter__(self):
    
    return self # This is an iterator object

  def __next__(self): # For python 3

    return self.next()
   
  def next(self): # For Python 2
  
    self._istep += 1
    
    if self._istep < self.numAnnealStages:
      return self.annealStages[self._istep]
      
    else:
      self._istep = 0
      raise StopIteration()
        
  
  def addAnnealStage(self, domainDict={}, bboneSep=[int(5e4), int(5e4)], tempStart=5000.0, tempEnd=10.0,
                     tempSteps=100, dynSteps=100, timeStep=0.001, useTrans=True):
    
    step = AnnealStage(self, domainDict, bboneSep, tempStart, tempEnd,
                       tempSteps, dynSteps, timeStep, useTrans)
    
    return step
    
  def getStep(self, i):
  
    if i <  self.numAnnealStages:
      return self.annealStages[i]
  

class Nucleus:
  """Top-level containment object, corresponds to one HDF5 format binary file"""
  
  notifiers = set()
  isModified = True
  
  def __init__(self, fileName=None, version=VERSION, experimentRef=None, genomeRef=None, mode='a'):
    
    if (fileName is None) or self._checkFileName(fileName):
      
      if fileName and self._filePathExists(fileName):
        try:
          self.root = root = File(fileName, mode=mode)
        
        except Exception as err:
          print('File "%s" not a valid %s file' % (fileName, PROGRAM))
          raise Exception(err)

        root.attrs['id'][2] = time.time()
        
      else:
        self.root = root = File(fileName, mode='a')
        
        hierarchy = (('contacts',   ('original', 'working')),
                     ('display',    ()),
                     ('dataTracks', (DERIVED, INNATE, EXTERNAL)),
                     ('sample',     ('protocol', 'organism', 'tissue')),
                     ('structures', ('chromosomes', 'restraints', 'transforms', 'coords')),
                     ('images',     ('foci', 'grid', 'coords'))
                     )
        
        for parent, children in hierarchy:
          group = root.create_group(parent)
        
          for child in children:
            group.create_group(child)
        
        now = int(time.time())
        random.seed(now)        
        
        data = array([random.random(), now, now], float32)
        root.attrs['id'] = data
        
        name = fileName if fileName else 'Unknown'
        root['sample'].attrs['name'] = string_(name)  
          
    else:
      raise Exception('%s file name "%s" not valid' % (PROGRAM, fileName))
    
    backwardCompatibility(self.root)
    
    self._setLinkAttrs()
    self._chromoLimitCache = {}
    self._dataTrackHistogramCache = {}
    self._contactsCache = {}
    self._chromosomes = None
    self._genomeRefNuc = None
    self._experimentRefNuc = None
    
    # Root attributes
    
    if experimentRef:
      root.attrs['experimentRef'] = string_(experimentRef)
    elif experimentRef is False:
      experimentRef = None
    elif 'experimentRef' in root.attrs:
      experimentRef = root.attrs['experimentRef']
      
      if experimentRef == 'None':
        experimentRef = None
      
    if genomeRef:
      root.attrs['genomeRef'] = string_(genomeRef)
    elif genomeRef is False:
      genomeRef = None
    elif 'genomeRef' in root.attrs:
      genomeRef = root.attrs['genomeRef']
     
      if genomeRef == 'None':
        genomeRef = None
       
    self.version = version # Runtime
    self.fileName = fileName
    self.setExperimentRef(experimentRef) # Runtime - depends on local context
    self.setGenomeRef(genomeRef) # Runtime - depends on local context

  
  def updateProxyDataTracks(self, remoteNuc, source):
  
    rNuc = remoteNuc
    rDataGroup = rNuc._getDataTrackGroup(source)
    lDataGroup = self._getDataTrackGroup(source)
    
    for code in rDataGroup:
      if code not in lDataGroup:
        rLayer = rDataGroup[code]
        lLayer = self._getGroup(code, lDataGroup)
        
        lLayer.attrs['stranded'] = rLayer.attrs['stranded']
        lLayer.attrs['options']  = rLayer.attrs['options'] 
        lLayer.attrs['display']  = rLayer.attrs['display'] 
        
        for chromo in rLayer:
          subGroup = self._getGroup(chromo, lLayer)
          
          self._setData('regions', subGroup, uint32, [])
          self._setData('values', subGroup, float32, [])
          break # only one chromo required
              
    return rNuc
    
      
  def importNucDataTrack(self, remoteNuc, source, code):
  
    rNuc = remoteNuc    
    rDataGroup = rNuc._getDataTrackGroup(source)
    lDataGroup = self._getDataTrackGroup(source)
    
    if code in rDataGroup:
      rLayer = rDataGroup[code]
      
      if code not in lDataGroup:
        lLayer = self._getGroup(code, lDataGroup)
        
        lLayer.attrs['stranded'] = rLayer.attrs['stranded']
        lLayer.attrs['options']  = rLayer.attrs['options'] 
        lLayer.attrs['display']  = rLayer.attrs['display'] 
      
      else:
        lLayer = self._getGroup(code, lDataGroup)
        
      for chromo in rLayer:
        lSubGroup = self._getGroup(chromo, lLayer)
        rSubGroup = rLayer[chromo]
        
        regionArray = array(rSubGroup['regions'], uint32)
        valueArray = array(rSubGroup['values'], float32)
        
        self._setData('regions', lSubGroup, uint32, regionArray)
        self._setData('values', lSubGroup, float32, valueArray)
        
        if ('annotations' in rSubGroup) and len(rSubGroup['annotations']):
          self._setData('annotations', lSubGroup, VL_str, string_(rSubGroup['annotations']))
        
        if 'models' in rSubGroup:
          modelArray = array(rSubGroup['models'], uint32)
          self._setData('models', lSubGroup, uint32, modelArray)
    
        key = (source, code, chromo)
        if key in self._dataTrackHistogramCache:
          del self._dataTrackHistogramCache[key]
      
      self._setDataTrackIndicesMb(source, code)     
      
        
  def setExperimentRef(self, filePath=None):
    """Fills local proxy datasets"""    
    
    if filePath:
      if not self._filePathExists(filePath):
        self._warning('Experiment data reference %s file "%s" missing' % (PROGRAM, filePath) )
        self.experimentRef = None
        return
      
      if not self._checkFileName(filePath):
        self._warning('%s file name "%s" not valid' % (PROGRAM, filePath) )
        self.experimentRef = None
        return
      
      self.experimentRef = filePath
      
      # runs back-compatibility
      rNuc = self.getExperimentRefNuc('a')
      rNuc.save()
      
      rNuc = self.getExperimentRefNuc()
 
      if rNuc:
        self.updateProxyDataTracks(rNuc, EXTERNAL)
        self._experimentRefNuc = rNuc
       
      else:
        self.experimentRef = None        
        self._experimentRefNuc = None
        return
    
    else:
      self.experimentRef = None
      self._experimentRefNuc = None
    
    self._setAttr(self.root, 'experimentRef', string_(filePath))
    self.root.flush()
    
    return self.experimentRef
     

  def setGenomeRef(self, filePath=None):
    """Fills local proxy datasets"""    
    
    if filePath:
      if not self._filePathExists(filePath):
        self._warning('Genome data reference %s file "%s" missing' % (PROGRAM, filePath) )
        self.genomeRef = None
        return
      
      if not self._checkFileName(filePath):
        self._warning('%s file name "%s" not valid' % (PROGRAM, filePath) )
        self.genomeRef = None
        return
      
      self.genomeRef = filePath
      
      # runs back-compatibility
      rNuc = self.getGenomeRefNuc('a')
      rNuc.save()
      
      rNuc = self.getGenomeRefNuc()
  
      if rNuc:
        self.updateProxyDataTracks(rNuc, EXTERNAL)
        self._genomeRefNuc = rNuc
        
      else:
        self.genomeRef = None        
        self._genomeRefNuc = None
        return
    
    else:
      self.genomeRef = None
      self._genomeRefNuc = None
     
    self._setAttr(self.root, 'genomeRef', string_(filePath))
    self.root.flush()

    return self.genomeRef
   
   
  def getExperimentRefNuc(self, mode='r'):

    if (mode=='r') and self._experimentRefNuc and \
       (self._experimentRefNuc.fileName == self.experimentRef):
      return self._experimentRefNuc
     
    if self.experimentRef and self._filePathExists(self.experimentRef):
      try:
        rNuc = Nucleus(self.experimentRef, experimentRef=False, genomeRef=False, mode=mode)
 
      except Exception() as err:
        msg = 'Reference %s file "%s" could not be loaded.\nPython error:\n%s'
        self._warning(msg % (PROGRAM, self.experimentRef, err))
        return
      
      if mode == 'r': # Only read-only HDF files persist in session
        self._experimentRefNuc = rNuc
      
      return rNuc
  
  
  def getGenomeRefNuc(self, mode='r'):
    
    if (mode=='r') and self._genomeRefNuc and \
       (self._genomeRefNuc.fileName == self.genomeRef):
      return self._genomeRefNuc
    
    if self.genomeRef and self._filePathExists(self.genomeRef):
      try:
        rNuc = Nucleus(self.genomeRef, experimentRef=False, genomeRef=False, mode=mode)
 
      except Exception() as err:
        msg = 'Reference %s file "%s" could not be loaded.\nPython error:\n%s'
        self._warning(msg % (PROGRAM, self.genomeRef, err))
        return
      
      if mode == 'r': # Only read-only HDF files persist in session
        self._genomeRefNuc = rNuc

      return rNuc
    
  
  def _delete(self, name, parent):
  
    if name in parent:
      section = parent[name].name
      self.notifiers.add(section)
      self.isModified = True
      
      del parent[name]
      
  
  def _setAttr(self, obj, name, value):
  
    if name in obj.attrs:
      if isinstance(value, ndarray):
        value = list(value)
         
      if value == array(obj.attrs[name]).tolist():
        return
    
    obj.attrs[name] = value
    
    section = obj.name + '/attrs/' + name
    self.notifiers.add(section)
    self.isModified = True
    
    return True
  
  
  def _getGroup(self, name, parent):
  
    if name in parent:
      return parent[name]
      
    else:
      group = parent.create_group(name)
      section = group.name
      self.notifiers.add(section)
      self.isModified = True
      
      return  group
  

  def _getContactGroup(self, name, parent):
  
    group = self._getGroup(name, parent)
    
    if 'selected' not in group.attrs:
      group.attrs['selected'] = 1
      group.attrs['color'] = (1.0, 1.0, 0.0, 1.0)
    
    return group
  
  
  def _setData(self, name, parent, dataType, dataArray, compression=None):
    
    if name in parent:
      
      #if parent[name].shape != dataArray.shape:
      
      try:
        del parent[name] # Any alternative to this?
      except:
        pass
      
      parent[name] = dataArray
      dataset = parent[name]
      
    else:
      dataset = parent.create_dataset(name, dtype=dataType, data=dataArray,
                                      compression=compression)
    
    section = dataset.name
    self.notifiers.add(section)
    self.isModified = True
    
    return dataset
  
  
  def _warning(self, msg):
  
    print(("* * * * WARNING: %s * * * * " % (msg,)))
   
        
  def _setLinkAttrs(self):
  
    root = self.root
    
    self.attrs = root.attrs
    
    # Convenience group attributes
    
    self.contacts   = root['contacts']
    self.dataTracks = root['dataTracks']
    self.sample     = root['sample']
    self.structures = root['structures']
    self.images     = root['images']
    self.display    = root['display']
    
    self.origContacts = self.contacts['original']
    self.workContacts = self.contacts['working']
    self.externalData = self.dataTracks[EXTERNAL]
    self.derivedData  = self.dataTracks[DERIVED]
    self.innateData   = self.dataTracks[INNATE]
    self.chromosomes  = self.structures['chromosomes']
    self.restraints   = self.structures['restraints']
    
    self.foci = self.images['foci']
    self.grid = self.images['grid']
    
    for typ in self.dataTracks:
      for code in self.dataTracks[typ]:
        attrs = self.dataTracks[typ][code].attrs
        
        if 'options' not in attrs:
          attrs['stranded'] = 0
          attrs['options'] = (0, 0, 1, 0, 0, 0, 0, 0) # shown, labels, symbol, trackType
          attrs['display'] = (0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0) # r, g, b, scale, threshold
    
    for mainGroup in (self.origContacts, self.workContacts):
      for i, groupName in enumerate(mainGroup):
        attrs = mainGroup[groupName].attrs
        
        if 'selected' not in attrs:
          if i == 0:
            selected = 1
            r, g, b = (1.0, 1.0, 0.0)
          elif groupName == 'singleCell':
            selected = 1
            r, g, b = (1.0, 1.0, 0.0)
          else:
            selected = 0  
            r, g, b = (0.0, 1.0, 1.0)
        
          attrs['selected'] = selected
          attrs['color'] = (r, g, b, 1.0)
    
    displayAttrs = self.display.attrs
    
    if 'view3d' not in displayAttrs:
      # Structure    0..2; dx, dy, zoom
      # Contact map  3..5; dx, dy, zoom
      # Interactome  6..8; dx, dy, zoom
      displayAttrs['view3d'] = array([0.0,0.0,1.0,
                                      0.0,0.0,1.0,
                                      0.0,0.0,1.0], float)
                                      
    if 'detailLevel' not in displayAttrs:
      # balls, lines, tube line, tube detail, unused
      displayAttrs['detailLevel'] = array([1,2,1,1], int)
    
    if 'sizes' not in displayAttrs:
      # sphere, bond, line, tube, density, unused...
      displayAttrs['sizes'] = array([1.0,0.5,2.0,1.0,
                                     16.0,1.0,1.0,1.0], float)
   
    if 'options' not in displayAttrs:
      # colorMode, dispMode, cisRestr, transRestr, text labels, scalebar,
      # restr colouring, unused...
      displayAttrs['options'] = array([0, 1, 0, 0,
                                       0, 0, 0, 0], int)

    if 'interactome' not in displayAttrs:
      # displayMode, colorMode, renderMode, thickness, is3d, showCis, showTrans, unused..
      displayAttrs['interactome'] = array([1, 0, 1, 25, 1,
                                           1, 1, 0, 0, 0], int)
    
    if 'models' not in displayAttrs:
      displayAttrs['models'] = [0,]

    if 'colorSchemes' not in self.display:
      schemeGroup = self._getGroup('colorSchemes', self.display)
      schemeAttrs = schemeGroup.attrs
      schemeAttrs['density']    = array([[0.0,0.5,1.0],
                                         [1.0,0.0,0.0],
                                         [1.0,1.0,0.0]], float)
      schemeAttrs['sequence']   = array([[1.0,0.0,0.0],
                                         [1.0,1.0,0.0],
                                         [0.0,1.0,0.0],
                                         [0.0,1.0,1.0],
                                         [0.0,0.0,1.0],
                                         [1.0,0.0,1.0],
                                         ], float)
      schemeAttrs['chromosome'] = array([[1.0,0.0,0.0],
                                         [0.0,1.0,0.0],
                                         [0.0,0.0,1.0]], float)
      schemeAttrs['model']      = array([[1.0,0.0,0.0],
                                         [0.0,1.0,0.0],
                                         [0.0,0.0,1.0]], float)
      schemeAttrs['restraint']  = array([[0.5,0.5,0.5],
                                         [1.0,1.0,0.0],
                                         [1.0,0.0,0.0]], float)
    
    if 'rmsd' not in self.display['colorSchemes'].attrs:
      schemeGroup = self._getGroup('colorSchemes', self.display)
      schemeAttrs = schemeGroup.attrs
      schemeAttrs['rmsd']    = array([[0.5,0.5,0.7],
                                      [1.0,1.0,0.0],
                                      [1.0,0.0,0.0],
                                      [1.0,0.0,1.0]], float)
      
    
    calcGroup = self._getGroup('calculation', self.structures)
    calcAttrs = calcGroup.attrs
    
    for name in STRUC_CALC_DEFAULTS:
      if name not in calcAttrs:
        calcAttrs[name] = STRUC_CALC_DEFAULTS[name]
    
  
  def saveAs(self, filePath):
  
    from shutil import copy2
    
    self.root.flush()
    
    copy2(self.fileName, filePath)

    self.__init__(filePath, self.version, self.experimentRef, self.genomeRef)
  
  
  def getFileName(self):
  
    return self.root.filename
      
  
  def _getDataTrackGroup(self, source=EXTERNAL):
  
    if source.lower() == INNATE:
      return self.innateData
    
    elif source.lower() == DERIVED:
      return self.derivedData
    
    else:
      return self.externalData
      
    
  def _checkChromosome(self, chromosome):
  
    if chromosome not in self.chromosomes:
      raise Exception('Chromosome "%s" not present' % chromosome)

  
  def _checkFileName(self, fileName):
    
    head, tail = os.path.splitext(fileName)
    
    if tail == FILE_EXT:
      return True
    
    elif tail == '.temp' and head.endswith(FILE_EXT):
      return True
      
    else:
      return False
      
      
  def _filePathExists(self, fileName):
  
    if not os.path.exists(fileName):
      return False
    
    if not os.path.isfile(fileName):
      return False
    
    if os.stat(fileName).st_size == 0:
      return False
    
    return True
    
    
  def save(self, fileName=None):
    """Save/SaveAs nucleus data to disk"""
    
    if fileName:
      if self._checkFileName(fileName):
        if fileName != self.fileName:
          self.saveAs(fileName)
        else:
          self.root.flush()
        
      else:
        raise Exception('File "%s" not a valid %s file' % (fileName, PROGRAM))
    
    else:
      self.root.flush()
    
    self.isModified = False
      
    print(("NUC save", fileName))


  def checkValid(self):
    """Determine whether the nucleus file is complete and stored properly"""
    
    if 'contacts' not in self.root:
      return False
    elif 'sample' not in self.root:
      return False
    elif 'dataTracks' not in self.root:
      return False
    elif 'structures' not in self.root:
      return False
    elif 'images' not in self.root:
      return False
    elif 'chromosomes' not in self.root['structures']:
      return False
    
    # Full check required
    
    return True  
  
  
  def getId(self):
    """Get the unique identifying code for the nucleus file"""
    
    sid = self.root.attrs['id']
    return tuple(sid[:2])
  
  
  def getCreationTime(self):
    """Get the time the Nucleus file was created"""
  
    sid = self.root.attrs['id']
    return datetime.fromtimestamp(sid[1])


  def getAccessTime(self):
    """Get the time the Nucleus file was last acceseed"""
  
    sid = self.root.attrs['id']
    return datetime.fromtimestamp(sid[2])
  
  
  def getSampleInfo(self):
    """Get information relating to the experimental sample"""
    
    sample = self.root['sample']
    
    infoDict = {}
    for key, value in sample.attrs.iteritems():
      infoDict[key] = value
    
    return infoDict
  
  
  def getDisplayedChromosomes(self):
         
    chromosomes = [c for c in self.getChromosomes() if self.getChromoDisplayParams(c)[0]]
    
    return chromosomes
  
  
  def getDisplayedModels(self):
    
    models = list(self.display.attrs['models'])
    
    return models
  
  
  def setDisplayedModels(self, models):
    
    return self._setAttr(self.display, 'models', models)
      
  
  def getDisplayedDataTracks(self):
    
    group = self.dataTracks
    keys = []
    
    for typ in group:
      subGroup = group[typ]
      
      for code in subGroup:
        if 'options' not in subGroup[code].attrs:
          continue
                  
        if subGroup[code].attrs['options'][0]:
          refGroup = self.getRefDataTrackGroup(typ, code)
          
          for chromo in refGroup:
            if len(refGroup[chromo]['regions']):
              keys.append( (typ, code) )
              break
    
    return keys
  
  
  def setInteractomeParams(self, displayMode=None, colorMode=None, renderMode=None,
                           thickness=None, is3d=None, showCis=None, showTrans=None):
  
    vals = list(self.display.attrs['interactome'])
    
    for i, param in enumerate([displayMode, colorMode, renderMode,
                               thickness, is3d, showCis, showTrans]):
      if param is not None:
        vals[i] = int(param)
    
    return self._setAttr(self.display, 'interactome', vals)
  
  
  def setDisplayDetailLevels(self, ballDetail=None, lineSmooth=None, tubeSmooth=None, tubeDetail=None):
    
    opts = list(self.display.attrs['detailLevel'])
    
    for i, param in enumerate([ballDetail, lineSmooth, tubeSmooth, tubeDetail]):
      if param is not None:
        opts[i] = max(min(3, int(param)), 0)
    
    return self._setAttr(self.display, 'detailLevel', opts)
    
  
  def setDisplaySizes(self, ball=None, stick=None, line=None, tube=None, density=None):
  
    opts = list(self.display.attrs['sizes'])
   
    for i, param in enumerate([ball, stick, line, tube, density]):
      if param is not None:
        opts[i] = max(float(param), 0.0)
          
    return self._setAttr(self.display, 'sizes', opts)
    
  
  def setRestraintsDisplayed(self, cis=True, trans=True):
  
    attrs = self.display.attrs
    opts = list(attrs['options'])
    
    showCis = 1 if cis else 0
    showTrans = 1 if trans else 0
    
    opts[2] = showCis
    opts[3] = showTrans
    
    return self._setAttr(self.display, 'options', opts)
    
    
  def setTextDisplayed(self, isShown):
  
    attrs = self.display.attrs
    opts = list(attrs['options'])
    opts[4] = 1 if isShown else 0
    
    return self._setAttr(self.display, 'options', opts)


  def setScaleDisplayed(self, isShown):
  
    attrs = self.display.attrs
    opts = list(attrs['options'])
    opts[5] = 1 if isShown else 0
    
    return self._setAttr(self.display, 'options', opts)
  
  
  def setRestraintColoring(self, mode):
  
    attrs = self.display.attrs
    opts = list(attrs['options'])
    opts[6] = mode
    
    return self._setAttr(self.display, 'options', opts)
  
  
  def setModelsDisplayed(self, indices):
  
    return self._setAttr(self.display, 'models', sorted(indices))
    

  def setDataTrackDisplayed(self, typ, code, isShown):
  
     if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']
        val = 1 if isShown else 0
        
        if vals[0] != val:
          vals[0] = val
          self._setAttr(group[code], 'options', vals)
          return True 
          
          
  def setDataTrackLabelled(self, typ, code, isShown):
  
     if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']
        val = 1 if isShown else 0
        
        if vals[1] != val:
          vals[1] = val
          self._setAttr(group[code], 'options', vals)
          return True 


  def setDataTrackColor(self, typ, code, rgba):
    
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['display']
        
        if list(vals[:4]) != list(rgba):
          vals[:4] = rgba
          self._setAttr(group[code], 'display', vals)
          return True
  
  
  def setDataTrackPeakType(self, typ, code, peakType):
    
    if not isinstance(peakType, int):
      peakType = DATA_TRACK_PEAK_TYPES.index(peakType) 
     
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']
        
        if vals[3] != peakType:
          vals[3] = peakType
          self._setAttr(group[code], 'options', vals)
          return True 
  
  
  def getDataTrackPeakType(self, typ, code, asInt=True):
    
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']
        index = min(len(DATA_TRACK_PEAK_TYPES)-1, max(0, vals[3]))
        
        if asInt:
          return index
        else:
          return DATA_TRACK_PEAK_TYPES[index]
          
  
  def getDataTrackColor(self, typ, code, dType=float, alpha=True):
   
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['display']
        
        if dType is float:
          if alpha:
            return list(vals[:4])
          else:
            return list(vals[:3])
          
        elif dType is int:
          if alpha:
            return [int(255*rgb) for rgb in vals[:4]]
          else:
            return [int(255*rgb) for rgb in vals[:3]]
            
        else:
          if alpha:
            vals = [int(255*rgb) for rgb in vals[:4]]
            return '#%02X%02X%02X%02X' % tuple(vals)
            
          else:
            vals = [int(255*rgb) for rgb in vals[:3]]
            return '#%02X%02X%02X' % tuple(vals)
        
  
  def setDataTrackSymbol(self, typ, code, index):
    
    if not isinstance(index, int):
      index = DATA_TRACK_SYMBOLS.index(index)
    
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']
        
        if vals[2] != index:
          vals[2] = index
          self._setAttr(group[code], 'options', vals)
          return True
  
  
  def getDataTrackSymbol(self, typ, code, asInt=True):
       
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        index = group[code].attrs['options'][2]
        
        if asInt:
          return index
        else:  
          return DATA_TRACK_SYMBOLS[index]
   

  def setDataTrackScale(self, typ, code, val):

    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['display']
        
        if vals[4] != val:
          vals[4] = val
          self._setAttr(group[code], 'display', vals)
          return True


  def getDataTrackScale(self, typ, code):

    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        return group[code].attrs['display'][4]
        

  def setDataTrackThreshold(self, typ, code, val):
    
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['display']
        
        if vals[5] != val:
          vals[5] = val
          self._setAttr(group[code], 'display', vals)    
          return True


  def getDataTrackThreshold(self, typ, code):
    
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        return group[code].attrs['display'][5]
           
   
  def areContactsSingleCell(self, groupName):
  
    group = self.getContactGroup(groupName)
    
    if group and 'isSingleCell' in group.attrs:
      return group.attrs['isSingleCell'] == 1
      
    return True


  def getContactsBinSize(self, groupName):
  
    group = self.getContactGroup(groupName)
    
    if group and 'binSize' in group.attrs:
      return group.attrs['binSize']
    
    
  def getContactGroup(self, groupName, working=None):
    
    if working is None:
      if groupName in self.workContacts:
        return self.workContacts[groupName]
        
      elif groupName in self.origContacts:
        return self.origContacts[groupName]
   
    elif working:
      if groupName in self.workContacts:
        return self.workContacts[groupName]
    
    else:
      if groupName in self.origContacts:
        return self.origContacts[groupName]
  
  
  def getContactGroups(self, working=None):
    
    if working is None:
      groups = [(n, self.workContacts[n]) for n in self.workContacts]
        
      for groupName in self.origContacts:
        if groupName not in self.workContacts:
          groups.append((groupName, self.origContacts[groupName]) )
   
    elif working:
      groups = [(n, self.workContacts[n]) for n in self.workContacts]
    
    else:
      groups = [(n, self.origContacts[n]) for n in self.origContacts]
        
    return groups
     
  
  def getContactGroupColor(self, groupName, dType=float, alpha=True):
    
    group = self.getContactGroup(groupName)
    
    if group:
      vals = group.attrs['color']
 
      if dType is float:
        if alpha:
          return vals[:4]
        else:
          return vals[:3]
 
      elif dType is int:
        if alpha:
          return [int(255*rgb) for rgb in vals[:4]]
        else:
          return [int(255*rgb) for rgb in vals[:3]]
 
      else:
        if alpha:
          vals = [int(255*rgb) for rgb in vals[:4]]
          return '#%02X%02X%02X%02X' % tuple(vals)
 
        else:
          vals = [int(255*rgb) for rgb in vals[:3]]
          return '#%02X%02X%02X' % tuple(vals)
    
  
  def getSelectedContactGroups(self):
  
    names = [name for name, group in self.getContactGroups() if group.attrs['selected']]
    
    return names
    
    
  def setContactGroupColor(self, groupName, rgb):
  
    group = self.getContactGroup(groupName)
    
    if group:
      color = list(group.attrs['color'])
      color[:len(rgb)] = rgb
      return self._setAttr(group, 'color', rgb)    
      
 
  def setSelectedContactGroups(self, groupNames):  
    
    groupNames = set(groupNames)
     
    for name, group in self.getContactGroups():
      if name in groupNames:
        self._setAttr(group, 'selected', 1)
      else:
        self._setAttr(group, 'selected', 0)
        
        
  def setChromoGroupSelected(self, groupName, isSelected):
  
    group = self.getContactGroup(groupName)
    
    if group:
      val = 1 if isSelected else 0
      self._setAttr(group, 'selected', val)
      
    
  def getChromoDisplayParams(self, name, default=(1,1,1,1)):
  
    if name in self.chromosomes:
      if 'display' in self.chromosomes[name].attrs:
        isShown, useLabels, colMode, dispMode = self.chromosomes[name].attrs['display']
        
      else:
        self.chromosomes[name].attrs['display'] = array(default, int)
        isShown, useLabels, colMode, dispMode = default
      
    else:
      isShown, useLabels, colMode, dispMode = default
    
    return isShown, useLabels, colMode, dispMode
  
  
  #set chromosome section centre (relative, between 0.0 and 1.0) and the section length to be displayed
  #useful for displaying only a section of chromosomes
  def getChromoSectionDisplayParams(self, name, section_default=[0.5,0]):
  
    if name in self.chromosomes:
      if 'display_section' in self.chromosomes[name].attrs:
        sectionDisplay = self.chromosomes[name].attrs['display_section']
        
      else:
        self.chromosomes[name].attrs['display_section'] = section_default
        sectionDisplay = section_default
      
    else:
      sectionDisplay = section_default
    
    return sectionDisplay
  
  
  #set chromosome section centre (relative, between 0.0 and 1.0) and the section length to be displayed
  #useful for displaying only a section of chromosomes
  def setChromoSectionDisplayParams(self, chromo, section_centre=None, section_numres=None):

     #print "setting display params for chromosome ", chromo, "centre=", section_centre, "length=", section_numres  
     if chromo in self.chromosomes:
       opts = self.chromosomes[chromo].attrs['display_section']
       #print "stored section display parameters: ", opts

       if section_centre is None:
          centre = opts[0]
       else:
          centre = section_centre
       if section_numres is None:
          numres = opts[1]
       else:
          numres = section_numres

       #print "setting display params for chromosome ", chromo, "centre=", centre, "length=", numres  
       return self._setAttr(self.chromosomes[chromo], 'display_section', [centre, numres])
  
       
  def resetChromoColors(self):
    
    chromos = self.getChromosomes()
    n = float(len(chromos))
    h = -2/n
    
    for i, name in enumerate(chromos):
      
      if i % 2 == 1:
        s = 0.6
        v = 1.0
        
      else:
        s = 1.0
        v = 0.8
        h += 2/n
      
      r, g, b = colorsys.hsv_to_rgb(h, s, v)
      color = array([r, g, b, 1.0], float)
      self._setAttr(self.chromosomes[name], 'color', color)
   
   
  def setChromoColor(self, chromo, rgba):
    
    if chromo in self.chromosomes:
      return self._setAttr(self.chromosomes[chromo], 'color', rgba) 
    

  def getChromoColor(self, name, default=(1.0, 0.0, 0.0, 1.0), dType=float, alpha=True):
  
    if name in self.chromosomes:
      if 'color' in self.chromosomes[name].attrs:
        color = array(self.chromosomes[name].attrs['color'])
        
      else:
        self.resetChromoColors()
        color = array(self.chromosomes[name].attrs['color'])
      
    else:
      color = array(default, float)
      
    if dType is float:
      if alpha:
        return color
      else:
        return color[:3]
 
    elif dType is int:
      if alpha:
        return [int(255*rgb) for rgb in color]
      else:
        return [int(255*rgb) for rgb in color[:3]]
 
    else:
      if alpha:
        vals = [int(255*rgb) for rgb in color]
        return '#%02X%02X%02X%02X' % tuple(vals)
 
      else:
        vals = [int(255*rgb) for rgb in color[:3]]
        return '#%02X%02X%02X' % tuple(vals)
       
    
  def setChromoDisplayParams(self, chromo, isShown=None, useLabels=None, colorMode=None, displayMode=None):
  
     if chromo in self.chromosomes:
       #TODO set default values when 'display' doesn't exist
       opts = self.chromosomes[chromo].attrs['display']
  
       for i, param in enumerate([isShown, useLabels, colorMode, displayMode]):
         if param is None:
           continue
         
         if param is True:
           param = 1
         elif param is False:
           param = 0
            
         opts[i] = int(param)
         
       return self._setAttr(self.chromosomes[chromo], 'display', opts)
  
       
  def getChromoDisplayed(self, chromo):
      
    if chromo in self.chromosomes:
      return bool(self.chromosomes[chromo].attrs['display'][0])
  
        
  def setChromoDisplayed(self, chromo, isShown):
      
    if chromo in self.chromosomes:
      opts = self.chromosomes[chromo].attrs['display']
      opts[0] = 1 if isShown else 0
    
      return self._setAttr(self.chromosomes[chromo], 'display', opts)
  
  
  def setColorScheme(self, name, colors):
  
    schemeGroup = self._getGroup('colorSchemes', self.display)
    schemeAttrs = schemeGroup.attrs  
    schemeAttrs[name] = array(colors, float)
    
      
  def getColorScheme(self, name):    
      
    schemeGroup = self._getGroup('colorSchemes', self.display)
    schemeAttrs = schemeGroup.attrs
    
    if name in schemeAttrs:
      return array(schemeAttrs[name])
  
  
  def calcInterpolatedColors(self, schemeColors, nPoints, alpha=1.0):
    
    # Redo in cython?
     
    nBlock = len(schemeColors)-1
    blockSize = float(nPoints-1)/nBlock
    
    if alpha is None:
      colorArray = empty((nPoints, 3))
      colorArray[-1] = schemeColors[-1]
    else:
      colorArray = empty((nPoints, 4))
      colorArray[:,3] = alpha
      colorArray[-1,:3] = schemeColors[-1]
    
    if not blockSize:
      return colorArray
    
    for i in range(nBlock):
      p1 = int(i * blockSize)
      p2 = int((i+1) * blockSize)
      delta = p2-p1
      ramp = array(range(delta), float) 
      f = ramp / delta
      g = 1.0 - f
      color = outer(g, schemeColors[i]) + outer(f, schemeColors[i+1])
      colorArray[p1:p2,:3] = color
    
    return colorArray
        
  def getChromosomes(self):
    """Get a list of chromosome or DNA molecule representations"""
    
    if self._chromosomes is None:
      names = [c for c in self.chromosomes]
      names.sort(key=naturalKey())
      self._chromosomes = names
      return names
    
    else:
      return self._chromosomes

  
  def getChromosomeLimits(self, chromosome):
    """Get the first and last sequence positions for a chromsome"""
  
    chromoGroup = self.chromosomes
    
    if chromosome in self._chromoLimitCache:
      start, end = self._chromoLimitCache[chromosome]
    
    else:
      positions = chromoGroup[chromosome]['positions']
      start = positions[0]
      end = positions[-1]
      self._chromoLimitCache[chromosome] = (start, end)
 
    
    return start, end
    
    
  def getDomains(self, chromosome):
    """Get a list of compated topological domain regions"""
    
    if 'domains' in self.derivedData and chromosome in self.derivedData['domains']:
      return array(self.derivedData['domains'][chromosome])
    
  
  def getParticlePositions(self, chromosomes=None, cis=True, trans=True,
                           backbone=None):
    """Get the chromosomal sequence locations of model particles"""
   
    chromoGroup = self.chromosomes
    posDict = {}
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    for chromo in chromosomes:
      if chromo not in chromoGroup:
        continue
    
      positions = array(chromoGroup[chromo]['positions'])
      selection = ones(len(positions), uint32)
      
      if backbone is True:
        bools =  array(chromoGroup[chromo]['backbone'])
        selection *= bools
      
      elif backbone is False:
        bools =  1 - array(chromoGroup[chromo]['backbone'])
        selection *= bools
      
      if not cis and (chromo in self.restraints):
        indicesA = set(chromo[chromo][:,0])
        indicesB = set(chromo[chromo][:,1])
        indices = set(indicesA) | set(indicesB)
        selection[tuple(indices)] = 0
      
      if not trans:
        for chromoA in self.restraints:
          subGroup = self.restraints[chromoA]
 
          for chromoB in subGroup:
            if chromoA == chromoB:
              continue
 
            if chromo in (chromoA, chromoB):
              indicesA = set(subGroup[chromoB][:,0])
              indicesB = set(subGroup[chromoB][:,1])
              indices = set(indicesA) | set(indicesB)
              selection[tuple(indices)] = 0
      
      if selection is not None:
        positions = positions[selection.nonzero()]
            
      posDict[chromo] = positions
    
    return posDict
  
  
  def getNumModels(self):
    """Get number of 3D coordinate structural models"""
    
    if 'coords' in self.structures:
      coords = self.structures['coords']
      for chromo in coords:
         return coords[chromo].shape[0]
      
    return 0


  def getNumCoords(self):
    """Get number of 3D coordinates in each structural model"""
    
    if 'coords' in self.structures:
      coords = self.structures['coords']
      for chromo in coords:
         return coords[chromo].shape[1]
      
    return 0
  
  
  def getCachedContacts(self, groupName):
    
    # TBD forget at some point to save memory
    
    if groupName in self._contactsCache:
      return self._contactsCache[groupName]
    
    if groupName in self.workContacts:
      group = self.workContacts[groupName]
    
    elif groupName in self.origContacts:
      group = self.origContacts[groupName]
    
    else:
      return {}
    
    subDict = {}
    
    for chrA in group:
      subDict[chrA] = {}
      
      for chrB in group[chrA]:
        data = group[chrA][chrB]
        
        if data.shape[1]:
          subDict[chrA][chrB] = array(data, uint32)
          
    self._contactsCache[groupName] = subDict

    return subDict 
          
  
  def getNumContacts(self, groupName, chromosomes=None, cis=True,
                     trans=True, asFraction=False):
    """Get the number of chromosomal contacts (without loading them)"""
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    chromosomes = set(chromosomes)
    
    subGroup = self.getCachedContacts(groupName)
    
    if not subGroup:
      return 0

    n = 0
    t = 0
    
    for chrA in subGroup:
      for chrB in subGroup[chrA]:
        c = subGroup[chrA][chrB].shape[1]
        t += c
        
        if ((chrA == chrB) and cis) or ((chrA != chrB) and trans):
          if (chrA in chromosomes) or (chrB in chromosomes):
            n += c

    if asFraction:
      if t:
        return n/float(t)
    
    else:
      return n
  
  
  def getNumRestraints(self, chromosomes=None, cis=True,
                       trans=True, asFraction=False):
    """Get number of distance restraints"""

    allChromosomes = self.getChromosomes()
    
    if not chromosomes:
      chromosomes = allChromosomes
    
    chromosomes = set(chromosomes)
      
    n = 0
    t = 0
    
    for chromoA in self.restraints:
      subGroup = self.restraints[chromoA]
      
      for chromoB in subGroup:
        c = subGroup[chromoB].shape[1]
        t += c
        
        if (chromoA == chromoB) and cis:
          if chromoA in chromosomes:
            n += c
        
        elif (chromoA != chromoB) and trans:
          if (chromoA in chromosomes) or (chromoB in chromosomes):
            n += c 
    
    if asFraction:
      if t:
        return n/float(t)
    
    else:
      return n
    
  def getRestraints(self, chromosomes=None, cis=True, trans=True, usePositions=False,
                    backboneSep=None, bboneLower=0.8, bboneUpper=1.2):
    """Get a list of stored distance restraints.
       Backbone sep is distance for a base pair"""
        
    # Filter on number of observations?
    
    restraintDict = {}
    posCache = {}
    chromoGroup = self.chromosomes
      
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    chromosomes = set(chromosomes)
    
    # IndexA (uint32) IndexB (uint32) NumObs (uint32)
    # NumObs (float) Weight (uint8) [ Target (float) Lower (float) Upper (float) ]]
    
       
    for chromoA in self.restraints:
      if (chromoA not in chromoGroup) or (chromoA not in chromosomes):
        continue
      
      subGroup = self.restraints[chromoA]
      if chromoA in posCache:
        positionsA = posCache[chromoA]
        
      else:
        positionsA = array(chromoGroup[chromoA]['positions'], uint32)
        
        posCache[chromoA] = positionsA
       
      for chromoB in subGroup:
        if chromoA == chromoB:
          if not cis:
           continue

        elif not trans:
          continue
      
        if (chromoB not in chromoGroup) or (chromoB not in chromosomes):
          continue
          
        if usePositions:
          if chromoB in posCache:
            positionsB = posCache[chromoB]
            
          else:
            positionsB = array(chromoGroup[chromoB]['positions'], uint32)
            
          posCache[chromoB] = positionsB
 
        pointsA = array(subGroup[chromoB][0], uint32)
        pointsB = array(subGroup[chromoB][1], uint32)
     
        if not len(pointsA):
          continue      
        
        if usePositions:
          posA = positionsA[pointsA]
          posB = positionsB[pointsB]
          
        else:
          posA = pointsA
          posB = pointsB
                  
        restAB = subGroup[chromoB][2:]
        dataTrans = vstack([posA, posB, restAB])
        restraintDict[(chromoA, chromoB)] = dataTrans.T
              
                
    if backboneSep:
      for chromo in chromosomes:
        if chromo not in chromoGroup:
          continue
          
        if chromo in posCache:
          positions = posCache[chromo]
 
        else:
          positions = array(chromoGroup[chromo]['positions'], uint32)

        numRes = len(positions)-1
        pos0 = positions[:-1]
        pos1 = positions[1:]
        #numObs = zeros(numRes, float)
        weights = ones(numRes, float32)
        values = backboneSep * ones(numRes, float32)
        
        seqSeps = pos1-pos0

        values *= seqSeps
        lowers = values * bboneLower
        uppers = values * bboneUpper
        
        if not usePositions:
          pos0 = array(range(numRes), uint32)
          pos1 = pos0 + 1
        
        data = vstack([pos0, pos1, weights, values, lowers, uppers])
        restraintDict[chromo] = data.T      
    
    
    return restraintDict
  
  
  def getMapability(self, chromosome):
    """Get a list of sequence mapability for chromosomal particles"""
    
    chromoGroup = self.chromosomes
    self._checkChromosome(chromosome)
      
    mapability = array(chromoGroup[chromosome]['mapability'])
    positions = array(chromoGroup[chromosome]['positions'])
    
    return vstack([positions, mapability]).T
   
  
  def getContactGroupNames(self, working=None):
    """Get a list of all available contact groups
       Original, working or all available"""
  
    if working is None:
      groups = [self.workContacts, self.origContacts]
   
    elif working:
      groups = [self.workContacts,]
    
    else:
      groups = [self.origContacts,]
  
    names = set() # Working often mirror original
    for mainGroup in groups:
      names.update( [name for name in mainGroup]  )
 
    return sorted(names)
    
   
  def getContacts(self, groupName='singleCell', chromosomes=None,
                  cis=True, trans=True, model=None):
    """Get list of interacting chromosomal positions"""  

    contactDict = {}
    cacheDict = self.getCachedContacts(groupName)
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
        
    chromosomes = set(chromosomes)     
    
    if cis:
      for chromo in cacheDict:
        if chromo not in chromosomes:
          continue
        
        if chromo in cacheDict[chromo]:
          contacts = cacheDict[chromo][chromo]
 
          if len(contacts):
            contacts = array(contacts, uint32)
            
            # TBD Cython this filtering for speed
            if (model is not None) and (contacts.shape[0] == 4):
              indices = (contacts[3] == model).nonzero()[0]
              contacts = contacts[:,indices]
              
            contactDict[(chromo, chromo)] = contacts
 
    
    if trans:
      for chromoA in cacheDict:
        if chromoA not in chromosomes:
          continue

        subDict = cacheDict[chromoA]
        
        for chromoB in subDict:
          if chromoB not in chromosomes:
            continue
          
          if chromoA == chromoB:
            continue
            
          contacts = subDict[chromoB]
        
          if len(contacts):
            contacts = array(contacts, uint32)
           
            if (model is not None) and (contacts.shape[0] == 4):
              indices = (contacts[3] == model).nonzero()[0]
              contacts = contacts[:,indices]
              
            contactDict[(chromoA, chromoB)] = contacts
  
    return contactDict
   
   
  def getCloseContacts(self, groupName, position, threshold=int(1e6), cis=True, trans=True):
    """Get list of chromosomal contacts and their observations
       that are close to a given position"""  

    chromo, pos = position
    contactData = []

    cacheDict = self.getCachedContacts(groupName)
          
    if cis and (chromo in cacheDict) and (chromo in cacheDict[chromo]):
      data = array(cacheDict[chromo][chromo])
    
      pos0 = data[0]
      pos1 = data[1]
      numObs = data[2]
      
      deltas0  = abs(pos0 - pos)
      deltas1  = abs(pos1 - pos)
               
      indices0 = argwhere(deltas0 < threshold)
      indices1 = argwhere(deltas1 < threshold)
      
      # removed pos0[i][0] below...
      contactData += [(chromo, pos0[i], numObs[i]) for i in indices0]
      contactData += [(chromo, pos1[i], numObs[i]) for i in indices1]
    
    if trans:
      for chromoA in cacheDict:
        for chromoB in cacheDict[chromoA]:
          if chromoA == chromoB:
            continue
        
          data = array(cacheDict[chromoA][chromoA])
          
          posA = data[0]
          numObsA = data[2]
          deltasA = abs(posA - pos)
          indicesA = argwhere(deltasA < threshold)
          
          # removed pos0[i][0] below...
          contactData += [(chromo, posA[i], numObsA[i]) for i in indicesA]         
   
    return contactData


  def getShuffledContacts(self, groupName, chromosomes=None, cis=True, trans=True, exponent=-1.0):
    """Get list of shuffled interacting chromosomal positions"""  
    
    # Calculations will tend to use this on-the-fly

    if not chromosomes:
      chromosomes = self.getChromosomes()
      
    cacheDict = self.getCachedContacts(groupName)
    contactDict = {}
    
    if cis:
      for chromo in cacheDict:
        if chromo in cacheDict[chromo]:
          data = cacheDict[chromo][chromo]
          pos0 = data[0]
          pos1 = data[1]
          numObsA = data[2]
 
          pool = list(pos1)
          array(shuffle(pool))
          posShuff = []
 
          # Redo in Cython
 
          j = 0
          while pool:
            probs = pool - pos0[j]
            probs[(probs < 1).nonzero()] = 1.0
            probs = power(probs, exponent)
            point = uniform(0, sum(probs))
            csum = cumsum(probs)
            index = (csum >= point).nonzero()[0][0]
 
            pos = pool.pop(index)
            posShuff.append(pos)
            j += 1
 
          dataCis = vstack([pos0, array(posShuff), numObsA])
          contactDict[(chromo, chromo)] = dataCis.T
    
    if trans:
      for chromoA in cacheDict:
        subDict =  cacheDict[chromoA]
        
        for chromoB in subDict:
          data = subDict[chromoB]
          posA = data[0]
          posB = data[1]
          
          random.shuffle(posB) 
          numObsAB = data[2]
                          
          dataTrans = vstack([posA, posB, numObsAB])
          contactDict[(chromoA, chromoB)] = dataTrans.T
    
    return contactDict
   
   
  def getContactMatrix(self, chrA, chrB, binSize=int(1e6),
                      groupName='singleCell', model=None):
    """Get full grid of contacts"""  
    
    from cUtil.apiUtil import binContacts
    
    self._checkChromosome(chrA)
    self._checkChromosome(chrB)
    
    chromoGroup = self.chromosomes
    
    if chrA == chrB:
      cis = True
      trans = False
      chromoKeys = [(chrA, chrB),] 
      chromos = (chrA,)
    else:
      trans = True
      cis = False
      chromoKeys = [(chrA, chrB), (chrB, chrA)]  
      chromos = (chrA, chrB)
    
    contactDict = self.getContacts(groupName, chromos, cis, trans, model)
   
    # # # # Might need to make this work off the contact limits
    startA, endA = self.getChromosomeLimits(chrA)
    startB, endB = self.getChromosomeLimits(chrB)
    
    extentA = endA - startA
    extentB = endB - startB
    
    n = int(ceil(extentA/binSize)) + 1
    m = int(ceil(extentB/binSize)) + 1
    
    matrix = zeros((n,m), int32)
    binSize = int32(binSize)
    
    for chromoKey in chromoKeys:
      if chromoKey in contactDict:
        chr1, chr2 = chromoKey
        cData = array(contactDict[chromoKey], int32)
        
        if chrA == chr1:
          start1, start2 = int32(startA), int32(startB)
          transpose = False
        else:
          start1, start2 = int32(startB), int32(startA)
          transpose = True
        
        binContacts(cData, matrix, start1, start2,
                    binSize, cis, transpose)
    
    return matrix


  def getDistanceMatrix(self, chromosome, models=None, binSize=int(1e6), useMinVal=False):
    """Get full grid of model distances along the chromosome sequence.
       Averages over all models or selected models."""  
    
    t0 = time.time()
    
    self._checkChromosome(chromosome)
    
    strucGroup = self.structures
    
    if 'coords' in strucGroup:
      coords = strucGroup['coords']
      
      if chromosome not in coords:
        return
      
      chromoCoords = coords[chromosome]
      chromoGroup = self.chromosomes
      
      start, end = self.getChromosomeLimits(chromosome)
      extent = end - start
      
      n = int(extent//binSize) + 1 
      m = self.getNumModels()
      
      if models:
        models = [i for i in models if i < m]
      else:
        models = range(m)
      
      # Do quicker in Cython
      
      positions = array(chromoGroup[chromosome]['positions'])
      bins = array((positions-start)/binSize, int)
      nPos = len(positions)
      coords = array(chromoCoords[models])
      matrix = zeros((m, nPos,nPos), float)
      
      for i in range(nPos):
        for k in models:
          deltas = coords[k,i]-coords[k]
          matrix[k,i] = sqrt((deltas*deltas).sum(axis=1))
      
      if useMinVal: # Min value over models
        matrix = matrix.min(axis=0)
        
      else: # Mean value over models  
        matrix = matrix.mean(axis=0)
      
      distMatrix = zeros((n,n), float)
      if useMinVal:
        for i, a in enumerate(bins):
          for j, b in enumerate(bins[i:], i):
            val = matrix[i,j]
            prev = distMatrix[a,b]
            
            if prev:
              if val < prev:
                distMatrix[a,b] = val
            
            else: 
              distMatrix[a,b] = val
      
      else:
        counts = zeros((n,n), float)
        for i, a in enumerate(bins):
          for j, b in enumerate(bins[i:], i):
            counts[a,b] += 1.0
            distMatrix[a,b] += matrix[i,j]
        
        indices = counts.nonzero()
        distMatrix[indices] /= counts[indices]
      
      distMatrix += distMatrix.T
        
      return distMatrix
        
        
  def getCisDistances(self, position, models=None, step=int(1e6)):
    """Get list of model distances to a specified point along the chromosome sequence"""
    
    chromo, pos = position
    self._checkChromosome(chromo)
    
    m = self.getNumModels()
    if models:
      models = [i for i in models if i < m]
    else:
      models = range(m)
    
    strucGroup = self.structures
    
    if 'coords' in strucGroup:
      coords = strucGroup['coords']
      
      if chromo in coords:
        chromoGroup = self.chromosomes
        start, end = self.getChromosomeLimits(chromo)
        extent = end - start
      
        n = int(extent//step) + 1
        posB = [(chromo, step * i) for i in range(n)]
        posA = [position] * len(posB)
        
        return self.getPositionDistances(posA, posB, models)
  
  
  def getModelSize(self, model, chromosomes=None):
    """Get the size of chromosomal model, in microns, along its longest and shortest orthogonal axes"""
    
    from util.Cluster import principleComponentAnalysis
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    strucGroup = self.structures
    if 'coords' in strucGroup:
      allCoords = self.getModelCoords(model, chromosomes)
      basis, energy = principleComponentAnalysis(allCoords)
      pc1 = basis[:,0]
      pc2 = basis[:,1]
      pc3 = basis[:,2]
      origin = allCoords.mean(axis=0)
      
      allCoords -= origin
      
      proj1 = dot(allCoords, pc1) # might have to transpose
      proj2 = dot(allCoords, pc3)
      
      extent1 = proj1.max(axis=0) - proj1.min(axis=0)
      extent2 = proj2.max(axis=0) - proj2.min(axis=0)
      
      return extent1, extent2
      
      
  def getGlobalTransform(self):
    """Get affine transformation associated with all structure models"""
    
    tformGroup = self._getGroup('transforms', self.structures)
    
    if 'global' in tformGroup:
      return array(tformGroup['global'])
    
    else:
      return eye(4)  
  
  
  def getModelTransform(self, model):
    """Get affine transformation associated with structure model"""
    
    tformGroup = self._getGroup('transforms', self.structures)
    
    if 'models' in tformGroup:
      if model < len(tformGroup['models']):
        return array(tformGroup['models'][model])
    
    
  def getModelCoords(self, model, chromosomes=None, backbone=None):
    """Retrieve the 3D structure model coordinates of one or more chromosomes.
       Backbone True:only, False:exclude"""

    if chromosomes:
      chromosomes = sorted(chromosomes)
    else:
      chromosomes = self.getChromosomes()
    
    
    chromoGroup = self.chromosomes
    strucGroup = self.structures
    
    if 'coords' in strucGroup:
      coordsGroup = strucGroup['coords']
      
      allCoords = None
      for chromo in chromosomes:
        if chromo in coordsGroup:
          chromoCoords = coordsGroup[chromo]
          nModels = len(chromoCoords)
          coords = array(chromoCoords[min(nModels-1,model)])
          
          if backbone is not None:
            if backbone:
              indices = array(chromoGroup['backbone']).nonzero()
              coords = coords[indices]
            
            else:
              indices = (chromoGroup['backbone'] < 0).nonzero()
              coords = coords[indices]
          
          if len(coords):
            if allCoords is None:
              allCoords = coords
            else:
              allCoords = append(allCoords, coords, axis=0)
      
      return allCoords
      
    
  def getModelBoundingBox(self, model, chromosomes=None):
    """Get the maximum extent of the structural model along the coordinate axes"""

    allCoords = self.getModelCoords(model, chromosomes)
    
    minPos = allCoords.min(axis=0).tolist()
    maxPos = allCoords.max(axis=0).tolist()
    
    return zip(minPos, maxPos)
    
    
  def getPositionCoords(self, model, positions, chromo=None):
    """Get the (interpolated) 3D coordinates for a list of sequential positions"""

    from cUtil.apiUtil import getInterpolatedCoords
    
    strucGroup = self.structures
    if 'coords' in strucGroup:
      coordsGroup = strucGroup['coords']
      chromoGroup = self.chromosomes
      
      
      if chromo is None:
        posDict = {}
        coordDict = {}
        coordData = []
        coordDataAppend = coordData.append
        
        for chromo, pos in positions:
          if chromo in posDict:
            coords = coordDict[chromo]
            points = posDict[chromo]
 
          elif chromo in coordsGroup:
            coords = coordDict[chromo] = array(coordsGroup[chromo][model])
            points = posDict[chromo] = array(chromoGroup[chromo]['positions'], int32)
 
          else:
            coordData.append(None)
            continue
 
          p = getInterpolatedCoords(points, coords, array([pos], int32))
          coordDataAppend(p[0])
        
        coordData = array(coordData)
        
      else:
        coords = array(coordsGroup[chromo][model])
        points = array(chromoGroup[chromo]['positions'], int32)
        pos = array(positions, int32)
        coordData = getInterpolatedCoords(points, coords, pos)
      
      return coordData
      
      
  def getPositionDistances(self, positionsA, positionsB, models=None, chromoA=None, chromoB=None):
    """Get the 3D model distances between a given sets of sequence positions"""
    
    nModels = self.getNumModels()
    if models:
      models = [i for i in models if i < nModels ]
    else:
      models = range(nModels )
    
    nA = len(positionsA)
    nB = len(positionsB)
    
    if nA > nB:
      positionsA = positionsA[:nB]
    elif nB > nA:
      positionsB = positionsB[:nA]  
    
    modelDists = []  
    for m in models:
      coordData1 = self.getPositionCoords(m, positionsA, chromoA)
      coordData2 = self.getPositionCoords(m, positionsB, chromoB)
      
      deltas = coordData1 - coordData2
      dists = sqrt((deltas*deltas).sum(axis=1))
      modelDists.append(dists)
  
    meanDists = array(modelDists).T.mean(axis=1)
    
    return meanDists
  
  
  def getCoordDistances(self, chromoA, chromoB, models=None):
    """Get the 3D model distances between particle coords of a given pair of chromosomes"""
    
    coordsGroup = self.structures['coords']
    
    if chromoA not in coordsGroup:
      return
      
    if chromoB not in coordsGroup:
      return  
    
    nModels = self.getNumModels()
    if models:
      models = [i for i in models if i < nModels ]
    else:
      models = range(nModels )
    
    coordsA = array(coordsGroup[chromoA])
    coordsB = array(coordsGroup[chromoB])
    
    modelDists = []  
    for m in models:
      coordData1 = coordsA[m]
      coordData2 = coordsB[m]
      
      deltas = coordData1 - coordData2
      dists = sqrt((deltas*deltas).sum(axis=1))
      modelDists.append(dists)
  
    meanDists = array(modelDists).T.mean(axis=1)
    
    return meanDists
    
  
  def getRandomCoords(self, model, chromosome, num):
    """Get 3D coordinates for random sequence positions along a backbone path"""
    
    self._checkChromosome(chromosome)
    
    chromoGroup = self.chromosomes
    start = chromoGroup[chromosome]['positions'][0]
    end = chromoGroup[chromosome]['positions'][-1]
    
    positions = random.randint(start, end, num)
    
    return self.getPositionCoords(model, positions, chromosome)
  
  
  def getBackboneSpacing(self):
    """Get the maximum sequence separation between backbone particles"""
  
    chromosomes = self.getChromosomes()
    
    if chromosomes:
      chromoGroup = self.chromosomes
              
      for chromo in chromosomes:
        if chromo in chromoGroup:
          positions = array(chromoGroup[chromo]['positions'])
          seps = positions[1:] - positions[:-1]
          return seps.max()
          
  
  def calcNormDataTrack(self, code, source=EXTERNAL, chromosomes=None, method='max'):
    """Normalise data track values, preserving any original data.
       Allowed methods are 'untity', 'max', 'orig', 'log', 'quantile' """    
             
    dataGroup = self._getDataTrackGroup(source)

    if code in dataGroup:
      dataLayer = dataGroup[code]
      refNuc = self.getRefDataTrackNuc(source, code, 'a')
      refLayer = refNuc._getDataTrackGroup(source)[code]
      
      if not chromosomes:
        chromosomes = self.getChromosomes()

      if not chromosomes:
        chromosomes = dataLayer.keys()
      
      maxVals = []
      
      for chromo in chromosomes:
        if chromo not in refLayer:
          continue
        
        subGroup = refLayer[chromo]
        valueData  = array(subGroup['values']).T # origValue, normValue
        
        if not len(valueData):
          continue
        
        origValues = valueData[0]
                
        if method == 'unity':
          values = ones(len(origValues), float)
          maxVals.append(1.0)
          
        elif method == 'max':
          values = origValues - origValues.min()
          values /= values.max() or 1.0
          maxVals.append(1.0)

        elif method == 'sqrt':
          values = origValues - origValues.min()
          values /= values.max() or 1.0
          values = sqrt(values)
          maxVals.append(1.0)
         
        elif method == 'orig':
          values = origValues
          maxVals.append(values.max())

        elif method == 'log':
          values = origValues - origValues.min()
          values = log(1.0 + values)
          vMax = values.max() or 1.0
          values /= vMax
          maxVals.append(1.0)
          
        elif method == 'quantile':
          order = origValues.argsort()
          step = 1.0/len(origValues)
          refValues = arange(0.0, 1.0+step, step)
          values = refValues[order.argsort()]
          maxVals.append(1.0)
  
        valueData[1] = values
        refNuc._setData('values', subGroup, float, valueData.T)
      
      dispVals = list(dataLayer.attrs['display'])
      dispVals[7] = max(maxVals or [0.0])
      dataLayer.attrs['display'] = dispVals
      
      key = (source, code, chromo) 
      if key in self._dataTrackHistogramCache:
        del self._dataTrackHistogramCache[key]            
      
      if refNuc is not self:
        refNuc.save()      

  
  def getRefDataTrackGroup(self, source, code, mode='r'):    
    
    dataGroup = self._getDataTrackGroup(source)
    if code in dataGroup:
      group = dataGroup[code]
      for chromo in group:
        if len(group[chromo]['regions']):
          return group
    
    if source == EXTERNAL:
      tNuc = self.getExperimentRefNuc(mode)
      if tNuc:
        dataGroup = tNuc._getDataTrackGroup(source)
        if code in dataGroup:
          return dataGroup[code]
    
    if source == INNATE:
      gNuc = self.getGenomeRefNuc(mode)
      if gNuc:
        dataGroup = gNuc._getDataTrackGroup(source)
        if code in dataGroup:
          return dataGroup[code]
    
    # Empty, but no ref
    if code in dataGroup:
      return dataGroup[code]
   
    
  def getRefDataTrackNuc(self, source, code, mode='r'):
  
    dataGroup = self._getDataTrackGroup(source)
    if code in dataGroup:
      group = dataGroup[code]
      for chromo in group:
        if len(group[chromo]['regions']):
          return self
    
    if source == EXTERNAL:
      tNuc = self.getExperimentRefNuc(mode)
      if tNuc:
        dataGroup = tNuc._getDataTrackGroup(source)
        if code in dataGroup:
          return tNuc
    
    if source == INNATE:
      gNuc = self.getGenomeRefNuc(mode)
      if gNuc:
        dataGroup = gNuc._getDataTrackGroup(source)
        if code in dataGroup:
          return gNuc
          
    # Empty, but no ref
    return self
  
  
  def getDataTrackCodes(self, source):
  
     dataGroup = self._getDataTrackGroup(source)
     
     return sorted(dataGroup.keys())
  
  
  def getDataTrack(self, source, code, chromosomes=None, model=None):
    """Get a list of superposed data values and annotations"""
  
    dataGroup = self._getDataTrackGroup(source)

    if code in dataGroup:
      dataLayer = dataGroup[code]
      refLayer  = self.getRefDataTrackGroup(source, code)
    
      stranded = dataLayer.attrs['stranded'] # 0, 1
      
      if not chromosomes:
        chromosomes = self.getChromosomes()
      
      dataDict = {}
      
      for chromo in chromosomes:
        if chromo not in refLayer:
          continue
        
        regionData = array(refLayer[chromo]['regions']) # start, end
        valueData  = array(refLayer[chromo]['values']) # origValue, normValue
        
        if (model is not None) and ('models' in refLayer[chromo]):
          sIndices = (refLayer[chromo]['models'] == model).nonzero()
          regionData = regionData[sIndices]
          valueData = valueData[sIndices]
       
        else:
          sIndices = None
        
        if 'annotations' in refLayer[chromo]:
          annotations = refLayer[chromo]['annotations']
          
          if sIndices:
            annotations = annotations[sIndices]
          
        else:
          annotations = None
 
        if stranded:
          starts = array(regionData[:,0])
          ends = array(regionData[:,1])
          opposites = (ends < starts).nonzero()
 
          strands = ones(len(regionData), int)
          strands[opposites] = -1
 
          regionData[:,0][opposites] = ends
          regionData[:,1][opposites] = starts
 
        else:
          strands = None
 
        dataDict[chromo] = (regionData, valueData, strands, annotations)
        
      return dataDict
  
  
  def _setDataTrackIndicesMb(self, source, code, chromosomes=None):
    """
    Sort regions and vales of data track according to seq start position
    Store indices for first region of every integral Mb
    """
    
    refNuc = self.getRefDataTrackNuc(source, code, 'a')
    refLayer = refNuc._getDataTrackGroup(source)[code]
   
    if chromosomes is None:
      chromosomes = [c for c in refLayer]
   
    for chromo in chromosomes:
      if chromo not in refLayer:
        continue
      
      if chromo not in self.chromosomes:
        continue
      
      subGroup = refLayer[chromo]
      regions = array(subGroup['regions'], int32)   
      values = array(subGroup['values'], float)
      
      starts = regions[:,1]
      order = argsort(starts)
      
      regions = regions[order]
      values = values[order]
      
      refNuc._setData('values', subGroup, float, values)    
      refNuc._setData('regions', subGroup, int32, regions)    
      
      if 'annotations' in subGroup:
        annotations = array(subGroup['annotations'])[order]
        refNuc._setData('annotations', subGroup, VL_str, annotations)
 
      if 'models' in subGroup:
        modelArray = array(subGroup['models'], uint32)[order]
        refNuc._setData('models', subGroup, uint32, modelArray)
         
      start, end = self.getChromosomeLimits(chromo)
      n = int(ceil(end/1e6))
      
      starts /= 1000000
      
      indices = zeros(n, int32)
      prev = -1
      for i, j in enumerate(starts):
        if j != prev:
          indices[prev+1:j+1] = i
          prev = j

      if j < len(indices):
        indices[j:] = indices[j]
        
      subGroup.attrs['indicesMb'] = indices
    
    if refNuc is not self:
      refNuc.save()
  
  
  def _setDataTrackHistogramCache(self, refDataGroup, source, code, chromo, threshold, maxNorm):
    
    from cUtil.dataLayer import regionBinValues
     
    start, end = self.getChromosomeLimits(chromo)
    regions = array(refDataGroup[chromo]['regions'], int32)    # start, end
    values = array(refDataGroup[chromo]['values'], float)[:,1] # origValue, normValue
    
    binValues = regionBinValues(regions, values, int32(DATA_TRACK_CACHE_BIN),
                                int32(start), int32(end), maxNorm, 1.0, threshold)
    
    indices = binValues.nonzero()
    
    starts = arange(0, DATA_TRACK_CACHE_BIN*len(binValues), DATA_TRACK_CACHE_BIN, dtype=int32)
    starts += start
    binRegions = vstack([starts, starts+(DATA_TRACK_CACHE_BIN-1)]).T

    binValues = binValues[indices]
    binRegions = binRegions[indices]
    
    key = (source, code, chromo)
    self._dataTrackHistogramCache[key] = (threshold, binValues, binRegions)   
    
    return binValues, binRegions
    

  def getDataTrackPixmap(self, source, code, chromo, start, end, width=800, height=200):
    """Get an RGB pixmap array (height, width RGB) representing a data track
       for a chromosome region, for use in genome browsers etc."""
    
    from cUtil.dataLayer import regionBinValues, addPixmapHistogram
    
    dataGroup = self._getDataTrackGroup(source)
    refLayer = self.getRefDataTrackGroup(source, code)
    pixmap = zeros((height, width, 4), dtype=uint8)
    
    if (code in dataGroup) and (chromo in refLayer):
      dataLayer = dataGroup[code]
      stranded = dataLayer.attrs['stranded']          # 0, 1
      optAttrs = list(dataLayer.attrs['options'])
      asDensity = True if optAttrs[3] else False
      dispAttrs = list(dataLayer.attrs['display'])
      color = [int(255*rgb) for rgb in dispAttrs[:4]] # 8-bit RGB
      scale = dispAttrs[4]                            
      thrsh = dispAttrs[5]
      numBins = float(width)
      binSize = int32((end-start)/numBins)
      maxOrig, maxNorm = self.getDataTrackMaxValues(source, code)
      
      
      
      if binSize > 3 * DATA_TRACK_CACHE_BIN:
        key = (source, code, chromo)
 
        if key in self._dataTrackHistogramCache:
          cacheThrsh, cacheValues, cacheRegions = self._dataTrackHistogramCache[key]
 
          if cacheThrsh == thrsh:
            values = cacheValues
            regions = cacheRegions
 
          else:
            values, regions = self._setDataTrackHistogramCache(refLayer, source, code,
                                                               chromo, thrsh, maxNorm)
 
        else:
          values, regions = self._setDataTrackHistogramCache(refLayer, source, code,
                                                             chromo, thrsh, maxNorm)
      
      else:      
        
        if not len(refLayer[chromo].attrs['indicesMb']):
          self._setDataTrackIndicesMb(source, code)
          
        indicesMb = refLayer[chromo].attrs['indicesMb']
        i = indicesMb[int(start/1e6)]
        j = indicesMb[min(int(end/1e6)+1, len(indicesMb)-1)]
        
        if i == j:
          return pixmap
        
        regions = array(refLayer[chromo]['regions'][i:j], int32)  # start, end
        values = array(refLayer[chromo]['values'][i:j,1], float)  #  normValue  
        
      hist = regionBinValues(regions, values, binSize, int32(start), int32(end), maxNorm, scale, thrsh)
                
      colors = array(color + [0, 0, 0, 255] + color, int32) # pos, zero, neg
      pixmap = addPixmapHistogram(pixmap, hist[:width], colors, asDensity)
    
    return pixmap
    

  def getChromosomeDataValues(self, code, chromosome, source=EXTERNAL, model=None):
    """Get a list of chromsome data values in stored order, without position information"""
  
    dataGroup = self._getDataTrackGroup(source)

    if code in dataGroup:
      dataLayer = self.getRefDataTrackGroup(source, code)

      if chromosome not in dataLayer:
        return 
     
      valueData  = array(dataLayer[chromosome]['values']) # origValue, normValue
      
      if (model is not None) and ('models' in dataLayer[chromosome]):
        sIndices = (dataLayer[chromosome]['models'] == model).nonzero()
        valueData = valueData[sIndices]
       
      return valueData


  def getDataTrackRegions(self, code, chromosome, source=EXTERNAL, model=None):
    """For a given chromosomes get a list of superposed regions with data values"""
   
    dataGroup = self._getDataTrackGroup(source)

    if code not in dataGroup:
      return
      
    dataLayer = self.getRefDataTrackGroup(source, code)
    if chromosome not in dataLayer:
      return  
    
    chromoData = dataLayer[chromosome]
    
    regionData = array(chromoData['regions']) # start, end
      
    if (model is not None) and ('models' in dataLayer[chromosome]):
      sIndices = (dataLayer[chromosome]['models'] == model).nonzero()
      regionData = regionData[sIndices]
    
    return regionData
      
      
  def getImageGrid(self, code):
    """Get a full volumetric grid"""

    if code in self.grid:
      gridData = self.grid[code]
      origin = gridData.attrs['origin']
      gridSize = gridData.attrs['gridSize']
      gridData = array(gridData) # Values[i,j,k]
      
      return origin, gridSize, gridData


  def getImageCoords(self, code):
    """Get coordinates of volumetric signal intensities (generally for sparse data)"""

    imageGroup = self.images

    if code in imageGroup['coords']:
      return array(imageGroup['coords'][code]) # [(x, y, z, (value))]

  
  def getFoci(self, code):
    """Get lists of spatial signal peaks/foci"""
    
    fociGroup = self.foci
    
    if code in fociGroup:
      data = array(fociGroup[code]['peaks']) # x, y, z, height, volume, dx, dy, dz
      
      if 'annotations' in fociGroup[code]:
        annotations = list(array(fociGroup[code]['annotations']))
      else:
        annotations = None
        
      return data, annotations
  
  
  def importCoordFile(self, fileName, chromosomes=None):
    """Import 3D chromosome structure data from a file with 'chromo pos X Y Z' data format.
       Comment lines starting with # are ignored.  Columns can be separated by any white space character.
       If chromosomes is given, only those chromosomes will be read in."""
    
    alc = ' '
    ins = ' '

    #open file to be read
    fileObj = open(fileName, 'r')

    positions = {}
    coordinates = {}

    #read lines
    for line in fileObj:

      if line == "\n" or line[0] == "#":
        continue

      a = line.rstrip().split()
      #check number of entries
      if len(a) != 5:
        sys.stderr.write("ERROR importCoordFile: found %d (!= 5) entries in the file\n%s" % (len(a), line))
        #TODO: what to do?
        #raiseException()
        sys.exit(1)

      #read entries
      chromo, pos, x, y, z = a

      #if we only want a given list of chromosomes
      if chromosomes:
        if chromo not in chromosomes:
          continue

      #collect positions and coordinates
      if chromo not in positions.keys():
        positions[chromo] = []
        coordinates[chromo] = []
      positions[chromo].append(int(pos)) #TODO: longint?
      coordinates[chromo].append(map(float,[x,y,z]))


    #save entries
    #chromosomes and positions
    current = self.getChromosomes()
    self.removeChromosomes(current)
    self.addChromosomes(positions)
    #coordinates
    for chromo in positions.keys():
      self.setAllCoords(array(coordinates[chromo]),[chromo])

    fileObj.close()

    return

  
  def exportPdbFile(self, fileName, chromosomes=None, scale=0.15, extended = False):
    """Save any 3D chromosome structure data in a pseudo PDB format"""
    
    alc = ' '
    ins = ' '
    prefix = 'HETATM'
    tlc = 'CHR'
    lFormat = '%-80.80s\n'
    if extended:
      pdbFormat = '%-6.6s%5.1d %4.4s%s%3.3s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  %10d\n'
      terFormat = '%-6.6s%5.1d      %s %s%4.1d%s                                                     %10d\n'
    else:
      pdbFormat = '%-6.6s%5.1d %4.4s%s%3.3s %s%4.1d%s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2.2s  \n'
      terFormat = '%-6.6s%5.1d      %s %s%4.1d%s                                                     \n'

    fileObj = open(fileName, 'w')
    
    sample = self.sample
    infoDict = self.getSampleInfo()
    
    line = 'TITLE     %sd' % (sample.attrs['name'])
    fileObj.write(lFormat % line)
    fileObj.write(lFormat % 'REMARK 210')
    
    keys = infoDict.keys()
    keys.sort()
    
    for key in keys:
      fileObj.write(lFormat % 'REMARK 210 %s:%s' % (key, infoDict[key])) 
 
    fileObj.write(lFormat % 'REMARK 210 Atom type C are backbone nodes')
    fileObj.write(lFormat % 'REMARK 210 Atom type N are restrained nodes')
    
    seqSep = self.getBackboneSpacing() // 1e6
    
    fileObj.write(lFormat % 'REMARK 210 Residue number increases every %.2f Mb node' % seqSep)
    fileObj.write(lFormat % 'REMARK 210 Chain letter is chromosome')
    fileObj.write(lFormat % 'REMARK 210 B-factor field is sequence mapability')
    fileObj.write(lFormat % 'REMARK 210 Coordinate scale is 1 unit to %.2f microns' % scale)
    if extended:
      fileObj.write(lFormat % 'REMARK 210 Extended PDB format with the genomic position of the bead in the extra column')
    fileObj.write(lFormat % 'REMARK 210')
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    numModels = self.getNumModels()
    chromoGroup = self.chromosomes
    
    for m in range(numModels):
      line = 'MODEL     %4d' % (m+1)
      fileObj.write(lFormat  % line)
      
      c = 0
      j = 1
      seqPrev = None
      for chromo in chromosomes:
        coords = self.getModelCoords(m, [chromo,]) * scale
        
        if not len(coords):
          continue
        
        backbone = array(chromoGroup[chromo]['backbone'])
        positions = array(chromoGroup[chromo]['positions'])
        mapability = array(chromoGroup[chromo]['mapability'])
        
        for i, pos in enumerate(positions):
          c += 1
 
          seq = int(positions[i]//1e6) + 1
          
          if seq == seqPrev:
            j += 1
          else:
            j = 1
          
          if backbone[i]:
            el = 'C'
            a = 'C%d' % j
          else:
            el = 'N'
            a = 'N%d' % j
            
          aName = '%-3s' % a
          b = mapability[i]
          x, y, z = coords[i] #XYZ coordinates
          if extended:
            g = positions[i] #genomic position
          seqPrev = seq
 
          if extended:
            line  = pdbFormat % (prefix,c,aName,alc,tlc,chromo,seq,ins,x,y,z,0.0,b,el,g)
          else:
            line  = pdbFormat % (prefix,c,aName,alc,tlc,chromo,seq,ins,x,y,z,0.0,b,el)
          fileObj.write(line)
 
        fileObj.write(lFormat  % 'ENDMDL')
 
    for i in range(c-2):
       line = 'CONECT%5.1d%5.1d' % (i+1, i+2)
       fileObj.write(lFormat  % line)
 
    fileObj.write(lFormat  % 'END')
    fileObj.close()
  
  
  def exportContactFile(self, fileName, groupName, format='tsv', chromosomes=None,
                        cis=True, trans=True, model=None):
    """Export chromosomal contacts as formatted text file"""
    
    fileObj = open(fileName, 'w')
    write = fileObj.write

    contactDict = self.getContacts(groupName, chromosomes, cis, trans, model)
    keys = contactDict.keys()
    keys.sort()
    
    if format in ('tsv', 'csv'):
      heads = ('chromoA', 'posA', 'chomoB', 'posB', 'numObs')
    
      if format == 'csv':
        separator = ','
        template = '%s,%d,%s,%d,%d\n'
      else:
        separator = '\t'
        template = '%s\t%d\t%s\t%d\t%d\n'
      
      line = separator.join(heads) + '\n'
      write(line)
      
      for key in keys:
        chrA, chrB = key
 
        contacts = contactDict[key].T
 
        for posA, posB, numObs in contacts:
          line = template % (chrA, posA, chrB, posB, numObs)
          write(line)
 
    
    fileObj.close()

    # More formats....
  
  
  def _checkImageFileName(self, fileName, format):
  
    root, fex = os.path.splitext(fileName)
    
    if format == 'PNG':
     if fex.upper() != 'PNG':
       fileName = root + '.png'
    
    elif format == 'JPG':
     if fex.upper() not in ('JPG','JPEG','JPE', 'JIF', 'JFIF', 'JFI'):
       fileName = root + '.jpg'
       
    return fileName
  
  
  def _getPixmapImage(self, pixmap):
    
    # TBD add more olour options
    from util.Image import pixmapToImage
    
    maxVal = pixmap.max()
    if maxVal:
      pixmap /= maxVal
    
    if pixmap.ndim == 2:
      red = array(pixmap) ** 1.5
      green = pixmap*255.0
      blue = pixmap*255.0
      pixmap = dstack([red, green, blue])
      
    image = pixmapToImage(pixmap)

    return image
  
  def getContactImage(self, chromosomes=None, binSize=int(1e6), gamma=0.5, model=None):
  
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    matrices = []
    for chrA in chromosomes:
      row = []
      
      for chrB in chromosomes:
        matrix = array(self.getContactMatrix(chrA, chrB, binSize, model), float)
        maxVal = matrix.max()
        if maxVal:
          matrix /= maxVal
        
        if gamma:
          matrix = power(matrix, gamma)
        
        row.append(matrix)
      matrices.append(row)
      
    n = sum([len(m[0]) for m in matrices[0]])
    nPix = len(matrices) + n + 1
    pixmap = ones((nPix, nPix), float) * 0.1
    
    y = 1
    for i, row in enumerate(matrices):
      x = 1
      for j, matrix in enumerate(row):
        dy, dx = matrix.shape
        pixmap[y:y+dy, x:x+dx] = matrix
        matrices[i][j] = None
        x += dx + 1
      y += dy + 1
    
    return self._getPixmapImage(pixmap[::-1])
    
    
  def exportContactImage(self, fileName, format='PNG', chromosomes=None,
                         binSize=int(1e6), gamma=0.5, model=None):
    """Export chromosomal contact maps as an image file"""
    
    from util.Image import pixmapToImage
    
    fileName = self._checkImageFileName(fileName, format)
    
    image = self.getContactImage(chromosomes, binSize, gamma, model)
    image.save(fileName, format)
    
    return image
    
    
  def exportRmsdMatrixImage(self, fileName, format='PNG', models=None, backboneOnly=False):
    """Export a pairwise model RMSD matrix as an image file"""
    
    if not models:
      models = range(self.getNumModels())
    
    fileName = self._checkImageFileName(fileName, format)
  
    matrix = self.calcModelRmsdMatrix(models, backboneOnly)
    image = self._getPixmapImage(matrix)
    image.save(fileName, format)
    
    return image
    

  def exportDistanceMatrixImage(self, fileName, chromosome, format='PNG',
                                models=None, binSize=int(1e6), useMinVal=True):
    """Export a chromosome distance matrix as an image file"""
 
    if not models:
      models = range(self.getNumModels())
    
    fileName = self._checkImageFileName(fileName, format)
    from itil.Image import pixmapToImage
  
    matrix = self.getDistanceMatrix(chromosome, models, binSize, useMinVal)
    maxVal = matrix.max()
    if maxVal:
      matrix /= maxVal
    
    image = self._getPixmapImage(matrix)
    image.save(fileName, format)

    return image
    
    
  def exportDataTrack(self, source, code, fileName, format='tsv', chromosomes=None):
    """Export a layer of chromosomal data that may be superposed on structures to an external file"""
    
    # TBD more formats WIG, BED, DAT, BIOMART?, SAM, BAM, .peaks
    
    dataDict = self.getDataTrack(source, code, chromosomes)
    
    if dataDict:
      format = format.lower()
      fileObj = open(fileName, 'w')
      write = fileObj.write
    
      chromos = dataDict.keys()
      chromos.sort()
      
      if format in ('tsv', 'csv'):
        heads = ('chr','start','end','origValue','value','label')
 
        if format == 'csv':
          separator = ','
          template = '%s,%d,%d,%.6f,%.6f,%s\n'
        else:
          separator = '\t'
          template = '%s\t%d\t%d\t%.6f\t%.6f\t%s\n'
 
        line = separator.join(heads) + '\n'
        write(line)
      
        for chromo in chromos:
          regionData, valueData, strands, annotations = dataDict[chromo]
        
          for i, region in enumerate(regionData):
            start, end = region
            
            if strands is not None:
              if strands[i] < 0:
                start, end = end, start
               
            if annotations is not None:
              annotation = annotations[i]
            else:
              annotation = ''
            
            origValue, value = valueData[i]
            
            line = template % (chromo, start, end, origValue, value, annotation)
            write(line)
      
      elif format == 'bed':
        template = 'chr%s\t%d\t%d\t%s\t%d\t%s\n' # chr, start, end, label, score, strand
        
        for chromo in chromos:
           regionData, valueData, strands, annotations = dataDict[chromo]
 
           for i, region in enumerate(regionData):
             start, end = sorted(region)
             
             if strands is not None:
               if strands[i] > 0:
                 strand = '+'
               else:
                 strand = '-'
             
             else:
               strand = '+'
 
             if annotations is not None:
               label = annotations[i]
             
             else:
               label = '%d' % i
 
             origValue, value = valueData[i]
             
             score = int(value * 1000)
  
             line = template % (chromo, start, end, label, score, strand)
             write(line)

      elif format == 'bedgraph':
        write('track type=bedGraph name="%s"\n' % code)
        template = 'chr%s\t%d\t%d\t%.5f\n' # chr, start, end, value
        
        for chromo in chromos:
           regionData, valueData, strands, annotations = dataDict[chromo]
 
           for i, region in enumerate(regionData):
             start, end = sorted(region)
             origValue, value = valueData[i] 
             line = template % (chromo, start, end, origValue)
             write(line)

      elif format == 'broadpeak':
        template = 'chr%s\t%d\t%d\t%s\t%d\t%s\t%.5f\t-1\t-1\n'
        
        for chromo in chromos:
           regionData, valueData, strands, annotations = dataDict[chromo]
 
           for i, region in enumerate(regionData):
             start, end = sorted(region)
             
             if strands is not None:
               if strands[i] > 0:
                 strand = '+'
               else:
                 strand = '-'
             
             else:
               strand = '+'
 
             if annotations is not None:
               label = annotations[i]
             
             else:
               label = '%d' % i
 
             origValue, value = valueData[i]
             
             score = int(value * 1000)
  
             line = template % (chromo, start, end, label, score, strand, origValue) # , pValue, qValue)
             write(line)
     
      elif format == 'wig':
        template = 'variableStep chrom=chr%s span=%d\n%d %.3f\n'
        
        for chromo in chromos:
           regionData, valueData, strands, annotations = dataDict[chromo]
           
           for i, region in enumerate(regionData):
            start, end = sorted(region)
            origValue, value = valueData[i]
            delta = end-start
            line = template % (chromo, delta, start, origValue)
 
            write(line)
      
      fileObj.close()


  def importSamFile(self, filePath, hdfGroup, binSize=None, cis=True, trans=True):

    # TBD: Needs more checks to make sure it is a properly processed and paired BAM file
     
    from util.Io import isFileBinary
    from cUtil.samread import readPairedSam
    from cUtil.apiUtil import binContacts, intMatrixToSparse
    import time 
    
    t0 = time.time()
    
    if isFileBinary(filePath):
      fileMode = 'rb'
    else:
      fileMode = 'r'
    
    dataDict, positionDict = readPairedSam(filePath, fileMode)
    
    if binSize and binSize > 1:
      binSize = int32(binSize)
      
      for chrA in dataDict:
        sa, ea = positionDict[chrA]
        n = ceil((ea-sa)/float(binSize))
        if n <= 0:
          continue
 
        subGroup = self._getGroup(chrA, hdfGroup)

        for chrB in dataDict[chrA]:
          sb, eb = positionDict[chrB]
          m = ceil((eb-sb)/float(binSize))
          if m <=0:
            continue
          
          contacts = dataDict[chrA][chrB].T
          binMatrix = zeros((n,m), int32)
          binContacts(contacts, binMatrix, int32(sa), int32(sb), binSize)
          
          contacts = intMatrixToSparse(binMatrix, binSize, binSize, int32(sa), int32(sb))
          self._setData(chrB, subGroup, uint32, contacts, compression='gzip')
   
    else:
      for chrA in dataDict:
        subGroup = self._getGroup(chrA, hdfGroup)
 
        for chrB in dataDict[chrA]:
          contacts = dataDict[chrA][chrB].T
          self._setData(chrB, subGroup, uint32, contacts, compression='gzip')
    
    print('BAM read time: %.4f' % (time.time()-t0))
 
    return positionDict
    
  
  def importContactArray(self, dataArray, chromoA, chromoB, posA, posB,
                         groupName='singleCell', isSingleCell=True, updateChromos=True):
    
    sizeA, sizeB = dataArray.shape
    
    if len(posA) != sizeA:
      msg = 'Positions for chromosome A must be of length %d' % sizeA
      raise Exception(msg)
    
    if len(posB) != sizeB:
      msg = 'Positions for chromosome B must be of length %d' % sizeB
      raise Exception(msg)
    
    group = self._getContactGroup(groupName, self.origContacts)
    group.attrs['isSingleCell'] = 1 if isSingleCell else 0
    
    if isSingleCell:
      group.attrs['binSize'] = 1
    else:
      posA = array(posA)
      deltas = posA[1:] - posA[:-1]
      group.attrs['binSize'] = deltas.min()
    
    contacts = []
    contactsAppend = contacts.append
    
    for i in range(sizeA):
      for j in range(sizeB):
        numObs = dataArray[i,j]
        
        if numObs:
          contactsAppend( (posA[i], posB[j], numObs)  )    
  
    contacts = array(contacts, uint32).T
    subGroup = self._getGroup(chromoA, group)
    self._setData(chromoB, subGroup, uint32, contacts, compression='gzip')
    
    if groupName in self._contactsCache:
      del self._contactsCache[groupName]
      
    if updateChromos:    
      addDict = {}
      chromos = set(self.getChromosomes())
      
      if chromoA not in chromos:
        addDict[chromoA] = [min(posA), max(posA)]

      if chromoB not in chromos:
        addDict[chromoB] = [min(posB), max(posB)]
     
      if addDict:
        self.addChromosomes(addDict) 
        
           
  def importContactFile(self, filePath, format='text', groupName='singleCell',
                        binSize=None, chromosomes=None, cis=True, trans=True,
                        isSingleCell=True, updateChromos=True):
    """Import chromosomal contacts from a formatted text file or BAM/SAM file.
       Effectively defines the .nuc file, if needed."""
       
    # Removing existing data?
    
    if not self._filePathExists(filePath):
      return
    
    from util.Io import readListFile
    
    chromoGroups = self.chromosomes
    
    if binSize is not None:
      binSize = int(binSize)
    
    group = self._getContactGroup(groupName, self.origContacts)
    group.attrs['isSingleCell'] = 1 if isSingleCell else 0
    group.attrs['filePath'] = string_(filePath)

    if isSingleCell:
      group.attrs['binSize'] = 1
    else:
      group.attrs['binSize'] = binSize or 1
    
    if format == 'sam':
      posDict = self.importSamFile(filePath, group, binSize, cis, trans)
    
    else: # format == 'text':
      converters = (None, int, None, int)
      data = readListFile(filePath, converters, None, skipFirst=True)
      
      contactDict = {}
      posDict = {}
      for chrA, posA, chrB, posB in data:
 
        if (chrA == chrB) and not cis:
          continue
 
        if (chrA != chrB) and not trans:
          continue
 
        if binSize:
          posA = binSize * int(posA//binSize)
          posB = binSize * int(posB//binSize)
 
        key = tuple(sorted([(chrA, posA), (chrB, posB)]))
 
        if key in contactDict:
          contactDict[key] += 1
        else:
          contactDict[key] = 1
        
      dataDict = {}
      for key in contactDict:
        pairA, pairB = key
        chrA, posA = pairA
        chrB, posB = pairB
        numObs = contactDict[key]
        
        if chrA in posDict:
          s, e = posDict[chrA]
          posDict[chrA] = (min(s, posA), max(e, posA))
        else:
          posDict[chrA] = (posA, posA)       
        
        if chrB in posDict:
          s, e = posDict[chrB]
          posDict[chrB] = (min(s, posB), max(e, posB))
        else:
          posDict[chrB] = (posB, posB)       
 
        chromoPair = (chrA, chrB)
 
        if chromoPair in dataDict:
          dataDict[chromoPair].append( (posA, posB, numObs) )
        else:
          dataDict[chromoPair] = [(posA, posB, numObs)]
 
      chromoPairs = sorted(list(dataDict.keys()))
      for chromoPair in chromoPairs:
        contacts = array(dataDict[chromoPair]).T
        chromoA, chromoB = chromoPair
 
        subGroup = self._getGroup(chromoA, group)
        self._setData(chromoB, subGroup, uint32, contacts, compression='gzip')
    
    if updateChromos:    
      addDict = {}
      chromos = set(self.getChromosomes())
      
      for chrA, limits in posDict.iteritems():
        if chrA not in chromos:
          addDict[chrA] = limits
     
      if addDict:
        self.addChromosomes(addDict)
 
    if groupName in self._contactsCache:
      del self._contactsCache[groupName]
    
    
  def importDataTrack(self, filePath, source, code, format=None):
    """Import a layer of chromosomal data that may be superposed on structures from an external file"""
    
    # Will need to ask about merge or overwrite (del first) of existing codes
    
    if not self._filePathExists(filePath):
      return
     
    from util.Io import readListFile, guessTextFileFormat
    
    dataGroup = self.externalData
    
    if format == 'tsv':
      converters, skipFirst, haveEnds = guessTextFileFormat(filePath, '\t')
      dataList = readListFile(filePath, converters, '\t', skipFirst)
     
    elif format == 'csv':
      converters, skipFirst, haveEnds = guessTextFileFormat(filePath, ',')
      dataList = readListFile(filePath, converters, ',', skipFirst)
       
    else: # Whitespace separated
      converters, skipFirst, haveEnds = guessTextFileFormat(filePath, None)
      dataList = readListFile(filePath, converters, None, skipFirst)
    
    self.setDataTrackList(source, code, dataList, haveEnds, None)
      
    # formats WIG, BED, DAT, BIOMART?, SAM, BAM, .peaks
    # Must have existing chromosome representation?

  
  def calcStructure(self, contacts, numCpus=1, updateFunc=None, trans=True, bgCalc=False):
    """Calculate structure using stored parameters"""
    
    chromosomes = self.getDisplayedChromosomes()
    attrs = self.structures['calculation'].attrs
  
    numModels    = attrs['numModels']
    bboneReg     = attrs['bboneReg']
    bboneSpace   = attrs['bboneSpace'] * 1000    # Stored in kb
    isBinned     = bool(attrs['restrBinned'])
    powerLaw     = attrs['powerLaw']
    seqScale     = attrs['seqUnitScale'] * 1000  # Stored in kb
    restrDist    = attrs['restrDist']
    restrErr     = attrs['restrErr']
    tempMax      = attrs['tempMax']
    tempMin      = attrs['tempMin']
    tempSteps    = attrs['tempSteps']
    dynSteps     = attrs['dynSteps']
    hierProtocol = attrs['hierProtocol']
    hierStart    = attrs['hierStart'] * 1000000  # Stored in Mb
    hierSteps    = attrs['hierSteps']
    startStruc   = attrs['startStruc']
    randRad      = attrs['randRad']
    randSeed     = attrs['randSeed'] or None
          
    calcParams = StrucCalcParams(distPowerLaw=powerLaw,
                                 distLower=(1.0-restrErr) * restrDist,
                                 distUpper=(1.0+restrErr) * restrDist,
                                 bboneLower=0.1,
                                 bboneUpper=1.1,
                                 seqScale=seqScale,         # Does not apply to binned restraints
                                 randSeed=randSeed,
                                 randStart=startStruc != 2, # 2 => Keep prev coords
                                 randWalk=startStruc == 1,  # 1 => Use random walk
                                 randRad=randRad)
    
    # Any initial hierarchical stages
    if hierProtocol:
      bestRes = bboneSpace if isBinned else seqScale
      stages = [({}, (hierStart, hierStart))]
      
      for i in range(1,hierSteps):
        frac = (hierSteps-i)/float(hierSteps)
        bboneSep = int( exp(log(hierStart/bestRes) * frac) * bestRes)
        stages.append(({}, (bboneSep, bboneSep) ))
    
    else:
      stages = []
    
    # Final, ultimate resolution stage
    if isBinned: 
      stages.append( ({}, (bboneSpace, bboneSpace)) ) # Binned restraints, regular backbone
      
    elif bboneReg:
      stages.append( (None, (bboneSpace, bboneSpace)) ) # No binning, regular spacer backbone
      
    else:
      stages.append( (None, None) ) # No binning, no spacers
    
    # Setup an annealing stage objects
    tempStart = tempMax
    deltaT = (tempMax-tempMin)/float(len(stages))
    
    for domainDict, sepBounds in stages:
      calcParams.addAnnealStage(domainDict, sepBounds, tempStart, tempStart-deltaT,
                                tempSteps, dynSteps, timeStep=0.001, useTrans=trans)
      tempStart -= deltaT
                        
    # Run the calculation
    job = self.annealStructure(contacts, chromosomes, numModels, calcParams,
                               updateFunc, bgCalc, numCpus)
                        
    return job
    
  
  def annealStructure(self, groupName, chromosomes=None, numModels=1, calcParams=None,
                      callback=None, bgCalc=False, numCpus=1, contactModel=None):
                      
    """Calculate chromosome structures using distance restraints from
       contact data using simple annealing protocol"""
    
    from cUtil.apiUtil import calcRestraints, concatenateRestraints
    
    if not calcParams:
      calcParams = StrucCalcParams() # i.e. use defaults
      calcParams.addAnnealStage()    # One default, binned stage  

    # Future: consider masses and fixed particles    
     
    # Some initialisation 
     
    if calcParams.randSeed is None:
      seed()
    else:
      seed(calcParams.randSeed)
    
    if chromosomes:
      chromosomes = sorted(chromosomes)
    else:
      chromosomes = self.getChromosomes()
      
    models = list(range(numModels))
    contactDict  = self.getCachedContacts(groupName)    
    isSingleCell = self.origContacts[groupName].attrs['isSingleCell']
    contactModel = -1 if contactModel is None else contactModel
    
    
    # Loop through annealing stages, collate restraints and calc parameters
    
    restrIndices = []  # Particle indices of restrained pairs
    restrDists = []    # Distrance restraint limits
    restrAmbig = []    # Restraint ambiguity data
    posDicts = []      # Particle seq positions for each chromosome
    temps = []         # Temperature schedule
    repulsScales = []  # Repulsive scale schedule
    timeSteps = []     # Time interval between updating particle trajectories
    dynSteps = []      # Particle dynamics steps for each temperature
    stageCounts = []   # Num restraints in each stage
    masses = []        # Particle mass, propotional to seq content
    radii = []         # Spherical particle radius, based on mass (seq content)
    
    self._cacheRestrDict = [] # For display
    
    for i, annealStage in enumerate(calcParams.annealStages):   
      
      if annealStage.bboneSep is not None: # TBD: Think about what to do here in the long-run
        seqScale = annealStage.bboneSep[1]
      else:
        seqScale = calcParams.seqScale       
           
      # calc/bin restraints for chrosen chromos from contacts  annealStage.bboneSep
      restrDict, posDict, bboneDict, binLstDict = calcRestraints(chromosomes, contactDict, isSingleCell,
                                                     annealStage.bboneSep, annealStage.domainDict, 1.0, # scale not used
                                                     calcParams.distPowerLaw, calcParams.distLower, calcParams.distUpper,
                                                     calcParams.minNumObs, calcParams.maxPopDist, contactModel)
      
      if not bgCalc:
        self._cacheRestrDict.append(restrDict)
      
      if i == 0:
        # For first stage setup starting coords for seq positions
        
        self.addChromosomes(posDict, bboneDict, interpolateCoords=True)

        coords = []
        if calcParams.randStart:
          self.setRandomCoords(models, chromosomes, maxStep=calcParams.randRad/100.0,
                               randWalk=calcParams.randWalk, randSeed=calcParams.randSeed)
 
        for model in models:
          mCoords = self.getModelCoords(model, chromosomes)
 
          if (mCoords is None) or not len(mCoords):
            self.setRandomCoords([model], chromosomes,
                                 randWalk=calcParams.randWalk,
                                 randSeed=calcParams.randSeed)
            mCoords = self.getModelCoords(model, chromosomes)
 
          coords.append(mCoords)
 
        coords = array(coords)
      
      # add backbone restraints, get single concatenated restraint arrays      
      indices, dists, ambig = concatenateRestraints(restrDict, posDict, seqScale,
                                                    calcParams.bboneLower, calcParams.bboneUpper)
      
      stageMasses = self.massesFromBins(binLstDict)

      masses.append( stageMasses ) 
      radii.append( self.radiiFromMasses(stageMasses) )
      
      restrIndices.append(indices)
      restrDists.append(dists)
      restrAmbig.append(ambig)
      stageCounts.append(len(ambig))
      posDicts.append(posDict)
      
      # construct annealing schedule for this stage
      temps.append(annealStage.getTempSchedule())
      repulsScales.append(annealStage.getRepulsionSchedule())
      timeSteps.append(annealStage.timeStep)
      dynSteps.append(annealStage.dynSteps)

    masses = concatenate(masses, axis=0)
    radii = concatenate(radii, axis=0)
    restrIndices = concatenate(restrIndices, axis=0)
    restrDists = concatenate(restrDists, axis=0)
    restrAmbig = concatenate(restrAmbig, axis=0)

    
    # run the annealing via wrapper that handles temp schedule, repulsion and parallelisation
    job = self._runAnnealing(coords, posDicts, stageCounts, restrIndices, restrDists, restrAmbig,
                             temps, repulsScales, timeSteps, dynSteps, callback, bgCalc, numCpus,
                             masses=masses, radii=radii)
      
    # store final restraints
    restraintGroup = self.restraints
    for chrA in restrDict:
      rGroup = self._getGroup(chrA, restraintGroup)
    
      for chrB in restrDict[chrA]:
        dataArray = restrDict[chrA][chrB]
        self._setData(chrB, rGroup, float32, dataArray)
    
    # store seq positions for coords 
    self.addChromosomes(posDict, bboneDict, interpolateCoords=False)
    
    if not bgCalc:
      # store final coords
      coords = job.getResult()
      self.setAllCoords(coords, chromosomes)
 
      # superimpose ensemble
      if len(coords) > 1:
        self.modelAlign(chromosomes=chromosomes)
 
      if callback:
        callback()
 
      self.save()
        
    return job

  def massesFromBins(self, binLstDict):

    masses = []
    for key in sorted(binLstDict):
      masses += ([x[1] - x[0] for x in binLstDict[key]])

    masses = array(masses, float)
    masses /= sum(masses)
    masses *= masses.shape[0]

    return array(masses)

  def radiiFromMasses(self, massesArray):

    radii = []
    for mass in massesArray:
      radii.append(mass**(1/3))
      # radii.append(2.25)

    radii = array(radii, float)
    #print(radii)

    return radii
        
      
  def _strucUpdateCallback(self, callback, newCoords, posDict, chromosomes, stage, step):

    # update restraints
    if step == 0 and (stage < len(self._cacheRestrDict)):
      restrDict = self._cacheRestrDict[stage]
      restraintGroup = self.restraints
      for chrA in restrDict:
        rGroup = self._getGroup(chrA, restraintGroup)
 
        for chrB in restrDict[chrA]:
          dataArray = restrDict[chrA][chrB]
          self._setData(chrB, rGroup, float32, dataArray)
    
    # update chromsome positions
    if step == 0:
      self.addChromosomes(posDict, interpolateCoords=False)
    
    # update coordinates
    self.setAllCoords(newCoords, chromosomes)
    callback()
    
  
  def _runAnnealing(self, coords, posDicts, stageCounts, indices, dists, ambig, temps,
                    repulsions, timeSteps, dynSteps, callback=None, bgCalc=False,
                    numCpus=1, minDist=1.50, repDist=1.5, printInterval=100,
                    masses=None, radii=None):
    """
    Wrapper to calulate structures via simulated annealing using concatenated restraint arrays.
    Controls the temperature and repulsive schedule. Handles parallelisation.
    Handles the display callback.
    """
    
    from solve.SimAnnealJob import simAnnealjob
    from parallel.Engine import Engine
      
    nModels = len(coords)
       
    engine = Engine(numCpus)

    # print(minDist)
    sameArgs = {'posDicts':posDicts, 'stageCounts':stageCounts,
                'indices':indices, 'dists':dists, 'ambig':ambig,
                'temps':temps, 'repulsions':repulsions, 'timeSteps':timeSteps,
                'dynSteps':dynSteps, 'minDist':minDist,
                'repDist':repDist, 'masses':masses, 'radii':radii}
    
    printIntervals = zeros(nModels, int)
    #printIntervals[0] = printInterval  # First model prints progress     
    
    # An input arg to specify func for graphical update
    _callback = lambda v, w, x, y, z: self._strucUpdateCallback(callback, v, w, x, y, z)
    callbackArg = ('callback', _callback) if callback and not bgCalc else None

    
    diffArgs = {'coords':coords,
                'printInterval':printIntervals}
    
    job = engine.run(simAnnealjob, sameArgs, diffArgs,
                     wait=not bgCalc, combineFunc=array,
                     callbackArg=callbackArg)
    
    if not bgCalc:
      engine.stop()
    
    return job 
  
  
  def calcBootstrapRmsds(self, chromosomes, numPartitions, numResamples, rmsdWeightScale=10.0):
    """
    Calculate structure model RMSDS within and between bootstrap partition bundles and overall
    """
    #TBD: Confidence interval on coorindate mean .e.g width of 95% tail limits
    
    from util.Structure import superimposeCoordArray, superimposeCoordPair
    
    nModels = self.getNumModels()
    strucGroup = self.structures
    coordsGroup = strucGroup['coords']
    chromos = [c for c in chromosomes if c in coordsGroup]
      
    coords = None
    for chromo in chromos:
      cCoords = array(coordsGroup[chromo])
      
      if len(coords):
        if coords is None:
          coords = cCoords
        else:
          coords = append(coords, cCoords, axis=0)        
    
    sampleSize = int(nModels/numResamples)
    partitionSize = int(sampleSize/numPartitions)
    
    for i in range(nModels):
      sample = int(i//sampleSize)
      part = int((i % sampleSize)//partitionSize)
    
    m = 0    
    modelIdx = []
    for samp in range(numResamples):
      sampleIdx = []
      for part in range(numPartitions):
        sampleIdx.append(list(range(m, m+partitionSize)))  
        m += partitionSize
      
      modelIdx.append(sampleIdx)
    
    # Within partitions - mean pairwise RMSD
    rmdsWithin = []
    for sampIdx in modelIdx:
      for partIdx in sampIdx:
        rmsds = []
      
        for i in range(partitionSize-1):
          a = partIdx[i]
          
          for j in range(i+1, partitionSize):
            b = partIdx[j]
            
            coords2, rot, rmsd, atomRmsds = superimposeCoordPair(coords[a], coords[b], rmsdWeightScale)
            rmsds.append(rmsd)
      
        rmdsWithin.append(array(rmsds).mean())
    
    # Between partitions of the same sample - mean pairwise RMSD
    rmsdsBetween = []
    for sampIdx in modelIdx:
      
      for i in range(numPartitions-1):
        partIdxA = sampIdx[i]
        
        for j in range(i+1, numPartitions):
          partIdxB = sampIdx[j]
          
          pairRmsds = []
          
          for a in partIdxA:
            for b in partIdxB:
              coords2, rot, rmsd, atomRmsds = superimposeCoordPair(coords[a], coords[b], rmsdWeightScale)
              pairRmsds.append(rmsd)
       
          rmsdsBetween.append(array(pairRmsds).mean())      
    
    # Compare everything
    bundle, rotations, rmsdsAll, atomRmsds = superimposeCoordArray(coords, rmsdWeightScale)  
    
    return array(rmdsWithin), array(rmsdsBetween), rmsdsAll
    
      
  def calcBootstrapCrossValid(self, chromosomes, numPartitions, numResamples, testContacts):
    
    # Look at distance distributions of those left out contacts
    # random k(10)-fold cross validation
    
    pass
  
  
  def annealStructureBootstrap(self, groupName, chromosomes=None,
                               numModels=10, numPartitions=10, numResamples=10,
                               calcParams=None, numCpus=None):
    """
    Calculate chromosome structures by simulated annealing with resampled
    bootstrapping - divide contacts into partitions and leave out each in turn
    nuModels is per calculation/partion, e.g. numModels x numPartitions x numSamples gives total models 
    """

    if not calcParams:
      calcParams = StrucCalcParams() # i.e. use defaults
      calcParams.addAnnealStage()    # One default, binned stage  
    
    if calcParams.randSeed is None:
      seed()
    else:
      seed(calcParams.randSeed)
    
    if chromosomes:
      chromosomes = sorted(chromosomes)
    else:
      chromosomes = self.getChromosomes()
    
    if not numCpus:
      import multiprocessing
      numCpus = multiprocessing.cpu_count()
      
    models = list(range(numModels))
    contactDict  = self.getCachedContacts(groupName)    
    isSingleCell = self.origContacts[groupName].attrs['isSingleCell']
    contactModel = -1 
    
    # for each resampling 
    # partition restraints into N sets
    for samp in range(numResamples):
    
      # calculate models with each set left out
      # store as blocks in ensemble
      for partition in range(numPartitions):
      
        # Loop through annealing stages, collate restraints and calc parameters
        restrIndices = []  # Particle indices of restrained pairs
        restrDists = []    # Distrance restraint limits
        restrAmbig = []    # Restraint ambiguity data
        posDicts = []      # Particle seq positions for each chromosome
        temps = []         # Temperature schedule
        repulsScales = []  # Repulsive scale schedule
        timeSteps = []     # Time interval between updating particle trajectories
        dynSteps = []      # Particle dynamics steps for each temperature
        stageCounts = []   # Num restraints in each stage
 
        for i, annealStage in enumerate(calcParams.annealStages):
 
          if annealStage.isBinned:
            seqScale = annealStage.bboneSep
          else:
            seqScale = calcParams.seqScale
 
          # calc/bin restraints for chrosen chromos from contacts  annealStage.bboneSep
          restrDict, posDict, bboneDict = calcRestraints(chromosomes, contactDict, isSingleCell,
                                                         annealStage.bboneSep, annealStage.isBinned, 1.0, # scale not used
                                                         calcParams.distPowerLaw, calcParams.distLower, calcParams.distUpper,
                                                         calcParams.minNumObs, calcParams.maxPopDist, contactModel)

          if i == 0:
            # For first stage setup starting coords for seq positions
 
            self.addChromosomes(posDict, bboneDict, interpolateCoords=False)

            coords = []
            if calcParams.randStart:
              self.setRandomCoords(models, chromosomes, maxStep=calcParams.randRad/100.0,
                                   randWalk=calcParams.randWalk, randSeed=calcParams.randSeed)
 
            for model in models:
              mCoords = self.getModelCoords(model, chromosomes)
 
              if (mCoords is None) or not len(mCoords):
                self.setRandomCoords([model], chromosomes,
                                     randWalk=calcParams.randWalk,
                                     randSeed=calcParams.randSeed)
                mCoords = self.getModelCoords(model, chromosomes)
 
              coords.append(mCoords)
 
            coords = array(coords)
 
          # add backbone restraints, get single concatenated restraint arrays
          indices, dists, ambig = concatenateRestraints(restrDict, posDict, seqScale,
                                                        calcParams.bboneLower, calcParams.bboneUpper)
 
          restrIndices.append(indices)
          restrDists.append(dists)
          restrAmbig.append(ambig)
          stageCounts.append(len(ambig))
          posDicts.append(posDict)
 
          # construct annealing schedule for this stage
          temps.append(annealStage.getTempSchedule())
          repulsScales.append(annealStage.getRepulsionSchedule())
          timeSteps.append(annealStage.timeStep)
          dynSteps.append(annealStage.dynSteps)

        restrIndices = concatenate(restrIndices, axis=0)
        restrDists = concatenate(restrDists, axis=0)
        restrAmbig = concatenate(restrAmbig, axis=0)
 
        # run the annealing via wrapper that handles temp schedule, repulsion and parallelisation
        job = self._runAnnealing(coords, posDicts, stageCounts, restrIndices, restrDists, restrAmbig,
                                 temps, repulsScales, timeSteps, dynSteps, callback, bgCalc, numCpus)
      
    # store final restraints
    restraintGroup = self.restraints
    for chrA in restrDict:
      rGroup = self._getGroup(chrA, restraintGroup)
    
      for chrB in restrDict[chrA]:
        dataArray = restrDict[chrA][chrB]
        self._setData(chrB, rGroup, float32, dataArray)
    
    # store seq positions for coords 
    self.addChromosomes(posDict, bboneDict, interpolateCoords=False)
    
    if not bgCalc:
      # store final coords
      coords = job.getResult()
      self.setAllCoords(coords, chromosomes)
 
      # superimpose ensemble
      if len(coords) > 1:
        self.modelAlign(chromosomes=chromosomes)
 
      if callback:
        callback()
 
      self.save()
        
    return job
 
  
  def calcSurface(self, model, chromosome, code=None, sigma=1.4, cubeSize=1.0, level=0.001):
    """Create a volumetric representation of the structural surface.
       A voxellated coordinate list in 3D contour order"""
    
    from util.Structure import coordsToVoxels
    from memops.cNg.Contourer3d import contourer3d
    
    # Inactive/unmappable?
    
    # Nothing formally to say this is a surface
    
    if not code:
      code = 'surface_%s%d' % (chromosome, model)
    
    coords = self.getModelCoords(model, [chromosome,])

    coordsGroup = self._getGroup('coords', self.images)
    chromoGroup = self._getGroup(chromosome, self.structures)
    colors = ones( (len(coords), 4), float)
    
    posVoxels, colorVoxels = coordsToVoxels(coords, colors, cubeSize, sigma, 1.0)
    posVoxels = array(posVoxels, float32)
    #colorArray = array(colorArray, float32)

    meshCoords = contourer3d(posVoxels.T, level, 0.0, 0.0, 0.0, cubeSize, cubeSize, cubeSize)
    #meshCoords = findExtraTriangles(meshCoords)
    
    #self._setData(code, coordsGroup, float, meshCoords)
    
    meshCoords = meshCoords.reshape(len(meshCoords)/3,3)
    self.setVoxelCoords(code, meshCoords)
    
    # Colours not stored at the moment...
    
    
  def calcDepths(self, models, chromosomes):
    """(Re)calculate the depth below the surface for particles in a 3D chromosomes model"""
    
    from util.Structure import coordsToDepths
    
    regionDict = {}
    valueDict = {}
    
    chromoGroup = self.chromosomes
    
    
    for chromo in chromosomes:      
      if chromo not in chromoGroup:
        continue
      
      depthList = []
      
      for model in models:
        coords = self.getModelCoords(model, [chromo,])
        depths = coordsToDepths(coords)
        depthList.append(depths)
        
      depths = array(depthList)
      depths = depths.mean(axis=0)  
            
      positions = array(chromoGroup[chromo]['positions'])
      
      depthsNorm = depths / (depths.max() or 1.0)
      
      valueArray = vstack([depths, depthsNorm]).T
      regionArray = vstack([positions, positions+1.0]).T
    
      regionDict[chromo] = regionArray
      valueDict[chromo] = valueArray
      
      
    if valueDict:
      self.setDataTrack('depth', DERIVED, regionDict, valueDict)
    
  
  def calcDensity(self, model, chromosomes, radius=100.0):
    """(Re)calculate the spatial density of model particles"""
    
    from cUtil.apiUtil import calcCoordDensity
        
    chromoGroup = self.chromosomes
    valueDict = {}
    regionDict = {}
    allCoords = []
    chrRanges = []
    chromos = []
    
    n = 0
    for chromo in chromosomes:
      if chromo not in chromoGroup:
        continue
    
      coords = self.getModelCoords(model, [chromo,])
      nCoords = len(coords)
      chrRanges.append((n, n + nCoords))
      chromos.append(chromo)
      
      n += nCoords
      allCoords.append(coords)
    
    allCoords = vstack(allCoords) 
    allDensities = calcCoordDensity(allCoords, radius)
    
    std = allDensities.std()
    mean = allDensities.mean()
    
    minDensity = mean - 3*std
    maxDensity = mean + 3*std
    
    allDensities = clip(allDensities, minDensity, maxDensity)
    
    for i, chromo in enumerate(chromos):
      start, end = chrRanges[i]
      
      densities = allDensities[start:end]
      
      densNorm = (densities - minDensity) / (maxDensity - minDensity)
      positions = array(chromoGroup[chromo]['positions'])
      
      valueArray = vstack([densities, densNorm]).T
      regionArray = vstack([positions, positions+1.0]).T
    
      regionDict[chromo] = regionArray
      valueDict[chromo] = valueArray
      
    if valueDict:
      self.setDataTrack('density', DERIVED, regionDict, valueDict)
  
  
  def calcRestraintViolations(self, chromosomes, models=None, cis=True,
                              trans=True, upperOnly=False, reportAll=True,
                              usePositions=True):
    """Calculate 3D model distances of directly restrained particles and compare to restraint bounds"""
    
    if not models:
      models = range(self.getNumModels())
    
    restrDict = self.getRestraints(chromosomes, cis, trans, usePositions=usePositions)
    violDict = {}
    
    for key in restrDict: # pos0, pos1, weight, target, lower, upper
      chromoA, chromoB = key
      restraints = restrDict[key].T
      if not len(restraints[0]):
        continue
      
      posA = array(restraints[0], uint32)
      posB = array(restraints[1], uint32)
       
      if usePositions:
        dists = self.getPositionDistances(posA, posB, models, chromoA, chromoB)
        
      else:
        dists = self.getCoordDistances(chromoA, chromoB, models)
      
      targets = array(restraints[3])
      uppers = array(restraints[5])
      
      if upperOnly:
        deltas = dists-uppers
        
        if reportAll:
          deltas = clip(deltas, 0, deltas.max())
          
        else:
          indices = (deltas > 0).nonzero()
          deltas = deltas[indices]
          posA = posA[indices]
          posB = posB[indices]
          targets = targets[indices]
          
      else:
        lowers = array(restraints[4])
 
        deltasUpper = dists-uppers
        deltasLower = lowers-dists
        
        deltasUpper = clip(deltasUpper, 0, deltasUpper.max())
        deltasLower = clip(deltasLower, 0, deltasLower.max())
        
        deltas = deltasUpper
        indices = (deltasLower > deltasUpper).nonzero()
        deltas[indices] = deltasLower[indices]
 
        if not reportAll:
          indices = (deltas > 0).nonzero()
          deltas = deltas[indices]
          posA = posA[indices]
          posB = posB[indices]
          targets = targets[indices]
      
      violDict[key] = vstack([posA, posB, deltas, targets]).T
       
    # (chromoA, chromoB) : [posA, posA, delta]
      
    return violDict
    
    
  def calcRestraintDistances(self, chromosomes, models=None, cis=True, trans=True):
    """Calculate 3D model distances of directly restrained particles"""
    
    if not models:
      models = range(self.getNumModels())

    restrDict = self.getRestraints(chromosomes, cis, trans, usePositions=True)
    distDict = {}
    
    for key in restrDict: # pos0, pos1, weight, target, lower, upper
      chromoA, chromoB = key
      restraints = restrDict[key]
      posA = array(restraints[0], uint32)
      posB = array(restraints[1], uint32)
      
      chrPosA = [(chromoA, pos) for pos in posA]
      chrPosB = [(chromoB, pos) for pos in posB]
    
      dists = self.getPositionDistances(chrPosA, chrPosB, models)
      distDict[key] = dists
      
      
    return distDict
    
  def calcModelRmsds(self, models=None, chromosomes=None, backbone=None, weightThreshold=1.0):
    """Calculate RMSD of structural models as a measure of coordinate precision.
       Backbone True:only, False:exclude, None:everything
       """
   
    from util.Structure import superimposeCoordArray
 
    # Single value for bundle
   
    if not models:
      models = range(self.getNumModels())
 
    if not chromosomes:
      chromosomes = self.getChromosomes()
  
    ensemble = [self.getModelCoords(m, chromosomes, backbone) for m in models]
    ensemble, rotations, rmsds, atomRmsds = superimposeCoordArray(ensemble, weightThreshold)
 
    return rmsds, atomRmsds


  def calcModelRmsdDataTrack(self, chromosome, models=None, backbone=None, weightThreshold=1.0):
    """Calculate RMSD of structural models along the chromosomal sequence
       Backbone True:only, False:exclude, None:everything
       """
       
    from util.Structure import superimposeCoordArray
   
    if not models:
      models = range(self.getNumModels())
   
    chromoGroup = self.chromosomes
    positions = array(chromoGroup[chromosome]['positions'])

    ensemble = [self.getModelCoords(m, [chromosome], backbone) for m in models]
    ensemble, rotations, rmsds, atomRmsds = superimposeCoordArray(ensemble, weightThreshold)
   
    if backbone:
      indices = array(chromoGroup[chromosome]['backbone']).nonzero()
      atomRmsds = atomRmsds[indices]
      positions = positions[indices]
      
    elif backbone is False:
      indices = (chromoGroup[chromosome]['backbone'] < 1).nonzero()
      atomRmsds = atomRmsds[indices]
      positions = positions[indices]

    regionDict = {}
    valueDict = {}
    
    atomRmsdsNorm = atomRmsds / (atomRmsds.max() or 1.0)

    valueArray = vstack([atomRmsds, atomRmsdsNorm]).T
    regionArray = vstack([positions, positions+1.0]).T
    
    regionDict[chromosome] = regionArray
    valueDict[chromosome] = valueArray

    self.setDataTrack('RMSD', DERIVED, regionDict, valueDict)


  def calcModelRmsdMatrix(self, models=None, backbone=None, weightThreshold=1.0):
    """Calculate a pairwise model RMSD matrix
       Backbone True:only, False:exclude, None:everything"""
    
    from util.Structure import superimposeCoordPair
    
    if not models:
      models = range(self.getNumModels())
    
    chromosomes = self.getChromosomes()
    ensemble = [self.getModelCoords(m, chromosomes, backbone) for m in models]
    n = len(models)
   
    matrix = zeros((n,n), float)
    
    for i in range(n-1):
      coords1 = ensemble[i]
    
      for j in range(i+1, n):
        coords2 = ensemble[j]
        coords3, rot, rmsd, atomRmsds = superimposeCoordPair(coords1, coords2, weightThreshold)
        
        matrix[i,j] = rmsd
        matrix[j,i] = rmsd
    
    return matrix
    
    
  def calcModelClusters(self, k=4, backbone=None):
    """Group structural models into a given mumber of clusters"""
    
    from util.Cluster import kMedioids, hierarchicalCluster
    
    rmsdMatrix = matrixOrig = self.calcModelRmsdMatrix(backbone=backbone)
    
    rmsdMatrix /= rmsdMatrix.max() or 1.0
    
    indices = (rmsdMatrix == 0.0).nonzero()
    rmsdMatrix[indices] = 1.0
    
    rmsdMatrix = 1.0 - rmsdMatrix
    
    rmsdMatrix, rows, cols = hierarchicalCluster(rmsdMatrix)
    
    step = int(len(rows)/(k+1))
    centers = [rows[(i+1)] * step for i in range(k)]
    centers, clusters = kMedioids(matrixOrig, centers, k)
    
    representatives = [rows.index(c) for c in centers]
    clusterRmsds = []
    
    for cluster in clusters:
      n = len(cluster)
      rmsds = []
      for i in range(n-1):
        for j in range(i+1, n):
          rmsds.append(matrixOrig[i,j])
    
      clusterRmsds.append(array(rmsds))
    
    meanRmsds = [x.mean() for x in clusterRmsds]
    stdDevRmsds = [x.std(ddof=1) for x in clusterRmsds]
    
    return representatives, meanRmsds, stdDevRmsds
  
  
  def calcBinnedDataTrack(self, inCode, outCode, binSize, start=None, source=EXTERNAL):
    """Create a new data layer by grouping data from another into binned regions"""
    
    dataGroup = self._getDataTrackGroup(source)
    
    from cUtil import dataLayer
    
    if inCode in dataGroup:
      dataLayerA = self.getRefDataTrackGroup(source, inCode)
      dataLayerB = self._getGroup(outCode, dataGroup)
       
      dataLayerB.attrs['stranded'] = dataLayerA.attrs['stranded'] 
      dataLayerB.attrs['options'] = dataLayerA.attrs['options'] 
      dataLayerB.attrs['display'] = dataLayerA.attrs['display'] 
           
      for chromo in dataLayerA:
        regions = array(dataLayerA[chromo]['regions'], uint32) # start, end
        values  = array(dataLayerA[chromo]['values'], float32)[:,0] # origValue
        
        if start is None:
          startPoint = regions.min()
        else:
          startPoint = int32(start)
        
        endPoint = regions.max()
        
        binned = dataLayer.regionBinValues(regions, values, int32(binSize), startPoint, endPoint)
        nBins = len(binned)
        
        positions = binSize * array(range(nBins+1), uint32)
        
        valueArrayB  = vstack([binned, binned/binned.max()]).T
        regionArrayB = vstack([positions[:-1]+1,
                               positions[1:]]).T
       
        grp = self._getGroup(chromo, dataLayerB)
          
        self._setData('regions', grp, uint32, regionArrayB)
        self._setData('values', grp, float32, regionArrayB)
        
        if 'annotations' in dataLayerA[chromo]:
          annotationsB = [''] * nBins
          self._setData('annotations', grp, VL_str, annotationsB)

      return dataLayerB


  def getNumIsolatedContacts(self, groupName, chromosomes=None, threshold=int(2e6),
                             cis=True, trans=True, asFraction=False):
    """Count unsupported structural contacts"""

    
    cacheDict = self.getCachedContacts(groupName)
    
    if not cacheDict:
      if asFraction:
        return 0.0
      else:
        return 0
      
    chromoGroup = self.chromosomes
    unsupIndices = {}

    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    from cUtil.apiUtil import getIsolatedPairs
    t = 0
    n = 0
    
    if cis:
      for chromo in cacheDict:
        if chromo in cacheDict[chromo]:
          positions = array(cacheDict[chromo][chromo][:2], int32).T
          inactive = getIsolatedPairs(positions, int32(threshold))
          t += len(positions)
          n += len(inactive)
 
        
    if trans:
      for chromoA in cacheDict:
        subDict = cacheDict[chromoA]
        
        for chromoB in subDict:
          positions = array(subDict[chromoB][:2], int32).T
          inactive = getIsolatedPairs(positions, int32(threshold))
          t += len(positions)
          n += len(inactive)
          
    if asFraction:
      if t:
        return n/float(t)
      else:
        return 0.0
      
    else:
      return n

    
  def getIsolatedContacts(self, groupName, chromosomes=None, threshold=int(2e6)):
    """Get unsupported structural contacts"""
    
    chromoGroup = self.chromosomes
    cacheDict = self.getCachedContacts(groupName)
    unsupIndices = {}

    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    chromosomes = set(chromosomes)
    
    from cUtil.apiUtil import getIsolatedPairs
     
    for chromoA in cacheDict:
      #if chromoA not in chromosomes:
      #  continue
    
      subDict = cacheDict[chromoA]
        
      for chromoB in subDict:
        #if chromoB not in chromosomes:
        #  continue
         
        positions = array(subDict[chromoB][:2], int32).T
        
        if len(positions):
          inactive = getIsolatedPairs(positions, int32(threshold))
          unsupIndices[(chromoA, chromoB)] = inactive    
    
    return unsupIndices
          
  
  def modelCentre(self, models=None, chromosomes=None):
    """Centre chromosome model coordinates at the orgin"""
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    if not models:
      models = range(self.getNumModels())
    
    centers = []
    
    for model in models:
      coords = self.getModelCoords(model, chromosomes)
      center = coords.sum(axis=0)/float(len(coords))
      coords -= center
      
      self.setModelCoords(coords, model, chromosomes)
      centers.append(center)
      
    return centers
    
    
  def modelAlign(self, models=None, chromosomes=None, backbone=None, weightThreshold=10.0):
    """Superpose structural chromosome models using RMSD weighted iterative SVD
           Backbone True:only, False:exclude, None:everything"""

    from util.Structure import superimposeCoordArray
      
    if not models:
      models = range(self.getNumModels())
 
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    ensemble = [self.getModelCoords(m, chromosomes, backbone) for m in models]
    
    ensemble, rotations, rmsds, atomRmsds = superimposeCoordArray(ensemble, weightThreshold)
    
    for i, modelCoords in enumerate(ensemble):
      self.setModelCoords(modelCoords, models[i],  chromosomes)


  def modelAlignAxes(self, models=None, chromosomes=None):
    """Align largest orthogonal model dimensions with coordinate axes"""

    from util.Cluster import principleComponentAnalysis
    
    if not models:
      models = range(self.getNumModels())
      
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    strucGroup = self.structures
    
    if 'coords' in strucGroup:
      for model in models:
        modelCoords = self.getModelCoords(model, chromosomes)
        eigenMat, energy = principleComponentAnalysis(modelCoords)
        modelCoords = dot(modelCoords, eigenMat)
        self.setModelCoords(modelCoords, model, chromosomes)


  def modelScale(self, factor, models=None, chromosomes=None):
    """Scale 3D chromosome model coordinates by specified factor"""
  
    matrix = factor * eye(3)
    
    self.modelTransform(matrix, models, chromosomes)
    
    
  def modelRotate(self, angle, axis=(0,0,1), models=None, chromosomes=None):
    """Rootate chromosome model coordinates by specified angle around an axis"""
    
    from util.Structure import getRotationMatrix
     
    matrix = getRotationMatrix(axis, angle)
    
    self.modelTransform(matrix, models, chromosomes)
    
    
  def modelTranslate(self, vector, models=None, chromosomes=None):
  
    vector = array(vector) 
     
    for model in models:
      modelCoords = self.getModelCoords(model, chromosomes)
      modelCoords += vector
      self.setModelCoords(modelCoords, model, chromosomes)
     
     
  def modelMirror(self, models=None, chromosomes=None):
    """Make 3D chromosome model coordinates a mirror image"""

    matrix = eye(3)
    matrix[0,0] = -1
    
    self.modelTransform(matrix, models, chromosomes)
    
    
  def modelTransform(self, matrix, models=None, chromosomes=None):
    """Apply a 3D matrixs tranformation to the coordinates of chromosomes in structural models"""

    if not models:
      models = range(self.getNumModels())
      
    if not chromosomes:
      chromosomes = self.getChromosomes()
   
    matrix = array(matrix)

    if matrix.shape == (3,3):
      for model in models:
        modelCoords = self.getModelCoords(model, chromosomes)
        modelCoords = dot(modelCoords, matrix)
        self.setModelCoords(modelCoords, model, chromosomes)

    else:
      matShape = 'x'.join(['%d' % x for x in matrix.shape])
      raise Exception('Transformation matrix must be 3x3 not %s' % matShape)
  
  
  def setSpiralCoords(self, radius=10.0, chromosomes=None):
    """Set stucture to a spriral, given chromosomal positions"""
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
      
    chromoGroup = self.chromosomes
    nPoints = sum([len(chromoGroup[c]['positions']) for c in chromosomes])
    dAngle = 25*TAU/float(nPoints or 1.0)
    dZ = 2 * radius/float(nPoints or 1.0)
    
    modelCoords = [None] * nPoints
    angle = 0.0
    z = 0.0
    
    for i in range(nPoints):
      x = radius * cos(angle)
      y = radius * sin(angle)
      
      modelCoords[i] = (x,y,z)
      angle += dAngle
      angle = angle % TAU
      z += dZ
    
    if modelCoords:
      modelCoords = array([modelCoords,])
      self.setAllCoords(modelCoords, chromosomes)
    
  
  def setRandomCoords(self, models=None, chromosomes=None, centre=(0.0,0.0,0.0),
                      randWalk=True, randSeed=None, maxStep=1.0):
    """Set the coordinates of chromsome structural models to random positions"""
    
    uniform = random.uniform
    
    if randSeed:
      seed(randSeed)
    else:
      seed()
    
    if not models:
      models = range(self.getNumModels())
      
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    
    centre = array(centre)
    chromoGroup = self.chromosomes
    nPoints = sum([len(chromoGroup[c]['positions']) for c in chromosomes])
    allCoords = []
    
    for model in models:
      modelCoords = [None] * nPoints
        
      if randWalk:
        x = 0.0
        y = 0.0
        z = 0.0
        phi = 0.0
        psi = 0.0
        
        for i in range(nPoints):
          phi += uniform() * TAU
          psi += uniform() * TAU
 
          phi = phi % TAU
          psi = psi % TAU
 
          x += maxStep * sin(phi) * sin(psi)
          y += maxStep * sin(phi) * cos(psi)
          z += maxStep * cos(phi)
 
          modelCoords[i] = [x,y,z]
      
      else:
        limit = maxStep*maxStep
        
        for i in range(nPoints):
          x = y = z = maxStep
 
          while x*x + y*y + z*z >= limit:
            x = maxStep * (2*uniform(0,1) - 1)
            y = maxStep * (2*uniform(0,1) - 1)
            z = maxStep * (2*uniform(0,1) - 1)
          
          modelCoords[i] = [x,y,z]
      
      modelCoords = centre + array(modelCoords)
      allCoords.append(modelCoords)
    
    allCoords = array(allCoords)  
    self.setAllCoords(allCoords, chromosomes)
  
  
  def setRandomContacts(self, groupName, model, numRestraints, chromosomes=None,
                        numNeighbours=10, replace=True):
    """Make restraints based upon random close points, 
       e.g. given a synthetic structure"""
    
    from cUtil.apiUtil import getClosestPoints
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    chromoGroup = self.chromosomes
    coords = self.getModelCoords(model, chromosomes)
    
    start = 0
    chromoStarts = []
    posDict = {}
    
    for chromo in chromosomes:
      chromoStarts.append(start)
      pos = array(chromoGroup[chromo]['positions'])
      posDict[chromo] = pos
      start += len(pos)
    
    n = 0
    contactDict = {}
    append = {}
    numCoords = len(coords)
    done = set()
    doneAdd = done.add
    numObs = 1
    numRestraints = min(numCoords/2, numRestraints)
    nn = int32(numNeighbours)
    
    while n < numRestraints:
      i = randint(0, numCoords-1)
      coord = coords[i]

      indices = getClosestPoints(coord, coords, nn)

      k = randint(0, numNeighbours-1)
      j = indices[k]
      
      key = frozenset((j,i))
      if key in done:
        k = randint(0, numNeighbours-1) # Second chance
        j = indices[k] 
        key = frozenset((j,i))

        if key in done:
          continue # Could increase numObs 
      
      for iChr, start in enumerate(chromoStarts):
        if start > i:
          iChr -= 1
          break
 
      for jChr, start in enumerate(chromoStarts):
        if start > j:
          jChr -= 1
          break
      
      chrA = chromosomes[iChr]
      chrB = chromosomes[jChr]
 
      posA = posDict[chrA]
      posB = posDict[chrB]
      
      idxA = i-chromoStarts[iChr]
      idxB = j-chromoStarts[jChr]
             
      if chrB < chrA:
        chrKey = (chrB, chrA)
        datum = [posB[idxB], posA[idxA], numObs, model]
      
      else:
        chrKey = (chrA, chrB)
        datum = [posA[idxA], posB[idxB], numObs, model]
      
      if chrKey in contactDict:
        append[chrKey](datum)
      else:
        contactDict[chrKey] = [datum]
        append[chrKey] = contactDict[chrKey].append
        
      doneAdd(key)
      
      n += 1
    
    self.setContacts(groupName, contactDict, replace)
    
    return n
    
    
  def setTestCoordsCoil(self, chromosomes=None, numRestraints=1024, ppTurn1=30, ppTurn2=7,
                        rad1=0.4, rad2=0.2, contact=0.06, scale=45.0):
                    
    """Set the coordinates of chromsome structural models to a synthetic test structure"""
    
    from util.Structure import getCoiledCoilCoords, HEX_GRID
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    chromoGroup = self.chromosomes
    model = 0
    
    seed()
 
    shuffle(chromosomes)

    spacing = 2.0 * (rad1 + rad2) + contact
    hexGrid = HEX_GRID * spacing
 
    dAngle1=TAU/ppTurn1
    dAngle2=TAU/ppTurn2
 
    allCoords = []
 
    nodes = []
    for i, chromo in enumerate(chromosomes):
      numPoints = len(chromoGroup[chromo]['positions'])
    
      x0, y0 = hexGrid[i]
      baseCoord = [x0, y0, 0.0]
      coords = getCoiledCoilCoords(baseCoord, numPoints, dAngle1,
                                   dAngle2, rad1, rad2, contact)
      coords *= scale
      
      self.setModelChromosomeCoords(coords, chromo, model)
 
    self.setRandomContacts('singleCell', model, numRestraints, chromosomes)
  
  
  def setTestCoordsHilbert(self, chromosomes=None, centre=(0.0,0.0,0.0), scale=5.0):
    """Set the coordinates of chromsome structural models to a synthetic test structure"""
    
    from util.Structure import indexToHilbert3D
    
    chromoGroup = self.chromosomes
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    model = 0
    chromoSizes = [len(chromoGroup[c]['positions']) for c in chromosomes]
    nPoints = sum(chromoSizes)
    
    coords = array([indexToHilbert3D(i) for i in range(nPoints)]) * scale
    coords -= coords.mean(axis=0) - array(centre) 
 
    self.setAllCoords(coords, chromosomes)
    self.setRandomContacts('singleCell', model, chromosomes)
  
  
  def setContacts(self, groupName, contactDict, replace=True,
                  selected=None, color=None):
    """Sets working contacts, never originals; they are only loaded
        By default any previous contacts will be replaced
    """
    
    group = self._getContactGroup(groupName, self.workContacts)
    
    if color is not None:
      group.attrs['color'] = color
        
    if selected is not None:
      group.attrs['selected'] = selected
      
    for chromoKey in contactDict:
      chromoA, chromoB = chromoKey
      contacts = array(contactDict[chromoKey]).T
      subGroup = self._getGroup(chromoA, group)
      
      if (chromoB in subGroup) and not replace:
        prev = array(subGroup[chromoB])
        prevDim = len(prev)
        contDim = len(contacts) 
        
        if prevDim < contDim:
          msg  = 'Attempt to add contacts with model indices to a dataset '
          msg += 'without model indices. Model indices will be removed.'
          self._warning(msg)
          contacts = contacts[:3]
          
        elif contDim < prevDim:
          msg =  'Attempt to add contacts without model indices to a dataset '
          msg += 'which requires them'
          raise Exception(msg)
          
        contacts = hstack([prev, contacts])
        
      ds = self._setData(chromoB, subGroup, uint32, contacts, compression='gzip')
    
    if groupName in self._contactsCache:
      del self._contactsCache[groupName]
  
  
  def revertContacts(self, groupName):
    """Revert a contact group back to the orginal, loaded data"""
    
    sourceGroup = self._getContactGroup(groupName, self.origContacts)
    
    if groupName in self.workContacts:
      del self.workContacts[groupName]
      group = self._getGroup(groupName, self.workContacts)
    
      for chrA in sourceGroup:
        subGroup = self._getGroup(chrA, group)
      
        for chrB in sourceGroup[chrA]:
          data = array(sourceGroup[chrA][chrB])
          self._setData(chrB, subGroup, uint32, data, compression='gzip') 

  
  def removeContacts(self, groupName):

    if groupName in self.workContacts:
      self._delete(groupName, self.workContacts)
    
    if groupName in self.origContacts:
      self._delete(groupName, self.origContacts)
   
  
  def removeViolatedContacts(self, groupName, chromosomes=None, threshold=4.0, models=None):
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    n = 0
    if models and len(models) == 1:
      model = models[0]
    else:
      model = None
    
    contDict = self.getContacts(groupName, chromosomes, cis=True, trans=True, model=model)
    group = None
    
    for chromoKey in contDict:
      chromoA, chromoB = chromoKey
      data = contDict[chromoKey]
      
      posA = data[0]
      posB = data[1]
      
      chrPosA = [(chromoA, pos) for pos in posA]
      chrPosB = [(chromoB, pos) for pos in posB]
      dists = self.getPositionDistances(chrPosA, chrPosB, models)
      meanDist = dists.mean()
      
      indices = (dists <= meanDist*threshold).nonzero()
      
      if len(indices):
        n +=  len(posA) - indices[0].shape[0]
        data = data[:,indices[0]]
 
        if not group:
          group = self._getContactGroup(groupName, self.workContacts)
 
        subGroup = self._getGroup(chromoA, group)
        self._setData(chromoB, subGroup, uint32, data, compression='gzip')
    
    if groupName in self._contactsCache:
      del self._contactsCache[groupName]
    
    return n  
    
    
  def removeIsolatedContacts(self, groupName, chromosomes=None, threshold=int(2e6)):
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    cacheDict = self.getCachedContacts(groupName)
    unsupIndices = self.getIsolatedContacts(groupName, chromosomes, threshold) or {}
    
    n = 0
    
    if unsupIndices:
      group = self._getContactGroup(groupName, self.workContacts)
    
      for chromoKey in unsupIndices:
        chromoA, chromoB = chromoKey
        inactive = unsupIndices[chromoKey]
        n += len(inactive)
 
        contacts = cacheDict[chromoA][chromoB]
        
        selected = ones(contacts.shape[1], uint32)
        selected[inactive] = 0
        
        indices = selected.nonzero()[0]
        contacts = contacts[:,indices]
        
        subGroup = self._getGroup(chromoA, group)
        self._setData(chromoB, subGroup, uint32, contacts, compression='gzip')
    
    if groupName in self._contactsCache:
      del self._contactsCache[groupName]
      
    self.getCachedContacts(groupName)
       
    return n
    
    
  def setModelChromosomeCoords(self, coords, chromosome, model=None):
    """Set all the 3D coordinates of one chromosome.
       Num models comes from size fo coords array if model not set"""
   
    strucGroup = self.structures
    chromoGroup = self.chromosomes
    allChromos = self.getChromosomes()
    nPoints = len(chromoGroup[chromosome]['positions'])
    
    if model is None:
      if coords.ndim == 2:
        n = len(coords)
        coords = coords.reshape(1, n, 3)
 
      else:
        n = len(coords[0])
 
    elif coords.ndim == 3:
      if len(coords) == 1:
        coords = coords[0]
        
      else:
        coords = coords[model]
      
      n = len(coords)
      
    else:
      n = len(coords)
      
    if n != nPoints:
      msg = 'Model coordinates must be an array of length %d or have shape numModels x %d not "%d"'
      raise(Exception(msg % (nPoints, nPoints, n)))
       
    coordsGroup = self._getGroup('coords', strucGroup)

    if model is None:
      self._setData(chromosome, coordsGroup, float, coords)
    
    elif chromosome in coordsGroup:
      allCoords = array(coordsGroup[chromosome])
      
      if model < len(allCoords):
        allCoords[model] = coords
      
      else:
        while model >= len(allCoords):
          allCoords = vstack([allCoords, coords])
          
      self._setData(chromosome, coordsGroup, float, allCoords)
      
    else:
      allCoords =  vstack([coords] * (model+1))
      self._setData(chromosome, coordsGroup, float, allCoords)
      
      
  def setAllCoords(self, coords, chromosomes=None):
    """Set all the 3D model coordinates of specified chromosomes.
       Coords must be in chromosome order"""

    strucGroup = self.structures
    chromoGroup = self.chromosomes
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    elif not isinstance(chromosomes, (tuple, list, ndarray)):
      chromosomes = [chromosomes,]
    
    
    chromoSizes = [len(chromoGroup[c]['positions']) for c in chromosomes]
    nPoints = sum(chromoSizes)
    
    if coords.ndim == 2:
      coords = array([coords,])
    
    m, n, dims = coords.shape

    if n != nPoints:
      msg = 'Model coordinates must be an array of numModels x %d' % (nPoints,)
      raise(Exception(msg))
    
    coordsGroup = self._getGroup('coords', strucGroup)
    
    j = 0
    for i, chromo in enumerate(chromosomes):
      span = chromoSizes[i]
      chromoCoords = coords[:,j:j+span] # input models
      self._setData(chromo, coordsGroup, float, chromoCoords)
      
      j += span
    
         
  def setModelCoords(self, coords, model, chromosomes=None):
    """Set all the 3D coordinates for one model of selected chromosomes.
       Corods must be in chromosome order"""
    
    strucGroup = self.structures
    chromoGroup = self.chromosomes

    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    elif not isinstance(chromosomes, (tuple, list, ndarray)):
      chromosomes = [chromosomes,]
    
    else:
      chromosomes = sorted(chromosomes)     
    
    if coords.ndim != 2:
      msg = 'Model coordinates must be an array of rank 2'
      raise(Exception(msg))
      
    
    chromoSizes = [len(chromoGroup[c]['positions']) for c in chromosomes]
    nPoints = sum(chromoSizes)
     
    n, dims = coords.shape

    if n != nPoints:
      msg = 'Model coordinates must be an array of length %d' % (nPoints,)
      raise(Exception(msg))
    
    coordsGroup = self._getGroup('coords', strucGroup)
    
    j = 0
    for i, chromo in enumerate(chromosomes):
      span = chromoSizes[i]
      
      if (chromo in coordsGroup) and len(coordsGroup[chromo].shape) == 3:
        nModels, nPointsM, dims = coordsGroup[chromo].shape
        nPoints = chromoSizes[i]
         
        if nPointsM == span:
          if model < nModels:
            coordsGroup[chromo][model] = coords[j:j+span]
 
          else:
            allCoords = array(coordsGroup[chromo])
            
            while model >= len(allCoords):
              allCoords = append(allCoords, coords[j:j+span])
              
            self._setData(chromo, coordsGroup, float, allCoords)
        
        else:# Points have changed, wipe out old coords
          chromoCoords = coords[j:j+span]
          chromoCoords = array([chromoCoords] * (model+1)) # Fill implicit models
          self._setData(chromo, coordsGroup, float, chromoCoords)
     
      else:
        chromoCoords = coords[j:j+span]
        chromoCoords = array([chromoCoords] * (model+1)) # Fill implicit models
        self._setData(chromo, coordsGroup, float, chromoCoords)
      
      j += span
       
    
  def setSampleInfo(self, infoDict):
    """Set information relating to the experimental sample"""

    for key, value in infoDict.iteritems():
      self.sample.attrs[key] = value

    
  def addChromosomes(self, positionsDict, backboneDict=None, interpolateCoords=True):
    """Add new chromosome/DNA region identities"""
    
    
    from cUtil.apiUtil import interpolateChromoModelCoords
        
    strucGroup = self.structures
    coordsGroup = self._getGroup('coords', strucGroup)
    models = list(range(self.getNumModels()))
    prevPosDict = {}
    
    
    for chromo in positionsDict:
      positions = positionsDict[chromo]
      positions = array(positions, uint32)
      n = len(positions)
      
      if backboneDict and (chromo in backboneDict):
        backbone = array(backboneDict[chromo], uint32)
      else:  
        backbone = zeros(n, uint32)
     
      mapability = ones(n, float32)
      
      chromoGroup = self._getGroup(chromo, self.chromosomes)
      coordArray = None
      
      if interpolateCoords and ('positions' in chromoGroup):
        prevPosDict[chromo] = array(chromoGroup['positions'], int)
       
      self._setData('backbone', chromoGroup, uint32, backbone)
      self._setData('positions', chromoGroup, float32, positions)
      self._setData('mapability', chromoGroup, float32, mapability)
      
      self._chromoLimitCache[chromo] = (positions[0], positions[-1])
    
    nModels = self.getNumModels()
    
    for chromo in prevPosDict:
      pDictA = {chromo: array(positionsDict[chromo], int32)}
      pDictB = {chromo: array(prevPosDict[chromo], int32)} 
    
      for model in range(nModels):
        coords = array(coordsGroup[chromo][model])
        coords = interpolateChromoModelCoords(pDictA, pDictB, coords)
        self.setModelCoords(coords, model, chromosomes=[chromo,])

    
    names = set(self.chromosomes)
    names.update(positionsDict.keys())
    names = list(names)
    names.sort(key=naturalKey())
    self._chromosomes = names
        
    if not nModels:
      self.setSpiralCoords()
      
    
  def removeChromosomes(self, names):
    """Remove existing chromosome/DNA region identities and directly associated data"""
    
    for chromo in names:
      self._delete(chromo, self.chromosomes)
    
    for mainGroup in (self.externalData, self.innateData, self.derivedData):
      for dataLayer in mainGroup:
        dataGroup = mainGroup[dataLayer]
 
        for code in dataGroup:
          codeGroup = dataGroup[code]
 
          for chromo in names:
            if chromo in codeGroup:
              self._delete(chromo, codeGroup)
       
    if 'coords' in self.structures:
      coordGroup = self.structures['coords']
    
      for chromo in names:
        if chromo in coordGroup:
          self._delete(chromo, coordGroup)
    
    chromos = [x for x in self.restraints]
    for chromoA in chromos:
      if chromoA in names:
        self._delete(chromoA, self.restraints)
    
    for chromoA in self.restraints:
      chromos = [x for x in self.restraints[chromoA]]
      
      for chromoB in chromos:
        if chromoB in names:
          self._delete(chromoB, self.restraints[chromoA])      
      
    for mainGroup in (self.workContacts, self.origContacts):
      for groupName in mainGroup:
        group = mainGroup[groupName]
        chrsA = [x for x in group]
        
        for chrA in chrsA:
          if chrA in names:
            self._delete(chrA, group)
          
          chrsB = [x for x in group[chrA]]
          for chrB in chrsB:
            if chrB in names:
              self._delete(chrB, group[chrA])
              
        if groupName in self._contactsCache:
          for chromo in names:
            if chromo in self._contactsCache[groupName]:
              del self._contactsCache[groupName][chromo]
          
          for chrA in self._contactsCache[groupName]:
            for chrB in names:
              if chrB in self._contactsCache[groupName][chrA]:
                del self._contactsCache[groupName][chrA][chrB]
          
    for chromo in names:
      if chromo in self._chromoLimitCache:
        del self._chromoLimitCache[chromo]
       
      if chromo in self._chromosomes:
        self._chromosomes.remove(chromo)    
         
     
  def setChromosomes(self, names, limits, bbSep=None):
    """Set entirely new chromosome/DNA region identities.
       Destroys existing chromsome data, e.g. coords, dataLayers"""
    
    current = self.getChromosomes()
    self.removeChromosomes(current)
    
    positionsDict = {}
    for i, name in enumerate(names):
      start, end = limits[i]
    
      if bbSep:
        nPos = int(end/float(bbSep)) + 1
        nPos -= int(start)
        start = bbSep * int(start)
        positions = [start + (i*bbSep) for i in range(nPos)]
      else:
        positions = [start, end]
        
      positionsDict[name] = positions
    
    self.addChromosomes(positionsDict)    
    
    
  def setDomains(self, chromosome, regions):
    """Set topological domain geions"""
    
    # TBD: Check regions are within bounds
    
    domainGroup = self._getGroup('domains', self.derivedData)
    regions = array(regions, uint32) # [(start1, end1), (start2, end2), ...]
    
    self._setData(chromosome, domainGroup, float, regions)
    
  
  def setRestraints(self, chromosomes=None, groupName='singleCell', bboneSep=int(5e4),
                    binned=True, scale=1.0, exponent=-2, lower=0.8, upper=1.2,
                    minCount=2, maxPopDist=5.0, model=None):
    """Setup backbone and contact restraints accoring to specified parameters"""
    
    from cUtil.apiUtil import calcRestraints
    
    strucGroup = self.structures
    chromoGroup = self.chromosomes
    
    if 'restraints' in strucGroup:
      del strucGroup['restraints']
    
    self.restraints = restraintGroup = self._getGroup('restraints', strucGroup)
    contactDict = self.getCachedContacts(groupName)
        
    if chromosomes:
      chromosomes = set(chromosomes)
    
    else:
      chromosomes = set(self.getChromosomes())
    
    if groupName not in self.origContacts:
      raise Exception('Contact group name "%s" not present in .nuc file' % groupName)
    
    isSingleCell = int32(self.origContacts[groupName].attrs['isSingleCell'])
    binned = 1 if binned else 0
    model = -1 if model is None else model
    
    restraintsDict, posDict, bboneDict = calcRestraints(chromosomes, contactDict, isSingleCell, int32(bboneSep),
                                                        int32(binned), scale, exponent, lower, upper,
                                                        int32(minCount), maxPopDist, int32(model))

    for chrA in restraintsDict:
      rGroup = self._getGroup(chrA, restraintGroup)
    
      for chrB in restraintsDict[chrA]:
        dataArray = restraintsDict[chrA][chrB]
        self._setData(chrB, rGroup, float32, dataArray)

    self.addChromosomes(posDict, bboneDict) 
 
  
  def setBackboneSpacing(self, bbSep=int(5e4), offset=None):
    """Set the chromosomal sequence locations of backbone model particles"""
    
    chromos = self.getChromosomes()
    limits = [self.getChromosomeLimits(c) for c in chromos]
    
    if offset is not None:
      for i, lim in enumerate(limits):
        lim = list(lim)
        lim[0] = offset
        limits[i] = lim
    
    self.setChromosomes(chromos, limits, bbSep)
    

  def setRefMapability(self, chromosome, values):
    """Set the sequence mapability for organism reference"""
    
    chromoGroup = self.chromosomes[chromosome]
    dataArray = array(values, float)
    
    if dataArray.shape != chromoGroup['positions'].shape:
      s1 = ','.join(['%d' % x for x in dataArray.shape])
      s2 = ','.join(['%d' % x for x in chromoGroup['positions'].shape])
      raise Exception('Mapability data position size mismatch %s != %s' % (s1, s2))
    
    self._setData('mapability', chromoGroup, float, dataArray)
  
  
  def setMapability(self, chromosomes=None):
    """Set the sequence mapability for particle positions"""
    
    if not self.genomeRef:
      return
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    for chromosome in chromosomes:
      # Depends on genome build reference of file specified
      posMaps = self.genomeRef.getMapability(chromosome)
      posMaps = posMaps.tolist().reverse()
 
      chromoGroup = self.chromosomes[chromosome]
      positions = array(chromoGroup['positions'])
 
      # the average mapability of all the values for ref points nearest to this one
      # or at least the closest
 
      valueLists = [[] for i in range(len(positions))]
 
      for i, posI in enumerate(positions[:-1]):
        j = i + 1
        posJ = positions[j]
 
        while posMaps:
          posRef, mapRef = posMaps.pop()
          deltaI = abs(posI-posRef)
          deltaJ = abs(posJ-posRef)
 
          if deltaJ < deltaI:
            if not valueLists[i]:
              valueLists[i].append(mapRef) # At least nearest
 
            valueLists[j].append(mapRef)
            break
 
          else:
            valueLists[i].append(mapRef)
 
      for i, values in enumerate(valueLists):
        if values:
          valueLists[i] = sum(values)/float(len(values))
        else:
          valueLists[i] = 0.0
 
      dataArray = array(valueLists)
 
      self._setData('mapability', chromoGroup, float, dataArray)
  
  
  def setStructureView(self, viewPos):
    """Set 3D viewer global positioning"""
    
    viewPos = array(viewPos, float)
    
    viewArray = self.display.attrs['view3d']
    viewArray[:3] = viewPos
    
    self.root.attrs['view'] = viewArray    
  
  
  def setGlobalTransform(self, matrix):
    """Set affine transformation associated with all models"""
    
    tformGroup = self._getGroup('transforms', self.structures)
    
    self._setData('global', tformGroup, float, array(matrix))
    
     
  def setModelTransform(self, model, matrix):
    """Set affine transformation associated with structure.
       This doesn't change the stored coordinates.
       This is for display and comparison only."""
    
    tformGroup = self._getGroup('transforms', self.structures)
    
    if 'models' in tformGroup:
      transforms = array(tformGroup['models'])
      
      while len(transforms) <= model:
        append(transforms, eye(3), axis=0)
    
    else:
      numModels = self.getNumModels()
      transforms = array([eye(3)] * numModels)
    
    transforms[model] = array(matrix)

    self._setData('models', tformGroup, float, transforms)


  def setDataTrackList(self, source, code, dataList, haveEnds=None, fixedValue=None):
    """Set data values for an existing data layer, or make a new one should none exist.
       Using simple list of (chr, start, end, value) etc"""
    
    # Auto-normalisation is per chromosome
    
    dataList.sort()
    normalise = False
    regionDict = {}
    valueDict = {}
    annoDict = {}
    
    chromos = set(['%s' % x[0] for x in dataList])
    for chromo in chromos:
      regionDict[chromo] = []
      valueDict[chromo] = []
      annoDict[chromo] = []
      
    numTypes = tuple([isinstance(x, (float, int)) for x in dataList[0]])
    haveEnds = type(dataList[0][2]) is int

    if numTypes == (False, True, True, True, True, False):
      for chromo, start, end, origVal, normVal, anno in dataList:
        regionDict[chromo].append((start, end))
        valueDict[chromo].append((origVal, normVal))
        annoDict[chromo].append(anno)
        
    elif numTypes == (False, True, True, False, False, True):
      for chromo, start, end, anno, null, origVal in dataList:
        regionDict[chromo].append((start, end))
        valueDict[chromo].append((origVal, origVal))
        annoDict[chromo].append(anno)
                
    elif numTypes == (False, True, True, True, True):
      for chromo, start, end, origVal, normVal in dataList:
        regionDict[chromo].append((start, end))
        valueDict[chromo].append((origVal, normVal))
   
    elif numTypes == (False, True, True, True, False):
      if haveEnds:
        normalise = True
        for chromo, start, end, value, anno in dataList:
          regionDict[chromo].append((start, end))
          valueDict[chromo].append(( value, value))
          annoDict[chromo].append(anno)
      
      else:
        for chromo, start, origVal, normVal, anno in dataList:
          regionDict[chromo].append((start, start+1))
          valueDict[chromo].append((origVal, normVal))
          annoDict[chromo].append(anno)

    elif numTypes == (True, True, True, True, False):
      if haveEnds:
        normalise = True
        for chromo, start, end, value, anno in dataList:
          chromo = '%d' % chromo
          regionDict[chromo].append((start, end))
          valueDict[chromo].append(( value, value))
          annoDict[chromo].append(anno)
      
      else:
        for chromo, start, origVal, normVal, anno in dataList:
          chromo = '%d' % chromo
          regionDict[chromo].append((start, start+1))
          valueDict[chromo].append((origVal, normVal))
          annoDict[chromo].append(anno)

    elif numTypes == (False, True, True, True):
      if haveEnds:
        normalise = True
        for chromo, start, end, value in dataList:
          regionDict[chromo].append((start, end))
          valueDict[chromo].append((value, value))
      
      else:
         for chromo, start, origVal, normVal in dataList:
          regionDict[chromo].append((start, start+1))
          valueDict[chromo].append((origVal, normVal))
     
    elif numTypes == (False, True, True, False):
      if haveEnds:
        for chromo, start, end, anno in dataList:
          regionDict[chromo].append((start, end))
          valueDict[chromo].append((1.0, 1.0))
          annoDict[chromo].append(anno)
    
      else:
        normalise = True
        for chromo, start, value, anno in dataList:
          regionDict[chromo].append((start, start+1))
          valueDict[chromo].append((value, value))
          annoDict[chromo].append(anno)
    
    elif numTypes == (False, True, True):
      normalise = True
      
      if haveEnds:
        for chromo, start, end in dataList:
          regionDict[chromo].append((start, end))
          valueDict[chromo].append((1.0, 1.0))
      
      else:
        for chromo, start, value in dataList:
          regionDict[chromo].append((start, start+1))
          valueDict[chromo].append((value, value))

    elif numTypes == (False, True, False):
      for chromo, start, anno in dataList:
        regionDict[chromo].append((start, start+1))
        valueDict[chromo].append((1.0, 1.0))
        annoDict[chromo].append(anno)

    elif numTypes == (False, True):
      for chromo, start in dataList:
        regionDict[chromo].append((start, start+1))
        valueDict[chromo].append((1.0, 1.0))
    
    else:
      typeStr = ', '.join([str(type(x)) for x in dataList[0]])
      raise Exception('Data types (%s) not understood in data layer input' % (typeStr,))
    
    stranded = False
    for chromo in valueDict:
      valArray = array(valueDict[chromo], float)
      regionArray = array(regionDict[chromo], int)
       
      if normalise:
        maxVal = valArray[:,1].max()
        if maxVal:
          valArray[:,1] /= maxVal
      
      if not stranded:
        widths = regionArray[:,1] - regionArray[:,0]
        
        if (widths < 0).nonzero():
          stranded = 1
      
      valueDict[chromo] = valArray
      regionDict[chromo] = regionArray
      
      if annoDict[chromo]:
        annoArray = array(annoDict[chromo], dtype='O')
        annoDict[chromo] = annoArray
    
    self.setDataTrack(code, source, regionDict, valueDict, annoDict, stranded)
  
  
  def getDataTrackMaxValues(self, source, code, recalc=False):
    """Fetch the maximum original and normalised data track values
       using cached values where possible"""
    

    dataGroup = self._getDataTrackGroup(source)
    
    if code in dataGroup:
      refLayer = self.getRefDataTrackGroup(source, code)
      dataLayer = dataGroup[code]
      maxOrig, maxNorm = dataLayer.attrs['display'][6:8]
    
      if recalc or not (maxOrig and maxNorm):
        maxVals = []
 
        for chromo in refLayer:
          maxVals.append( array(refLayer[chromo]['values']).max(axis=0) )
 
        maxVals = array(maxVals).max(axis=0)
        maxOrig, maxNorm = maxVals
 
        dispVals = list(dataLayer.attrs['display'])
        dispVals[6] = maxOrig
        dispVals[7] = maxNorm
        dataLayer.attrs['display'] = dispVals
      
      return maxOrig, maxNorm
      
    
  def setDataTrack(self, code, source, regionDict, valueDict, annoDict=None, stranded=None,
                    modelDict=None, color=(1.0, 0.0, 0.0), scale=1.0, threshold=0.0,
                    showText=True, shape=0):
    """Set data values for an existing data layer, or make a new one should none exist.
       Regions, values and annotations set using chromsome keyed dictionaries."""
  
    dataGroup = self._getDataTrackGroup(source)
    annoDict = annoDict or {}
    r, g, b = color or (1.0, 0.0, 0.0)
     
    if code in dataGroup:
      dataLayer = dataGroup[code]
        
      if stranded is not None:
        dataLayer.attrs['stranded'] = 1 if stranded else 0
      
      dispAttrs = list(dataLayer.attrs['display'])
      dispAttrs[6:8] = [0.0, 0.0]
      dataLayer.attrs['display'] = dispAttrs
      
    else:
      dataLayer = self._getGroup(code, dataGroup)
      stranded = 1 if stranded else 0
      showText = 1 if showText else 0
      
      dataLayer.attrs['stranded'] = stranded
      dataLayer.attrs['options'] = (0, showText, shape, 0, 0, 0, 0, 0)
      dataLayer.attrs['display'] = (r, g, b, 1.0, scale, threshold, 0.0, 0.0)

    for chromo in regionDict:
      regionArray = regionDict[chromo]
      n = len(regionArray)
      
      valueArray  = valueDict[chromo]
      if len(valueArray) != n:
        msg = 'Data track size mismatch: values (%d) must match number of regions (%d)'
        raise Exception(msg % (len(valueArray), n))
      
      subGroup = self._getGroup(chromo, dataLayer)
      self._setData('regions', subGroup, uint32, regionArray)
      self._setData('values', subGroup, float32, valueArray)
      
      if annoDict and chromo in annoDict:
        annoArray = annoDict[chromo]
        
        if len(annoArray):
          if len(annoArray) != n:
            msg = 'Data track size mismatch: annotations (%d) must match number of regions (%d)'
            raise Exception(msg % (len(annoArray), n))
          
          self._setData('annotations', subGroup, VL_str, annoArray)
      
      if modelDict and chromo in modelDict:
        modelArray = modelDict[chromo]
        
        if len(modelArray):
          if len(modelArray) != n:
            msg = 'Data track size mismatch: models (%d) must match number of regions (%d)'
            raise Exception(msg % (len(modelArray), n))
          
          self._setData('models', subGroup, uint32, modelArray)
      
      key = (source, code, chromo) 
      if key in self._dataTrackHistogramCache:
        del self._dataTrackHistogramCache[key]
   
    self._setDataTrackIndicesMb(source, code)      
  
  
  def setChromosomeData(self, chromo, code, source, regions, values,
                        annotations=None, models=None):
  
    dataGroup = self._getDataTrackGroup(source)
    dataLayer = self._getGroup(code, dataGroup)
    chromoGroup = self._getGroup(chromo, dataLayer)
    
    n = len(regions)
    
    if len(values) != n:
      msg = 'Data track size mismatch: values (%d) must match number of regions (%d)'
      raise Exception(msg % (len(values), n))
      
    regions = array(regions)
    if regions.ndim == 1:
      regions = vstack([regions, regions+1]).T
    
    values = array(values)
    if values.ndim == 1:
      values = vstack([values, values]).T
    
    self._setData('regions', chromoGroup, uint32, regions)
    self._setData('values', chromoGroup, float32, values)
    
    if annotations is not None:
      if len(annotations) != n:
        msg = 'Data track size mismatch: annotations (%d) must match number of regions (%d)'
        raise Exception(msg % (len(annotations), n))
    
      self._setData('annotations', chromoGroup, VL_str, annotations)
    
    if models is not None:
      if len(models) != n:
        msg = 'Data track size mismatch: models (%d) must match number of regions (%d)'
        raise Exception(msg % (len(models), n))
    
      self._setData('models', chromoGroup, uint32, models)

    key = (source, code, chromo) 
    if key in self._dataTrackHistogramCache:
      del self._dataTrackHistogramCache[key]      
    
    self._setDataTrackIndicesMb(source, code, [chromo])
    
    
    return chromoGroup
    

  def addDataTrack(self, code, regionDict, valueDict, source=EXTERNAL, annoDict=None, stranded=None):
    """Add a new data layer, with values for superimposing om structral models"""
    
    group = self._getDataTrackGroup(source)
    
    if code in group:
      msg = 'Attempt to add data layer "%s" which already exists.' % code
      raise(Exception(msg))
      
    else:
      self.setDataTrack(code, source, regionDict, valueDict, annoDict, stranded)
    

  def removeDataTrack(self, source, code):
    """Remove a data layer"""
  
    group = self._getDataTrackGroup(source)
  
    if code in group:
      self._delete(code, group)
  
  
  def setVoxelGrid(self, code, data, origin=(0,0,0), gridSize=(1,1,1)):
    """Set values for an existing full volumetric grid data"""

    gridGroup = self.grid
    
    gridData = self._setData(code, gridGroup, float, array(data), 'gzip')
    
    self._setAttr(gridData, 'origin', array(origin, float))
    self._setAttr(gridData, 'gridSize', array(gridSize, float))


  def addVoxelGrid(self, code, data, origin=(0,0,0), gridSize=(1,1,1)):
    """Add a new full volumetric grid data"""
    
    gridGroup = self._getGroup('grid', self.images)
    
    if code in gridGroup:
      raise Exception('Attempt to add existing voxel grid group "%s"' % (code))
    
    self.setVoxelGrid(code, data, origin, gridSize)
    

  def removeVoxelGrid(self, code):
    """Remove a full volumetric grid data"""
    
    if code in self.grid:
      self._delete(code, self.grid)
    
    
  def setVoxelCoords(self, code, data):
    """Set coordinates of existing group of volumetric signal intensities (generally for sparse data)"""

    coordsGroup = self._getGroup('coords', self.images)
    
    self._setData(code, coordsGroup, float, array(data), 'gzip')

    
  def addVoxelCoords(self, code, data):
    """Add new group of volumetric signal intensity coordinates (generally for sparse data)"""

    coordsGroup = self._getGroup('coords', self.images)
    
    if code in coordsGroup:
      raise Exception('Attempt to add existing voxel coords group "%s"' % (code))
    
    self.setVoxelCoords(code, data)
    
    
  def removeVoxelCoords(self, code):
    """Remove a group of volumetric signal intensity coordinates (generally for sparse data)"""

    if 'coords' in self.images:
      if code in self.images['coords']:
        self._delete(code, self.images['coords'])


  def setFoci(self, code, positions, heights=None, volumes=None, widths=None, annotations=None):
    """Set existing lists of spatial signal peaks/foci"""
    
    # height (float) volume (float) x (float) y (float) z (float) dx (float) dy (float) dz (float) Text (char32)

    fociGroup = self._getGroup(code, self.foci)
    numPeaks = len(positions)
   
    if not heights:
      heights = ones(numPeaks, float)
   
    if not volumes:
      volumes = ones(numPeaks, float)
   
    if not widths:
      widths = zeros((numPeaks, 3), float)
    
    positions = array(positions)
    
    dataArray = vstack([positions.T, heights, volumes, widths.T]).T # widths -  3 rows

   
    self._setData('peaks', fociGroup, float, dataArray)
    
    if annotations:
      self._setData('annotations', fociGroup, VL_str, annotations)

    
  def addFoci(self, code, positions, heights=None, volumes=None, widths=None, annotations=None):
    """Add a list of spatial signal peaks"""
     
    if code in self.foci:
      raise Exception('Attempt to add existing foci group "%s"' % (code))
    
    self.setFoci(code, positions, heights, volumes, widths, annotations)
    

  def removeFoci(self, code):
    """Remove a list of spatial signal peaks"""
    
    if code in self.foci:
      self._delete(code, self.foci)
    
   
  def calcContactDataTrack(self, groupName='singleCell', cis=True, trans=True, binSize=int(1e6)): 
    """Create a data track layer from contact information, per chromosome"""
    
    from cUtil import dataLayer
    
    contDict = self.getContacts(groupName, None, cis, trans)
    regionDict = {}
    valueDict = {}
    
    for key in contDict:
      chrA, chrB = key
      
      data = contDict[key]
      
      if chrA == chrB:
        pos = data[:2].ravel()
        regions = array([pos, pos+1], int32).T
        values = ones(len(regions), float)
        
        hist = dataLayer.regionBinValues(regions, values, int32(binSize))
        bins = arange(0.0, binSize*(1+len(hist)), binSize)
        regionDict[chrA] = vstack([bins[:-1], bins[1:]-1]).T
        valueDict[chrA] = vstack([hist, hist/hist.max()]).T
        
      else:
        posA = data[0]
        posB = data[1]
        regionsA = array([posA, posA+1], int32).T
        regionsB = array([posB, posB+1], int32).T
        values = ones(len(regionsA), float)
        
        hist = dataLayer.regionBinValues(regionsA, values, int32(binSize))
        bins = arange(0.0, binSize*(1+len(hist)), binSize)
        regsA = vstack([bins[:-1], bins[1:]-1]).T
        valsA = vstack([hist, hist/hist.max()]).T
         
        if chrA in regionDict:
          regionDict[chrA] = vstack([regionDict[chrA], regsA])
          valueDict[chrA]  = vstack([valueDict[chrA],  valsA]) 
        
        else:
          regionDict[chrA] = regsA
          valueDict[chrA]  = valsA
      
        hist = dataLayer.regionBinValues(regionsB, values, int32(binSize))
        bins = arange(0.0, binSize*(1+len(hist)), binSize)
        regsB = vstack([bins[:-1], bins[1:]-1]).T
        valsB = vstack([hist, hist/hist.max()]).T
        
        if chrB in regionDict:
          regionDict[chrB] = vstack([regionDict[chrB], regsB])
          valueDict[chrB]  = vstack([valueDict[chrB],  valsB])
        
        else:
          regionDict[chrB] = regsB
          valueDict[chrB]  = valsB
    
    for chrA in regionDict:
      regions = regionDict[chrA]
      values  = valueDict[chrA]
      
      valid = values[:,0].nonzero()      
      regionDict[chrA] = regions[valid]
      valueDict[chrA]  = values[valid]
 
     
    if cis == trans:
      name = groupName
      color = (1.0, 1.0, 0.0)
    elif cis:
      name = groupName + '[cis]'
      color = (1.0, 0.0, 1.0)
    else:
      name = groupName + '[trans]'
      color = (0.0, 1.0, 1.0)
        
    return self.setDataTrack(name, DERIVED, regionDict, valueDict, color=color)
    
   
  def calcTransContactGscore(self, groupName='singleCell', binSize=int(10e6)):
    """Calculate the G-test score (c.f. relative entropy) for
       single cell trans contacts."""
    
    # Cold get exp counts from cis/all contacts within bin
    
    obsCounts = {}
    contactDictTrans = self.getContacts(groupName, cis=False)
    chromosomes = self.getChromosomes()
    
    
    for i, chrA in enumerate(chromosomes[:-1]):
      startA, endA = self.getChromosomeLimits(chrA)
      binAe = int(endA/binSize)
      
      for chrB in chromosomes[i+1:]:
        startB, endB = self.getChromosomeLimits(chrB)
        binBe = int(endB/binSize)
        
        for binA in range(binAe+1):
          for binB in range(binBe+1):
          
            key = frozenset([(chrA, binA), (chrB, binB)])
            obsCounts[key] = 0
   
    for pairKey in contactDictTrans:
      chrA, chrB = pairKey
      contacts = contactDictTrans[pairKey]
            
      for posA, posB, count in contacts:
        binA = int(posA/binSize)
        binB = int(posB/binSize)
        
        key = frozenset([(chrA, binA), (chrB, binB)])
        obsCounts[key] += 1
    
    nCont = float(sum(obsCounts.values()))
    pExp = 1.0/float(len(obsCounts)) # isotropic
    
    gScore = 0.0
    for key in obsCounts:
      pObs = obsCounts[key]/nCont
       
      if pObs:  
        gScore += pObs * log(pObs/pExp)
        
    #from scipy.stats import chi2
    #degFree = len(obsCounts)-1
    #pValue = chi2.sf(2*gScore, degFree)
    #print(pValue, gScore, degFree)
 
    return 2*gScore

if __name__ == '__main__':
  
  pass
  
  # TBD import and run testing
  
  
