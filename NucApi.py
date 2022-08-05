import os, sys, shutil, time, colorsys, re
import loadLibs

from datetime import datetime
from inspect import getcallargs, getargspec
from types import FunctionType
from functools import wraps

from h5py import File, Group, special_dtype, h5o
from numpy import array, float32, uint32, int64, int32, int16, uint8, ones, zeros, vstack, arange, mean
from numpy import dstack, abs, sqrt, dot, append, random, empty, log, log2, log10, outer, concatenate, ascontiguousarray
from numpy import power, eye, argsort, argwhere, cumsum, sqrt, hstack, exp, clip, ndarray, string_, unique
from numpy import integer, floating, dtype, cov, diag, median, array_equal, searchsorted, full, repeat, delete

from random import seed, shuffle, randint, uniform
from math import sin, cos, atan, ceil
from string import hexdigits

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
ENGINE = None

PI = 3.14159265358979323846
TAU = 2.0 * PI

# For lots of MPI process pipes
#from resource import getrlimit, setrlimit, RLIMIT_NOFILE
#hardFileLimit = getrlimit(RLIMIT_NOFILE)[1]
#setrlimit(RLIMIT_NOFILE, (hardFileLimit, hardFileLimit))

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
INTERACTIONS = 'interactions'

COLOR_MODES = ['Seq. position','Chromsome ID', 'Data track',
               'Model number', 'RMSD', 'Faint'] # , 'Density']

DISPLAY_MODES = ['Ball and Stick', 'Line', 'Tube', 'Surface']

DATA_TRACK_SYMBOLS = ('None (color region)', 'Circle', 'Star', 'Square', 'Triangle', 'Hexagon')

DATA_TRACK_PEAK_TYPES = ('Histogram', 'Region strip', 'Contact box')

DATA_TRACK_SYMBOL_PARAMS = ((0.0, 0), (TAU/10.0, 10), (TAU/2.5, 5), (TAU/4.0, 4),
                            (TAU/3.0, 3), (TAU/6.0, 6)) # (dAngle, nPoints)

RESTRAINT_COLOR_MODES = ['Cis/Trans', 'Distance', 'Seq. separation', 'Chromosome ID']
  
STRUC_CALC_DEFAULTS = {'numModels':1, 'bboneReg':0, 'bboneSpace':100, 'restrBinned':0,
                       'powerLaw':-0.33, 'seqUnitScale':10, 'restrDist':1.0,
                       'restrErr':0.2, 'tempMax':5000, 'tempMin':10, 'tempSteps':256,
                       'dynSteps':256, 'hierProtocol':1, 'hierStart':4,
                       'hierSteps':1, 'startStruc':0, 'randRad':100.0, 'randSeed':7}

DATA_TRACK_CACHE_BIN = 100000

CROMO_LABEL_DTYPE = dtype('a32,i4,f8,f8,f8,i4,i4,i4')


def isWord(value):

  return len(value.split()) == 1 


def isDataSource(value):

  return value in (DERIVED, EXTERNAL, INNATE)


def isArrayN2(value):

  return array(value).shape == (2,)

  
def argTypeCheck(**argTypeChecks):
  # See http://legacy.python.org/dev/peps/pep-0318/
  
  def makeDecorator(func):
    
    @wraps(func)
    def wrapFunc(*args, **kw):
      
      valDict = getcallargs(func, *args, **kw)
      argSpec = getargspec(func)
      defaultDict = {}
      
      if argSpec.defaults:
        defaultingArgs = argSpec.args[-len(argSpec.defaults):]
        defaultDict = dict(zip(defaultingArgs, argSpec.defaults))
      
      for i, argName in enumerate(argSpec.args):
        if argName in argTypeChecks:
          val = valDict[argName]
          
          if (argName in defaultDict) and (val == defaultDict[argName]):
            continue
          
          checks = argTypeChecks[argName]
          typs = []
          
          if not isinstance(checks, (tuple, list)):
            checks = (checks,)
            
          for check in checks:
            if type(check) is FunctionType:
              assert check(val), 'Value %r failed check "%s"' % (val, check.__name__)
            else:
              typs.append(check)
          
          if typs:
            assert isinstance(val, type(typs)), "Argument %r does not match type(s) %s" % (val,typs)
        
      return func(*args, **kw)
    
    wrapFunc.__name__ = func.__name__

    return wrapFunc
    
  return makeDecorator

    
class NucApiException(Exception):

  pass
  

# Each class to separate file ?


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
        <chrA>/
          <chrB>/ - posA, posB, numObs 
          
    working/
      <groupName>/
        <chrA>/
          <chrB>/ - posA, posB, numObs, model
  
  chromosomes/
    <chrName> - display, color, limits
      labels - label, pos
      regions - [(start, end, display_mode, thickness, selected, r, g, b), 
    
  structures/
    <code>/ - name
      calculation/
 
      particles/
        <chrName>/ - display, color
          backbone  - 0/1
          positions - seqPos
          mapability - 0.0 .. 1.0
 
      restraints/
        <chrA>/
          <chrB>/ idxA, idxB, upper, lower, weight, ..
 
      transforms/
        global - x, y, z
        models - model, x, y, z
 
      coords/
        <chr>/ model, atom, [x, y, z]

  dataTracks/
    derived/ # From contacts/structures
      <code>/ - display, options
        <chr>/ - indicesMb
          regions - first, last
          values  - origVal, normVal
          annotations (optional)  -> "labels"
          models (optional)

    innate/  # From sequence/molecule
      <code>/
        <chr>/ ...
    
    external/
      <code>/
        <chr>/ ...
  
  
  interactions/
    <code>/ - display, options  
      <chrA>/
        <chrB>/ - 
          regions - firstA, lastA, firstB, lastB
          values  - origVal, normVal
          labels (optional)
     
    
  images/
    <code>/ - origin, gridSizes, transform
      grid/ - value [i,j,k]
      coords/  - x, y, z
      intensities/ - [a, (b)]
      shapes/ - [dx,dy,dz]
      labels
     

"""

def backwardCompatibility(root):
  # TBD: Move to separate file
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
            root['dataTracks'][group][code].attrs['display']  = (0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0) # r, g, b, a, scale, minVal, maxVal, maxNorm

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


  if 'interactions' not in root:
     root.create_group('interactions')
     
  for code in root['interactions']:
    if 'details' not in root['interactions'][code].attrs:
      root['interactions'][code].attrs['details'] = code
      
  for group in root['dataTracks']:
    for code in root['dataTracks'][group]:
      if 'details' not in root['dataTracks'][group][code].attrs:
        root['dataTracks'][group][code].attrs['details'] = code
        
      #for chromo in root['dataTracks'][group][code]:
      #  attrs = root['dataTracks'][group][code][chromo].attrs
      #  
      #  if 'indicesMb' not in attrs:
      #    attrs['indicesMb'] = []
        
  
  if 'display' not in root:
    root.create_group('display')
  
  
  if 'models' in root['display']:
    models = root['display']['models']
    del root['display']['models']
  
  else:
    models = [0,]

  if 'original' not in root['contacts']:
    root['contacts'].create_group('original')
    
  if 'working' not in root['contacts']:
    root['contacts'].create_group('working')
  
  if 'voxels' in root:
    del root['voxels']
    
  if 'images' not in root:
    group = root.create_group('images')
  
  if 'foci' in root['images']:
    del root['images']['foci']
    del root['images']['grid']
    del root['images']['coords']
    
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
  
  if 'chromosomes' in root['structures']:  #  Prior to multiple, separate structures
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
 
      if 'colors' in root['structures']['chromosomes'][chromo]:
        del root['structures']['chromosomes'][chromo]['colors']
 
      if 'annotations' in root['structures']['chromosomes'][chromo]:
        del root['structures']['chromosomes'][chromo]['annotations']
    
    # New top-level chromosomes group
    chromoGroup = root.create_group('chromosomes')
    
    # make structure the first of many
    if '0' in root['structures']:
      strucGroup = root['structures']['0']
    else:
      strucGroup = root['structures'].create_group('0')
    
    strucGroup.attrs['name'] = string_(str('1'))  
    strucGroup.attrs['displayModels'] = models
    strucGroup.attrs['isActive'] = 1
     
    for chromo in root['structures']['chromosomes']:
      oldGroup = root['structures']['chromosomes'][chromo]
      newGroup = chromoGroup.create_group(chromo)
      
      newGroup.attrs['display'] = array(oldGroup.attrs['display'])
      newGroup.attrs['color'] = array(oldGroup.attrs['color'])
      newGroup.attrs['limits'] = array([0,0])
      #newGroup.attrs['region_of_interest'] = array([0,0,0])
    
    names = ('calculation', 'chromosomes', 'restraints', 'transforms', 'coords')
    for name in names:
      if (name in root['structures']) and (name not in strucGroup):
        src = root['structures'][name]
        if name == 'chromosomes':
          name2 = 'particles'
        else:
          name2 = name
        
        root.copy(src, strucGroup, name2, shallow=False)    
      
        del root['structures'][name]
    
    for chromo in strucGroup['particles']:
      del strucGroup['particles'][chromo].attrs['display']
      del strucGroup['particles'][chromo].attrs['color']
    
    
  for section in root['contacts']:
    for groupName in root['contacts'][section]:
      attrs = root['contacts'][section][groupName].attrs
      
      if 'diagonalSide' not in attrs:
        attrs['diagonalSide'] = 1
      
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
  
  """            
  
  for group in root['dataTracks']:
    for code in root['dataTracks'][group]:
      chromos = list(root['dataTracks'][group][code].keys())
      
      for chrA in chromos:
        if 'regions' not in root['dataTracks'][group][code][chrA]:
           del root['dataTracks'][group][code][chrA]
      
        elif not len(root['dataTracks'][group][code][chrA]['regions']):
          del root['dataTracks'][group][code][chrA]
        
        elif chrA.lower().startswith('chr'):
          chrB = chrA[3:]
          
          if chrB in root['dataTracks'][group][code]:
            a = array(root['dataTracks'][group][code][chrA]['regions'])
            b = array(root['dataTracks'][group][code][chrB]['regions'])
            
            if array_equal(a, b):
              del root['dataTracks'][group][code][chrA]
  for chromo in root['chromosomes']:
    a, b = root['chromosomes'][chromo].attrs['limits']
    
    if 'region_of_interest' not in root['chromosomes'][chromo].attrs:
      root['chromosomes'][chromo].attrs['region_of_interest'] = (a, b, 0)
    
    elif len(root['chromosomes'][chromo].attrs['region_of_interest']) == 2:
      a, b = root['chromosomes'][chromo].attrs['region_of_interest']
      shown = int((a,b) != tuple(root['chromosomes'][chromo].attrs['limits']))
      root['chromosomes'][chromo].attrs['region_of_interest'] = (a, b, shown)
  """
    
           
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
        repScales.append(0.1 * repScale)
      
      #print repScales
      
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
    
    self._istep = 0
    return self # This is an iterator object

  def __next__(self): # For python 3

    return self.next()
   
  def next(self): # For Python 2
  
    self._istep += 1
    
    if self._istep < self.numAnnealStages:
      return self.annealStages[self._istep]
      
    else:
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
                     ('chromosomes',()),
                     ('dataTracks', (DERIVED, INNATE, EXTERNAL)),
                     ('sample',     ('protocol', 'organism', 'tissue')),
                     ('structures', ('0')),
                     ('images',     ())
                     )
        
        for parent, children in hierarchy:
          group = root.create_group(parent)
        
          for child in children:
            group.create_group(child)
        
        for child in ('particles', 'restraints', 'transforms', 'coords'):
          root['structures']['0'].create_group(child)
        
        now = int(time.time())
        random.seed(now)        
        
        data = array([random.random(), now, now], float32)
        root.attrs['id'] = data
        
        name = fileName if fileName else 'Unknown'
        root['sample'].attrs['name'] = string_(name)  
          
    else:
      raise Exception('%s file name "%s" not valid' % (PROGRAM, fileName))
    
    #try:
    if mode != 'r':
      backwardCompatibility(self.root)
 
    self._setLinkAttrs()
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
    
    #except Exception as err:
    #  self.root.close()
    #  raise err

  
  @argTypeCheck(source=isDataSource)
  def updateProxyDataTracks(self, remoteNuc, source):
    # TBD: private
    
    rNuc = remoteNuc
    rDataGroup = rNuc._getDataTrackGroup(source)
    lDataGroup = self._getDataTrackGroup(source)
    
    # Remove old, unused refs
    for code in list(lDataGroup.keys()):
      if code not in rDataGroup:
        locLayer = self._getGroup(code, lDataGroup)
        chromos = list(locLayer.keys())
        rem = 0
        
        for chromo in chromos:
          if not len(locLayer[chromo]['regions']):
            self._delete(chromo, locLayer)
            rem += 1
        
        if rem == len(chromos):
          self._delete(code, lDataGroup)
    
    # Add new refs
    for code in rDataGroup:
      if (code not in lDataGroup) or not len(lDataGroup[code].keys()):
        refLayer = rDataGroup[code]
        locLayer = self._getGroup(code, lDataGroup)
        
        locLayer.attrs['stranded'] = refLayer.attrs['stranded']
        locLayer.attrs['options']  = refLayer.attrs['options'] 
        locLayer.attrs['display']  = refLayer.attrs['display'] 
        
        if 'details' in refLayer.attrs:
          locLayer.attrs['details']  = refLayer.attrs['details'] 
        else:
          locLayer.attrs['details']  = code
              
        for chromo in refLayer:
          subGroup = self._getGroup(chromo, locLayer)
          
          self._setData('regions', subGroup, uint32, [])
          self._setData('values', subGroup, float32, [])
          break # only one chromo required
              
    return rNuc
    

  def updateProxyInteractions(self, remoteNuc):
    # TBD: private
    
    rNuc = remoteNuc
    
    rIntGroup = rNuc.interactions
    lIntGroup = self.interactions
    
    # Remove old, unused refs
    for code in list(lIntGroup.keys()):
      if code not in rIntGroup:
        locLayer = self._getGroup(code, lIntGroup)
        chromos = list(locLayer.keys())
        found_a = False
        
        for chrA in chromos:
          found_b = False
          
          for chrB in list(locLayer[chrA].keys()):
            if len(locLayer[chrA][chrB]['regions']):
              found_a = True
              found_b = True
            
            else:  
              self._delete(chrB, locLayer[chrA])
        
          if not found_b:
            self._delete(chrA, locLayer)
    
        if not found_a:
          self._delete(code, lIntGroup)
    
    
    # Add new refs
    for code in rIntGroup:
      if (code not in lIntGroup) or not len(lIntGroup[code].keys()):
        refLayer = rIntGroup[code]
        locLayer = self._getGroup(code, lIntGroup)
        
        locLayer.attrs['options']  = refLayer.attrs['options'] 
        locLayer.attrs['display']  = refLayer.attrs['display'] 
        locLayer.attrs['details']  = refLayer.attrs['details'] 
        
        for chrA in refLayer:
          chromoGroup = self._getGroup(chrA, locLayer)
          
          for chrB in chromoGroup:
            subGroup = self._getGroup(chrB, chromoGroup)
            self._setData('regions', subGroup, uint32, [])
            self._setData('values', subGroup, float32, [])
            break # only one chromo pair required
 
    return rNuc
    
      
  def importNucDataTrack(self, remoteNuc, source, code):
    # TBD: Public, type check
  
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
        lLayer.attrs['details']  = rLayer.attrs['details'] 
     
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
      
      #self._setDataTrackIndicesMb(source, code)     
      
        
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
      #rNuc = self.getExperimentRefNuc('a')
      #rNuc.save()
      
      rNuc = self.getExperimentRefNuc()
 
      if rNuc:
        self.updateProxyDataTracks(rNuc, EXTERNAL)
        self.updateProxyInteractions(rNuc)
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
      #rNuc = self.getGenomeRefNuc('a')
      #rNuc.save()
      
      rNuc = self.getGenomeRefNuc()
  
      if rNuc:
        self.updateProxyDataTracks(rNuc, INNATE)
        self.updateProxyInteractions(rNuc)
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
  
  def _checkLoadRefNuc(self, filePath, mode):
  
    try:
      rNuc = Nucleus(filePath, experimentRef=False, genomeRef=False, mode=mode)

    except IOError as err:
      if mode == 'r':
        rNuc = Nucleus(filePath, experimentRef=False, genomeRef=False, mode='a')
        rNuc.save() #By this point any back-compatibility worked
        rNuc.root.close()
        rNuc = Nucleus(filePath, experimentRef=False, genomeRef=False, mode=mode)
        self._warning('Experiment reference file "%s" was updated to the latest version so it loads properly' % filePath)
      else:
        raise err
 
    except Exception as err:
      msg = 'Reference %s file "%s" could not be loaded.\nPython error:\n%s'
      self._warning(msg % (PROGRAM, filePath, err))
      return  
    
    return rNuc
   
   
  def getExperimentRefNuc(self, mode='r'):

    if (mode=='r') and self._experimentRefNuc and \
       (self._experimentRefNuc.fileName == self.experimentRef):
      return self._experimentRefNuc
     
    if self.experimentRef and self._filePathExists(self.experimentRef):
      
      if self._experimentRefNuc: # Mode not 'r' but it was previously
        self._experimentRefNuc.root.close()
        self._experimentRefNuc = None
      
      rNuc = self._checkLoadRefNuc(self.experimentRef, mode)
      
      if mode == 'r': # Only read-only HDF files persist in session
        self._experimentRefNuc = rNuc
      
      return rNuc
  
  
  def getGenomeRefNuc(self, mode='r'):
    
    if (mode=='r') and self._genomeRefNuc and \
       (self._genomeRefNuc.fileName == self.genomeRef):
      return self._genomeRefNuc
    
    if self.genomeRef and self._filePathExists(self.genomeRef):
      if self._genomeRefNuc:
        self._genomeRefNuc.root.close()
        self._genomeRefNuc = None

      rNuc = self._checkLoadRefNuc(self.genomeRef, mode)
      
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
    
    name = str(name)
    
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
      if (parent is self.workContacts) and (name in self.origContacts):
        origGroup = self._getGroup(name, self.origContacts)
        
        for attrib in ('selected','color','binSize','isSingleCell','diagonalSide'):
          group.attrs[attrib] = origGroup.attrs[attrib]
      
      else:
        group.attrs['selected'] = 1
        group.attrs['color'] = (1.0, 1.0, 0.0, 1.0)
        group.attrs['binSize'] = 1
        group.attrs['isSingleCell'] = 1
        group.attrs['diagonalSide'] = 1
      
    return group
  
  
  def _setData(self, name, parent, dataType, dataArray, compression=None):
    
    #print("_setData", name, parent, dataType, dataArray.shape)
    
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
    # TBD: Harmonise better with __init__
  
    root = self.root
    
    self.attrs = root.attrs
    
    # Convenience group attributes
    
    self.contacts     = root['contacts']
    self.chromosomes  = root['chromosomes']
    self.dataTracks   = root['dataTracks']
    self.sample       = root['sample']
    self.structures   = root['structures']
    self.images       = root['images']
    self.display      = root['display']
    self.interactions = root['interactions']
    
    self.origContacts = self.contacts['original']
    self.workContacts = self.contacts['working']
    self.externalData = self.dataTracks[EXTERNAL]
    self.derivedData  = self.dataTracks[DERIVED]
    self.innateData   = self.dataTracks[INNATE]
        
    for chromo in self.chromosomes:
      attrs = self.chromosomes[chromo].attrs
      
      if 'display' not in attrs:
        attrs['display'] = array([1, 1, 1, 1]) # isShown, useLabels, colMode, dispMode

      if 'color' not in attrs:
        attrs['color'] = array([1.0, 0.0, 0.0, 1.0]) # rgba
        
    for typ in self.dataTracks:
      for code in self.dataTracks[typ]:
        attrs = self.dataTracks[typ][code].attrs
        
        if 'options' not in attrs:
          attrs['stranded'] = 0
          attrs['options'] = (0, 0, 1, 0, 0, 0, 0, 0) # shown, labels, symbol, trackType
          attrs['display'] = (0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0) # r, g, b, scale, threshold
    
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
    
    #if 'roi_color' not in displayAttrs:
    #  displayAttrs['roi_color'] = array([1.0, 1.0, 1.0], float)
    #  displayAttrs['roi_options'] = array([0, 0, 0], int) # Color mode, thickness, unused

    if 'chromoAlpha' not in displayAttrs:
      displayAttrs['chromoAlpha'] = 1.0
    
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
    
    if 1: # 'colorSchemes' not in self.display:
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
                                         [1.0,0.5,0.0],
                                         [1.0,1.0,0.0],
                                         [0.0,1.0,0.0],
                                         [0.0,1.0,1.0],
                                         [0.0,0.0,1.0],
                                         [0.5,0.0,1.0],
                                         [1.0,0.0,1.0],
                                         [0.5,0.5,0.5],
                                         [0.5,0.0,0.0],
                                         [0.5,0.25,0.0],
                                         [0.5,0.5,0.0],
                                         [0.0,0.5,0.0],
                                         [0.0,0.5,0.5],
                                         [0.0,0.0,0.5],
                                         [0.25,0.0,0.5],
                                         [0.5,0.0,0.5]], float)
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
      
    for struc in self.structures:
      calcGroup = self._getGroup('calculation', self.structures[struc])
      calcAttrs = calcGroup.attrs
 
      for name in STRUC_CALC_DEFAULTS:
        if name not in calcAttrs:
          calcAttrs[name] = STRUC_CALC_DEFAULTS[name]
    
    active = None
    struc_codes = list(self.structures.keys())
    
    for code in struc_codes:
      strucGroup = self.getStructureGroup(code)
      
      if active:
        strucGroup.attrs['isActive'] = 0
        
      elif strucGroup.attrs['isActive']:
        self.structure = strucGroup
        active = code
    
    if struc_codes:
      if not active:
        strucGroup = self.structures[struc_codes[-1]]
        strucGroup.attrs['isActive'] = 1
        self.structure = strucGroup
    
    else:
      self.structure = self.getStructureGroup(0, '1: Default')
   
  
  def saveAs(self, filePath):
  
    from shutil import copy2
    
    self.root.flush()
    
    if self.fileName == filePath:
      msg = 'Nuc file "Save As" attempted to current file..'
      self._warning(msg)
    
    else:
      copy2(self.fileName, filePath)
      self.__init__(filePath, self.version, self.experimentRef, self.genomeRef)

  
  def getFileName(self):
    # TBD : add property
 
    return self.root.filename
      
  
  def _getParallelEngine(self):
    from parallel.Engine import Engine

    global ENGINE
    if not ENGINE:
      ENGINE = Engine()
  
    return ENGINE
    
  
  def _getDataTrackSource(self, code, source=None):
  
    if source:
      if source.lower() == INNATE:
        dataGroup, typ = self.innateData, INNATE
 
      elif source.lower() == DERIVED:
        dataGroup, typ = self.derivedData, DERIVED
 
      else:
        dataGroup, typ = self.externalData, EXTERNAL

    else:
      for dataGroup, typ in ((self.externalData, EXTERNAL), (self.innateData, INNATE), (self.derivedData, DERIVED)):
        if code in dataGroup:
          break
    
    if code in dataGroup:
      return dataGroup, typ
      
    else:
      return None, None
      
      
  def _getDataTrackGroup(self, source):
  
    if source.lower() == INNATE:
      return self.innateData
    
    elif source.lower() == DERIVED:
      return self.derivedData
    
    else:
      return self.externalData
      
    
  def _checkChromosome(self, chromosome):
    # TBD : Change to a decorator
    
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
    # TBD : Full check required
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
    
    
    return True  
  
  
  def getId(self):
    # TBD : Property
    """Get the unique identifying code for the nucleus file"""
    
    sid = self.root.attrs['id']
    return tuple(sid[:2])
  
  
  def getCreationTime(self):
    # TBD : Property
    """Get the time the Nucleus file was created"""
  
    sid = self.root.attrs['id']
    return datetime.fromtimestamp(sid[1])


  def getAccessTime(self):
    # TBD : Property
    """Get the time the Nucleus file was last acceseed"""
  
    sid = self.root.attrs['id']
    return datetime.fromtimestamp(sid[2])
  
  
  def getSampleInfo(self):
    # TBD : Property
    """Get information relating to the experimental sample"""
    
    sample = self.root['sample']
    
    infoDict = {}
    for key, value in sample.attrs.iteritems():
      infoDict[key] = value
    
    return infoDict
  
  
  def getDisplayedChromosomes(self, structure=None):
    # TBD : Return objects? - New func for names
    """
    Get a list of currently displayed chromsome names,
    optionally filtering on what has coordinates in a particular structure.

    .. describe:: Input
 
    int
 
    .. describe:: Output
    
    str[]
    """
    
    chromos = self.getChromosomes()
    
    if structure is None:
      chromos = [c for c in chromos if self.getChromoDisplayParams(c)[0]]
    
    else:
      avail = set(self._getCoordsGroup(structure).keys())
      chromos = [c for c in chromos if c in avail and self.getChromoDisplayParams(c)[0]]
      
    return chromos
  
  
  def getDisplayedModels(self, structure=None):
    
    if not structure:
      structure = self.structure.name
    
    nModels = self.getNumModels(structure)
    models = [x for x in self.structures[structure].attrs['displayModels'] if x < nModels]
    
    return models
  
  
  def setDisplayedModels(self, models):
    
    nModels = self.getNumModels()
    models = [x for x in models if x < nModels]
    return self._setAttr(self.structure, 'displayModels', models)
      
  
  def getDisplayedDataTracks(self):
    # TBD : Return objects
    
    group = self.dataTracks
    keys = []
    
    for typ in group:
      subGroup = group[typ]
      
      for code in subGroup:
        if 'options' not in subGroup[code].attrs:
          continue
                  
        if subGroup[code].attrs['options'][0]:
          refGroup = self.getRefDataTrackGroup(typ, code)
          
          if refGroup:
            for chromo in refGroup:
              if len(refGroup[chromo]['regions']):
                keys.append( (typ, code) )
                break
    
    return keys
 
 
  def getDisplayedInteractions(self):
  
    group = self.interactions
    codes = set()

    for code in group:
                
      if group[code].attrs['options'][0]:
        refGroup = self.getRefInteractionsGroup(code)
        
        if refGroup: # Check there's some data
          for chrA in refGroup:
            for chrB in refGroup[chrA]:
              if len(refGroup[chrA][chrB]['regions']):
                codes.add(code)
                break
    
    return list(codes)
    
  
  def setInteractomeParams(self, displayMode=None, colorMode=None, renderMode=None,
                           thickness=None, is3d=None, showCis=None, showTrans=None):
  
    vals = list(self.display.attrs['interactome'])
    
    for i, param in enumerate([displayMode, colorMode, renderMode,
                               thickness, is3d, showCis, showTrans]):
      if param is not None:
        vals[i] = int(param)
    
    return self._setAttr(self.display, 'interactome', vals)
  
  
  def setDisplayDetailLevels(self, ballDetail=None, lineSmooth=None,
                             tubeSmooth=None, tubeDetail=None):
    # TBD : Separate properties, e.g. nuc.ballDetail = 1
    
    opts = list(self.display.attrs['detailLevel'])
    
    for i, param in enumerate([ballDetail, lineSmooth, tubeSmooth, tubeDetail]):
      if param is not None:
        opts[i] = max(min(3, int(param)), 0)
    
    return self._setAttr(self.display, 'detailLevel', opts)
    
  
  def setDisplaySizes(self, ball=None, stick=None, line=None,
                      tube=None, density=None):
    # TBD : Separate properties, e.g. ballSize, densityRadius
    # Simplify in HDF5 ?  nuc/display/sizes/ball = 1 ?
  
    opts = list(self.display.attrs['sizes'])
   
    for i, param in enumerate([ball, stick, line, tube, density]):
      if param is not None:
        opts[i] = max(float(param), 0.0)
          
    return self._setAttr(self.display, 'sizes', opts)
    
  
  def setRestraintsDisplayed(self, cis=True, trans=True):
    # TBD : Separate properties
  
    attrs = self.display.attrs
    opts = list(attrs['options'])
    
    showCis = 1 if cis else 0
    showTrans = 1 if trans else 0
    
    opts[2] = showCis
    opts[3] = showTrans
    
    return self._setAttr(self.display, 'options', opts)
    
    
  def setChromoTextDisplayed(self, isShown):
    # TBD : Property
   
    attrs = self.display.attrs
    opts = list(attrs['options'])
    opts[4] = 1 if isShown else 0
    
    return self._setAttr(self.display, 'options', opts)


  def setTextDisplayed(self, isShown):
    # TBD : Property
   
    attrs = self.display.attrs
    opts = list(attrs['options'])
    opts[5] = 1 if isShown else 0
    
    return self._setAttr(self.display, 'options', opts)


  def setScaleDisplayed(self, isShown):
    # TBD : Property
   
    attrs = self.display.attrs
    opts = list(attrs['options'])
    opts[6] = 1 if isShown else 0
    
    return self._setAttr(self.display, 'options', opts)
  
  
  def setRestraintColoring(self, mode):
    # TBD : Property, using enumerateion names
  
    attrs = self.display.attrs
    opts = list(attrs['options'])
    opts[6] = mode
    
    return self._setAttr(self.display, 'options', opts)
      

  def setDataTrackDisplayed(self, typ, code, isShown):
    # TBD: Data track object
   
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']
        val = 1 if isShown else 0
        
        if vals[0] != val:
          vals[0] = val
          self._setAttr(group[code], 'options', vals)
          return True 
  
  def getDataTrackDisplayed(self, typ, code):
  
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']        
        return bool(vals[0]) 
          
  def setDataTrackDetails(self, typ, code, details):
   
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:   
        self._setAttr(group[code], 'details', details)
        return True 
   

  def getDataTrackDetails(self, typ, code):
   
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:     
        return group[code].attrs['details']      
          
          
  def setDataTrackLabelled(self, typ, code, isShown):
    # TBD: Data track object
  
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']
        val = 1 if isShown else 0
        
        if vals[1] != val:
          vals[1] = val
          self._setAttr(group[code], 'options', vals)
          return True 

  def getDataTrackLabelled(self, typ, code):
  
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']        
        return bool(vals[1]) 


  def setDataTrackColor(self, typ, code, rgba):
    # TBD: Data track object
    
    if typ in self.dataTracks:
      rgba = Color(rgba).rgba
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['display']
        
        if list(vals[:4]) != list(rgba):
          vals[:4] = rgba
          self._setAttr(group[code], 'display', vals)
          return True
  
  
  def setDataTrackPeakType(self, typ, code, peakType):
    # TBD: Data track object
    
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
  
  def setInteractionsDetails(self, code, details):
   
    if code in self.interactions:
      group = self.interactions[code]
      self._setAttr(group, 'details', details)
      return True 
   

  def getInteractionsDetails(self, code):

    if code in self.interactions:
      group = self.interactions[code]
      return group.attrs['details']
      
      
  def setInteractionsDisplayed(self, code, isShown):
    # TBD: Data track object
   
    if code in self.interactions:
      group = self.interactions[code]
      vals = group.attrs['options']
      val = 1 if isShown else 0
      
      if vals[0] != val:
        vals[0] = val
        self._setAttr(group, 'options', vals)
        return True 
          
  def getInteractionsDisplayed(self, code):
  
    if code in self.interactions:
      vals = self.interactions[code].attrs['options']
      return bool(vals[0])
  
          
  def setInteractionsLabelled(self, code, isShown):
  
    if code in self.interactions:
      group = self.interactions[code]
      vals = group.attrs['options']
      val = 1 if isShown else 0
      
      if vals[1] != val:
        vals[1] = val
        self._setAttr(group, 'options', vals)
        return True 
 
 
  def getInteractionsLabelled(self, code):
  
    if code in self.interactions:
      vals = self.interactions[code].attrs['options']
      return bool(vals[1])


  def setInteractionsColor(self, code, rgba):
    
    if code in self.interactions:
      rgba = Color(rgba).rgba
      group = self.interactions[code]
      vals = group.attrs['display']
      
      if list(vals[:4]) != list(rgba):
        vals[:4] = rgba
        self._setAttr(group, 'display', vals)
        return True
  
  
  def getInteractionsColor(self, code):
   
    if code in self.interactions:
      vals = self.interactions[code].attrs['display']
      return Color(vals[:4])

  
  def setInteractionsStyle(self, code, lineType=None, lineWidth=None):
    
    mod = False
    
    if code in self.interactions:
      group = self.interactions[code]
       
      if lineType is not None:
        vals = group.attrs['options']
        
        if lineType != vals[2]:
          vals[2] = lineType
          self._setAttr(group, 'options', vals)
          mod = True
 
      if lineWidth is not None:
        vals = group.attrs['display']
        
        if lineWidth != vals[4]:
          vals[4] = lineWidth
          self._setAttr(group, 'display', vals)
          mod = True
    
    return mod

  
  def getInteractionsStyle(self, code):
   
    if code in self.interactions:
      attrs = self.interactions[code].attrs
      vals1 = attrs['options']
      vals2 = attrs['display']
      
      return vals1[2], vals2[4]
 

  def setInteractionsThresholds(self, code, minVal=None, maxVal=None):

    group = self.interactions
    
    if code in group:
      vals = group[code].attrs['display']
      setAttr = False
      
      if (minVal is not None) and (vals[5] != minVal):
        vals[5] = minVal
        setAttr = True

      if (maxVal is not None) and (vals[6] != maxVal):
        vals[6] = max(maxVal, vals[5])
        setAttr = True
      
      if setAttr:
        self._setAttr(group[code], 'display', vals)    
        return True


  def getInteractionsThresholds(self, code):

    group = self.interactions
    
    if code in group:
      attrs = group[code].attrs['display']
    
      return attrs[5], attrs[6]
 
 
  def getDataTrackPeakType(self, typ, code, asInt=True):
    # TBD: Data track object
    
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['options']
        index = min(len(DATA_TRACK_PEAK_TYPES)-1, max(0, vals[3]))
        
        if asInt:
          return index
        else:
          return DATA_TRACK_PEAK_TYPES[index]
          
  
  def getDataTrackColor(self, typ, code):
    # TBD: Data track object
   
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['display']
        return Color(vals[:4])
                
  
  def setDataTrackSymbol(self, typ, code, index):
    # TBD: Data track object
    
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
    # TBD: Data track object
       
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        index = group[code].attrs['options'][2]
        
        if asInt:
          return index
        else:  
          return DATA_TRACK_SYMBOLS[index]
   

  def setDataTrackScale(self, typ, code, val):
    # TBD: Data track object

    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['display']
        
        if vals[4] != val:
          vals[4] = val
          self._setAttr(group[code], 'display', vals)
          return True


  def getDataTrackScale(self, typ, code):
    # TBD: Data track object

    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        return group[code].attrs['display'][4]
        

  def setDataTrackThresholds(self, typ, code, minVal=None, maxVal=None):
    # TBD: Data track object
    
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        vals = group[code].attrs['display']
        setAttr = False
        
        if (minVal is not None) and (vals[5] != minVal):
          vals[5] = minVal
          setAttr = True

        if (maxVal is not None) and (vals[6] != maxVal):
          vals[6] = max(maxVal, vals[5])
          setAttr = True
        
        if setAttr:
          self._setAttr(group[code], 'display', vals)    
          return True


  def getDataTrackThresholds(self, typ, code):
    # TBD: Data track object
    
    if typ in self.dataTracks:
      group = self.dataTracks[typ]
      
      if code in group:
        attrs = group[code].attrs['display']
      
        return attrs[5], attrs[6] or 1.0
           
   
  def areContactsSingleCell(self, groupName):
    # TBD: ContactGroup object
  
    group = self.getContactGroup(groupName)
    
    if group and 'isSingleCell' in group.attrs:
      return group.attrs['isSingleCell'] == 1
      
    return True


  def getContactsBinSize(self, groupName):
   # TBD: ContactGroup object
  
    group = self.getContactGroup(groupName)
    
    if group and 'binSize' in group.attrs:
      return group.attrs['binSize']
    

  def getContactsDiagonalSide(self, groupName):
   # TBD: ContactGroup object
  
    group = self.getContactGroup(groupName)
    
    if group and 'diagonalSide' in group.attrs:
      return group.attrs['diagonalSide']
  
  
  def getDefaultContactGroup(self):
  
    for n in self.workContacts:
      return n
    
    for n in self.origContacts:
      return n
    
    return 'singleCell'
    
    
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
     
  
  def getContactGroupColor(self, groupName):
   # TBD: ContactGroup object
    
    group = self.getContactGroup(groupName)
    
    if group:
      vals = group.attrs['color']
      return Color(vals)
     
  
  def getSelectedContactGroups(self):
  
    names = [name for name, group in self.getContactGroups() if group.attrs['selected']]
    
    if not names:
      groups = self.getContactGroups()
      
      if groups:
        name, group = groups[0]
        names = [name,]
        self._setAttr(group, 'selected', 1)
    
    return names
    
    
  def setContactGroupColor(self, groupName, rgb):
   # TBD: ContactGroup object
  
    group = self.getContactGroup(groupName)
    
    if group:
      color = list(group.attrs['color'])
      rgb = Color(rgb).rgb
      color[:len(rgb)] = rgb
      return self._setAttr(group, 'color', color)    
      
 
  def setSelectedContactGroups(self, groupNames):  
   # TBD: ContactGroup object
    
    groupNames = set(groupNames)
     
    for name, group in self.getContactGroups():
      if name in groupNames:
        self._setAttr(group, 'selected', 1)
      else:
        self._setAttr(group, 'selected', 0)
        
        
  def setChromoGroupSelected(self, groupName, isSelected):
    # TBD: Rename?! ContactGroup object
  
    group = self.getContactGroup(groupName)
    
    if group:
      val = 1 if isSelected else 0
      self._setAttr(group, 'selected', val)
      
    
  def getChromoDisplayParams(self, name, default=(1,1,1,1)):
    # TBD: Chromosome object, separate properties
  
    if name in self.chromosomes:
      if 'display' in self.chromosomes[name].attrs:
        isShown, useLabels, colMode, dispMode = self.chromosomes[name].attrs['display']
        
      else:
        self.chromosomes[name].attrs['display'] = array(default, int)
        isShown, useLabels, colMode, dispMode = default
      
    else:
      isShown, useLabels, colMode, dispMode = default
    
    return isShown, useLabels, colMode, dispMode
  
       
  def resetChromoColors(self, stretchScheme=True):
    
    scheme = self.getColorScheme('chromosome')
     
    if scheme is None:
      scheme = array([[1.0,0.0,0.0],
                      [1.0,0.5,0.0],
                      [1.0,1.0,0.0],
                      [0.0,1.0,0.0],
                      [0.0,1.0,1.0],
                      [0.0,0.0,1.0],
                      [0.5,0.0,1.0],
                      [1.0,0.0,1.0]], float)
    
    chromos = self.getChromosomes()
   
    if stretchScheme:
      colorArray = self.calcInterpolatedColors(scheme, len(chromos), 1.0)
 
      for i, name in enumerate(chromos):
        self._setAttr(self.chromosomes[name], 'color', colorArray[i])
    
    
    else:
      n = len(scheme)
      
      for i, name in enumerate(chromos):
        r, g, b = scheme[i%n][:3]
        self._setAttr(self.chromosomes[name], 'color', [r,g,b,1.0])
   
   
  def setChromoColor(self, chromo, rgba):
    # TBD: Chromosome object
    
    if chromo in self.chromosomes:
      return self._setAttr(self.chromosomes[chromo], 'color', Color(rgba).rgba) 
    

  def getChromoColor(self, name, default=(1.0, 0.0, 0.0, 1.0)):
    # TBD: Chromosome object
    """
    Fetch the current display colour for a named chromosomes, optionally including alpha channel

    .. describe:: Input
 
    str, float[4], type or bool
 
    .. describe:: Output
    
    Color
    """    
    
    if name in self.chromosomes:
      if 'color' in self.chromosomes[name].attrs:
        color = array(self.chromosomes[name].attrs['color'])
        
      else:
        self.resetChromoColors()
        color = array(self.chromosomes[name].attrs['color'])
      
    else:
      color = array(default, float)
    
    return Color(color)
  
  
  def getChromoLabels(self, chromo):
    """
    Get chromosome positional annotations
    """
    
    if chromo in self.chromosomes:
      chromoGroup = self.chromosomes[chromo]
      
      if 'labels' in chromoGroup:
        return array(chromoGroup['labels'], dtype=CROMO_LABEL_DTYPE) 

  
  def setChromoLabels(self, chromo, data):
    """
    data is  a sequence of (text, bpPosition, red, green, blue, optA, optB, optC)
    """
    
    if chromo in self.chromosomes:
      chromoGroup = self.chromosomes[chromo]
      data = array(data, dtype=CROMO_LABEL_DTYPE)
      self._setData('labels', chromoGroup, CROMO_LABEL_DTYPE, data)
  
  
  def addChromoLabel(self, chromo, text, pos, color, opt=None):
  
    if chromo in self.chromosomes:
      chromoGroup = self.chromosomes[chromo]
      dataArray = self.getChromoLabels(chromo)
      
      if not opt:
        opt = (0,0,0)
      
      datum = array([(text, pos,
                    color[0], color[1], color[2],
                    opt[0], opt[1], opt[2]),], dtype=CROMO_LABEL_DTYPE)
      
      if not len(dataArray):
        dataArray = datum
      
      else:
        dataArray = concatenate((dataArray, datum), axis=0)
    
      self._setData('labels', chromoGroup, CROMO_LABEL_DTYPE, dataArray)
    
    
  def setChromoDisplayParams(self, chromo, isShown=None, useLabels=None, colorMode=None, displayMode=None):
    # TBD: Chromosome object, separate properties
 
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
    # TBD: Chromosome object
      
    if chromo in self.chromosomes:
      return bool(self.chromosomes[chromo].attrs['display'][0])
  
        
  def setChromoDisplayed(self, chromo, isShown):
    # TBD: Chromosome object
      
    if chromo in self.chromosomes:
      opts = self.chromosomes[chromo].attrs['display']
      opts[0] = 1 if isShown else 0
    
      return self._setAttr(self.chromosomes[chromo], 'display', opts)
  

  def getChromoRegionsOfInterest(self, chromo):
    # TBD: Chromosome object
      
    if chromo in self.chromosomes:
      chromoGroup = self.chromosomes[chromo]
      
      if 'regions' in chromoGroup:
        rois = array(chromoGroup['regions'], dtype=int)
        if not len(rois):
          return
        
        return rois
        
        
  def setChromoRegionsOfInterest(self, chromo, regions, display_mode=0, colors=None, thickness=None, selected=0):
    # TBD: Chromosome object
      
    if chromo in self.chromosomes:
      
      n_regions = len(regions)
      
      if isinstance(display_mode, int):
        display_mode = full(n_regions, display_mode)
      else:
        display_mode = array(display_mode)
      
      if isinstance(selected, int):
        selected = full(n_regions, selected)
      else:
        selected = array(selected)
              
      if colors:
        colors = array(colors)
        
        if colors.ndim == 1:
          colors = repeat([colors], n_regions, axis=0)
      
      else:
        colors = repeat([0, 128, 255], n_regions, axis=0)
      
      if thickness is None:
        thickness = ones(n_regions)
        
      elif isinstance(thickness, int):
        thickness = full(n_regions, thickness)
      
      data = zeros((n_regions, 8), int)
      data[:,0:2] = regions
      data[:,2]   = display_mode
      data[:,3]   = thickness
      data[:,4]   = selected
      data[:,5:8] = colors
      
      chromoGroup = self.chromosomes[chromo]
      self._setData('regions', chromoGroup, int, data)
  
  
  def removeChromoRegionOfInterest(self, chromo, index):
  
   if chromo in self.chromosomes:
      chromoGroup = self.chromosomes[chromo]
      data = self.getChromoRegionsOfInterest(chromo)
      
      if (data is not None) and (index < len(data)):
        data = delete(data, index, axis=0)
        self._setData('regions', chromoGroup, int, data)
  
  
  def setChromoRegionOfInterest(self, chromo, index, start, end, display_mode=0, color=[0, 128, 255], thickness=1, selected=0):
  
    if chromo in self.chromosomes:
      chromoGroup = self.chromosomes[chromo]
      data = self.getChromoRegionsOfInterest(chromo)
      
      if (data is not None) and (index < len(data)):
        data[index] = array([start, end, display_mode, thickness, selected, color[0], color[1], color[2]], int)
        self._setData('regions', chromoGroup, int, data) 
  
  
  def addChromoRegionOfInterest(self, chromo, start, end, display_mode=0, color=[0, 128, 255], thickness=1, selected=1):
  
    if chromo in self.chromosomes:
      chromoGroup = self.chromosomes[chromo]
      data = self.getChromoRegionsOfInterest(chromo)
      datum = array([[start, end, display_mode, thickness, selected, color[0], color[1], color[2]],], int)
      
      if (data is None) or not len(data):
        data = datum
      
      else:
        data = concatenate((data, datum), axis=0)
    
      self._setData('regions', chromoGroup, int, data)

  
  def setColorScheme(self, name, colors):
  
    schemeGroup = self._getGroup('colorSchemes', self.display)
    schemeAttrs = schemeGroup.attrs  
    schemeAttrs[name] = array(colors, float)
    
      
  def getColorScheme(self, name):    
      
    schemeGroup = self._getGroup('colorSchemes', self.display)
    schemeAttrs = schemeGroup.attrs
    
    if name in schemeAttrs:
      scheme = array(schemeAttrs[name])
      
      if scheme.shape[1] == 3:
        n = len(scheme)
        scheme = concatenate([scheme, ones((n, 1))], axis=1)
      
      return scheme
  
  
  def calcInterpolatedColors(self, schemeColors, nPoints, alpha=None):
        
    from cUtil.apiUtil import getInterpolatedColors
    
    return getInterpolatedColors(schemeColors, nPoints, alpha)
    
        
  def getChromosomes(self, structure=None):
    """Get a list of chromosome or DNA molecule representations"""
    
    if structure is None:
      if self._chromosomes is None:
        names = [c for c in self.chromosomes]
        names.sort(key=naturalKey())
        self._chromosomes = names
 
      else:
        names = self._chromosomes
    
    else:
      particGroup = self._getParticleGroup(structure)
      names = [c for c in particGroup]
      names.sort(key=naturalKey())

    return names
  
  
  def getChromosomeLimits(self, chromosome, structure=None):
    # TBD: Chromosome object
    """
    Get the first and last sequence positions for a chromsome
    Option to consider the particle positions in a given structure
    """
   
    chromoGroup = self.chromosomes 
    start, end = chromoGroup[chromosome].attrs['limits']
    
    if (end == 0) or (structure is not None):
      particleGroup = self._getParticleGroup(structure) # None structure defaults to active/current
      
      if chromosome in particleGroup:
        positions = particleGroup[chromosome]['positions']
        p0 = positions[0]
        p1 = positions[-1]
 
 
        if (p0 < start) or (p1 > end):
          start = min(p0, start)
          end = max(p1, end)
          self.chromosomes[chromosome].attrs['limits'] = array([start, end])
          self.chromosomes[chromosome].attrs['region_of_interest'] = (0, 0, 0)
    
    return start, end
    
    
  def getDomains(self, chromosome):
    # TBD: Redo, returns multiple, checks flag
    """Get a list of compated topological domain regions"""
    
    if 'domains' in self.derivedData and chromosome in self.derivedData['domains']:
      return array(self.derivedData['domains'][chromosome])
    
    
  def setCurrentStructure(self, selected):

    selected = str(selected)
 
    if selected in self.structures:
      for code in self.structures:
        if selected == code:
          isActive = 1
          self.structure = self.structures[selected]
          
        else:
          isActive = 0
       
        self.structures[code].attrs['isActive'] = isActive
  
  
  def getStructureNames(self):
    """
    Get pairs of structure identifying codes and human readable names
    """
    
    codes = []
    names = []
    for code in self.structures:
      if 'name' not in self.structures[code].attrs:
        self.structures[code].attrs['name'] = string_(str(code))
    
      name = self.structures[code].attrs['name']
      codes.append(code)
      names.append(name)
    
    return codes, names
  
  
  def getNextStructureCode(self):
  
    codes = set(self.structures.keys())
    
    i = 0
    while str(i) in codes:
      if not self.getNumModels(i):
        return str(i)
    
      i += 1
    
    return str(i)
  
  
  def _getStructureSubGroup(self, structure, groupName):
    
    
    if structure is None:
      strucGroup = self.structure
    
    else:
      structure = str(structure)
       
      if structure in self.structures:
        strucGroup = self.structures[structure]
 
      else:
        strucGroup = self.getStructureGroup(code=structure, name=None, isActive=None)
      
    if groupName in strucGroup:
      group = strucGroup[groupName]
      
    else:
      group = self._getGroup(groupName, strucGroup)
    
    return group
  
  
  def _getParticleGroup(self, structure=None):
    
    return self._getStructureSubGroup(structure, 'particles')


  def _getCoordsGroup(self, structure=None):
    
    return self._getStructureSubGroup(structure, 'coords')
  
   
  def _getRestraintsGroup(self, structure=None):

    return self._getStructureSubGroup(structure, 'restraints')

  
  def _getCalculationGroup(self, structure=None):
    
    return self._getStructureSubGroup(structure, 'calculation')


  def _getTransformGroup(self, structure=None):
    
    return self._getStructureSubGroup(structure, 'transforms')
  
  
  def getParticlePositions(self, chromosomes=None, cis=True, trans=True,
                           backbone=None, structure=None):
    # TBD: Chromosome property
    """
    Get the chromosomal sequence locations of model particles.
    Defaults to consider the current structure unless otherwise specified
    """
   
    posDict = {}
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    particGroup = self._getParticleGroup(structure)
    restGroup = self._getRestraintsGroup(structure)
    
    for chromo in particGroup:
      if chromo not in chromosomes:
        continue
    
      positions = array(particGroup[chromo]['positions'])
      selection = ones(len(positions), uint32)
      
      if backbone is True:
        bools =  array(particGroup[chromo]['backbone'])
        selection *= bools
      
      elif backbone is False:
        bools =  1 - array(particGroup[chromo]['backbone'])
        selection *= bools
      
      if not cis and (chromo in restGroup):
        indicesA = set(chromo[chromo][:,0])
        indicesB = set(chromo[chromo][:,1])
        indices = set(indicesA) | set(indicesB)
        selection[tuple(indices)] = 0
      
      if not trans:
        for chromoA in restGroup:
          subGroup = restGroup[chromoA]
 
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
  
  def getNumStructures(self):
    """Get the number of structures/bundles"""
  
    return len(self.structures.keys())
  
  
  def getModels(self, structure=None):
    """Get a list of coorinate model indices"""
  
    if not structure:
      structure = self.structure.name
    
    nModels = self.getNumModels(structure)
    models = list(range(nModels))
    
    return models

  
  def getNumModels(self, structure=None):
    # TBD: Structure object too
    """Get number of 3D coordinate structural models"""   
    
    coordsGroup = self._getCoordsGroup(structure)
    
    if coordsGroup is not None:
      for chromo in coordsGroup:
         return coordsGroup[chromo].shape[0]
      
    return 0


  def getNumCoords(self, structure=None):
    # TBD: Structure object too
    """Get number of 3D coordinates in each structural model"""
    
    coordsGroup = self._getCoordsGroup(structure)
       
    if coordsGroup:
      for chromo in coordsGroup:
        return coordsGroup[chromo].shape[1]
      
    return 0
  
  
  def getCachedContacts(self, groupName):
    # TBD private, forget at some point to save memory?
    
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
                     trans=True, asFraction=False, maxSep=None):
    # TBD: ContactGroup object
    """Get the number of chromosomal contacts (without loading them)"""
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    chromosomes = set(chromosomes)
    
    subGroup = self.getCachedContacts(groupName)
    
    if not subGroup:
      return 0

    n = 0
    t = 0
    
    if maxSep:
      for chrA in subGroup:
        for chrB in subGroup[chrA]:
          t += subGroup[chrA][chrB].shape[1]
          
          if chrA != chrB:
            continue
          
          
          contacts = array(subGroup[chrA][chrB], int)
          diffs = abs(contacts[0] - contacts[1])
          comp = diffs <= maxSep
          
          idx = comp.nonzero()
          #print(idx[0])
          c = len(idx[0])
 
          if chrA in chromosomes:
            n += c
    
    else:
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
  
  
  def getNumRestraints(self, chromosomes=None, cis=True, trans=True,
                       asFraction=False, structure=None):
    # TBD: Structure object
    """Get number of distance restraints"""
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    chromosomes = set(chromosomes)
      
    n = 0
    t = 0
    
    restGroup = self._getRestraintsGroup(structure)
    
    for chromoA in restGroup:
      subGroup = restGroup[chromoA]
      
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
                    backboneSep=None, bboneLower=0.8, bboneUpper=1.2, structure=None):
    # TBD: Structure object
    """Get a list of stored distance restraints.
       Backbone sep is distance for a base pair"""
        
    # Filter on number of observations?
    
    restraintDict = {}
    posCache = {}
      
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    chromosomes = set(chromosomes)
    
    # IndexA (uint32) IndexB (uint32) NumObs (uint32)
    # NumObs (float) Weight (uint8) [ Target (float) Lower (float) Upper (float) ]]
    
    restGroup = self._getRestraintsGroup(structure)
    particGroup = self._getParticleGroup(structure)
      
    for chromoA in restGroup:
      if (chromoA not in particGroup) or (chromoA not in chromosomes):
        continue
      
      subGroup = restGroup[chromoA]
      if chromoA in posCache:
        positionsA = posCache[chromoA]
        
      else:
        positionsA = array(particGroup[chromoA]['positions'], uint32)
        
        posCache[chromoA] = positionsA
       
      for chromoB in subGroup:
        if chromoA == chromoB:
          if not cis:
           continue

        elif not trans:
          continue
      
        if (chromoB not in particGroup) or (chromoB not in chromosomes):
          continue
          
        if usePositions:
          if chromoB in posCache:
            positionsB = posCache[chromoB]
            
          else:
            positionsB = array(particGroup[chromoB]['positions'], uint32)
            
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
        if chromo not in particGroup:
          continue
          
        if chromo in posCache:
          positions = posCache[chromo]
 
        else:
          positions = array(particGroup[chromo]['positions'], uint32)

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
  
  
  def getMapability(self, chromosome, structure=None):
    """Get a list of sequence mapability for chromosomal particles"""
    
    particGroup = self._getParticleGroup(structure)
    self._checkChromosome(chromosome)
      
    mapability = array(particGroup[chromosome]['mapability'])
    positions = array(particGroup[chromosome]['positions'])
    
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
  
  
  def normaliseContacts(self, groupName):
    """Normalise population Hi-C data according to coverage/sensitivity.
       Assumes good coverage, i.e. matrix is not particlarly sparse in cis
    """
    
    from cUtil.apiUtil import outerProdRelEntropy
    
    group = self._getContactGroup(groupName, self.origContacts)
    
    if group.attrs['isSingleCell']:
      self._warning('Can only normalise population contact maps')
      return

    cacheDict = self.getCachedContacts(groupName)
    binSize = group.attrs['binSize']
    scaleDict = {}
    
    # get scalings and normalise cis
    for chromo in cacheDict:
      if chromo not in cacheDict[chromo]:
        continue
    
      matrix = array(self.getContactMatrix(chromo, chromo, binSize, groupName), float)
      scale = matrix.sum(axis=0) - diag(matrix)
      scale /= scale.sum() or 1.0
      scaleDict[chromo] = scale     
      
      n, matrix = outerProdRelEntropy(matrix, scale, scale)
      matrix *= n * 10
  
      #print(chromo, matrix.shape)
      
      self.setContactMatrix(groupName, chromo, chromo, matrix, binSize)
    
    #collect valid trans chromo pairs
    transPairs = set() 
    for chromoA in scaleDict:
      for chromoB in scaleDict:
        if chromoA != chromoB:
          transPairs.add(tuple(sorted([chromoA, chromoB])))
    
    #normalse trans using cached values
    for chromoA, chromoB in sorted(transPairs):
      matrix = array(self.getContactMatrix(chromoA, chromoB, binSize, groupName), float)

      n, matrix = outerProdRelEntropy(matrix, scaleDict[chromoA], scaleDict[chromoB])
      matrix *= n * 10
        
      self.setContactMatrix(groupName, chromoA, chromoB, matrix, binSize)
        
   
  def getContacts(self, groupName=None, chromosomes=None,
                  cis=True, trans=True, model=None):
    """Get list of interacting chromosomal positions"""  
    
    if not groupName:
      groupName = self.getDefaultContactGroup()
    
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
   
  
  def getRegionContacts(self, groupName, chromo, startPos, endPos, cis=True, trans=False):
    """
    Get list of chromosomal contacts and their observations
    that are within a specified region
    """  

    contactData = []

    cacheDict = self.getCachedContacts(groupName)
    
    pos = (startPos + endPos)/2.0
    threshold = abs(endPos-startPos)/2.0
          
    if cis and (chromo in cacheDict) and (chromo in cacheDict[chromo]):
      data = array(cacheDict[chromo][chromo])
    
      pos0 = data[0]
      pos1 = data[1]
      numObs = data[2]
      
      deltas0  = abs(pos0 - pos)
      deltas1  = abs(pos1 - pos)
               
      indices0 = argwhere(deltas0 < threshold)
      indices1 = argwhere(deltas1 < threshold)
      
      contactData += [(chromo, pos1[i], numObs[i]) for i in indices0]
      contactData += [(chromo, pos0[i], numObs[i]) for i in indices1]
    
    if trans:
      for chromoA in cacheDict:
        for chromoB in cacheDict[chromoA]:
          if chromoA == chromoB:
            continue
        
          data = array(cacheDict[chromoA][chromoB])
          
          posA = data[0]
          posB = data[1]
          numObs = data[2]
          
          if chromo == chromoA:
            deltas = abs(posA - pos)
            indices = argwhere(deltas < threshold)
            contactData += [(chromoB, posB[i], numObs[i]) for i in indices]
            
          else:
            deltas = abs(posB - pos)
            indices = argwhere(deltas < threshold)
            contactData += [(chromoA, posA[i], numObs[i]) for i in indices]
          
   
    return contactData
  
   
  def getCloseContacts(self, groupName, position, threshold=int(1e6), cis=True, trans=True):
    """Get list of chromosomal contacts and their observations
       that are close to a given position"""  

    print('Function getCloseContacts() is deprecated. Use getRegionContacts() instead')
    
    chromo, pos = position
    start = pos - threshold
    end = pos + threshold
    
    return self.getRegionContacts(groupName, chromo, start, end, cis, trans)
    

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
   
   
  def getContactMatrix(self, chrA, chrB, binSize=int(1e6), groupName=None,
                       regionA=None, regionB=None, model=None):
    """Get full grid of contacts"""  
    
    from cUtil.apiUtil import binContacts
    
    
    if not groupName:
      groupName = self.getDefaultContactGroup()
    
    self._checkChromosome(chrA)
    self._checkChromosome(chrB)
    
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
   
    if regionA:
      startA, endA = regionA
    else:
      # # # # Might need to make this work off the contact limits
      startA, endA = self.getChromosomeLimits(chrA)

    if regionB:
      startB, endB = regionB
    else:
      startB, endB = self.getChromosomeLimits(chrB)
    
    s_a = int(startA/binSize)
    s_b = int(startB/binSize)
    e_a = int(endA/binSize)
    e_b = int(endB/binSize)
    
    extentA = endA - startA
    extentB = endB - startB
    
    n = e_a - s_a + 1
    m = e_b - s_b + 1
    
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


  def getDistanceMatrix(self, chromosome, models=None, binSize=int(1e6),
                        useMinVal=False, structure=None):
    """Get full grid of model distances along the chromosome sequence.
       Averages over all models or selected models."""  
    
    t0 = time.time()
    
    self._checkChromosome(chromosome)
    
    coords = self._getCoordsGroup(structure)
    
    if coords:

      if chromosome not in coords:
        return
      
      chromoCoords = coords[chromosome]
      particGroup = self._getParticleGroup(structure)
      
      start, end = self.getChromosomeLimits(chromosome)
      extent = end - start
      
      n = int(extent//binSize) + 1 
      m = self.getNumModels(structure)
      
      if models:
        models = [i for i in models if i < m]
      else:
        models = range(m)
      
      # Do quicker in Cython
      
      positions = array(particGroup[chromosome]['positions'])
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
        
        
  def getCisDistances(self, position, models=None, step=int(1e6), structure=None):
    """Get list of model distances to a specified point along the chromosome sequence"""
    
    chromo, pos = position
    self._checkChromosome(chromo)
    
    m = self.getNumModels(structure)
    if models:
      models = [i for i in models if i < m]
    else:
      models = range(m)
    
    coords = self._getCoordsGroup(structure)
    
    if coords and (chromo in coords):
      start, end = self.getChromosomeLimits(chromo)
      extent = end - start
      
      n = int(extent//step) + 1
      posB = [(chromo, step * i) for i in range(n)]
      posA = [position] * len(posB)
      
      return self.getPositionDistances(posA, posB, models, structure=structure)
  
  
  def getModelSize(self, model, chromosomes=None, structure=None):
    """Get the size of chromosomal model, in microns, along its longest and shortest orthogonal axes"""
    
    from util.Cluster import principleComponentAnalysis
    
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    coordsGroup = self._getCoordsGroup(structure)
    
    if coordsGroup:
      allCoords = self.getModelCoords(model, chromosomes, structure)
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
      
      
  def getGlobalTransform(self, structure=None):
    """Get affine transformation associated with all models of a structure"""
    
    tformGroup = self._getTransformGroup(structure)
    
    if 'global' in tformGroup:
      return array(tformGroup['global'])
    
    else:
      return eye(4)  
  
  
  def getModelTransform(self, model, structure=None):
    """Get affine transformation associated with structure model"""
    
    tformGroup = self._getTransformGroup(structure)
    
    if 'models' in tformGroup:
      if model < len(tformGroup['models']):
        return array(tformGroup['models'][model])
  
  
  def getAllCoords(self, structure=None):
    """
    Get all 3D coords for a structure, concatenating all chromosomes together.
    Also reurns  dictionary specifying the indices ranges for each chromsome
    """
  
    coordsGroup = self._getCoordsGroup(structure)
    particGroup = self._getParticleGroup(structure)
    
    chromoRanges = {}
    allCoords = None
    
    a = 0
    for chromo in coordsGroup:
      
      coords = array(coordsGroup[chromo])
      n = coords.shape[-2]
      b = a + n
      
      if n:
        if allCoords is None:
          allCoords = coords
        else:
          allCoords = append(allCoords, coords, axis=1)
      
      chromoRanges[chromo] = (a,b)
      a = b
    
    return allCoords, chromoRanges
  
  
  def getChromoSeqPos(self, chromosome, structure=None):
  
    particGroup = self._getParticleGroup(structure)
    
    if chromosome in particGroup:
      
      return array(particGroup[chromosome]['positions'], int32)
  
  
  
  def getChromoCoords(self, chromosome, structure=None, model=None):
    """
    Get the 3D structure model coordinates for a chromosome
    If nodel is None all models are retrieved.
    """
    
    coordsGroup = self._getCoordsGroup(structure)
    particGroup = self._getParticleGroup(structure)
    
    coords = None
    
    if chromosome in coordsGroup:
      coords = array(coordsGroup[chromosome])
      
      if model is not None:
        if not (model < len(coords)):
          msg = 'Model index %d exceeds available models' % model
          raise Exception(msg)
      
        coords = coords[model]
    
    return coords
   
    
  def getModelCoords(self, model, chromosomes=None, structure=None, backbone=None):
    """Retrieve the 3D structure model coordinates of one or more chromosomes.
       Backbone True:only, False:exclude"""

    if chromosomes:
      chromosomes = sorted(chromosomes)
    else:
      chromosomes = self.getChromosomes()
    
    coordsGroup = self._getCoordsGroup(structure)
    particGroup = self._getParticleGroup(structure)

    allCoords = None
    for chromo in coordsGroup:
    
      if chromo not in chromosomes:
        continue
      
      chromoCoords = coordsGroup[chromo]
      nModels = len(chromoCoords)
      coords = array(chromoCoords[min(nModels-1,model)])
      
      if backbone is not None:
        if backbone:
          indices = array(particGroup['backbone']).nonzero()
          coords = coords[indices]
        
        else:
          indices = (particGroup['backbone'] < 0).nonzero()
          coords = coords[indices]
      
      if len(coords):
        if allCoords is None:
          allCoords = coords
        else:
          allCoords = append(allCoords, coords, axis=0)
    
    return allCoords
    
    
  def getModelBoundingBox(self, model, chromosomes=None, structure=None):
    """Get the maximum extent of the structural model along the coordinate axes"""

    allCoords = self.getModelCoords(model, chromosomes, structure)
    
    minPos = allCoords.min(axis=0).tolist()
    maxPos = allCoords.max(axis=0).tolist()
    
    return zip(minPos, maxPos)
  
  
  def getDataTrackCoords(self, source, code, chromosomes=None, structure=None,
                         offsetFrac=0.5, clipMin=False):
    """
    offsetFrac can be 0.0 for start, 1.0 for end, 0.5 for middle of regions etc.
    clip means to only use values above minimum threshold for the data track
    Returns a chromosome keyed dictionary of coordinates
    """
    # Test: getDataTrackCoords('innate', 'gene')
  
    coordsGroup = self._getCoordsGroup(structure)
    dataGroup = self._getDataTrackGroup(source)
    models = self.getModels(structure)

    if not coordsGroup and dataGroup and models:
      return
  
    refLayer = self.getRefDataTrackGroup(source, code)
    if not refLayer:
      return
    
    minVal, maxVal = self.getDataTrackThresholds(source, code)
      
    if not chromosomes:
      chromosomes = list(refLayer.keys())
    
    coordDict = {}
    
    for chromo in chromosomes:
      if chromo not in refLayer:
        continue

      if chromo not in coordsGroup:
        continue
 
      regions = array(refLayer[chromo]['regions'], int32) # start, end
      if not len(regions):
        continue
      
      starts = regions[:,0]
      ends = regions[:,1]
      
      if clipMin:
        values = array(refLayer[chromo]['values'])[:,1]
        idx = (values > minVal).nonzero()

        
        starts = starts[idx]
        ends = ends[idx]

      if not len(starts):
        continue
      
      positions = (1.0-offsetFrac)*starts + offsetFrac*ends
      coords = []
            
      for m in models:
        coordData = self.getPositionCoords(m, positions, chromo, structure, clipEnds=True)
        coords.append(coordData)

      coordDict[chromo] = array(coords)
    
    return coordDict   
  
  
  def getDataTrackDistances(self, source, code, chromo, position, radius=None,
                            offsetFrac=0.5, structure=None, models=None):
    """
    Get a list of (mean) distances from a chromosomal position to a data track, optionally within a given radius 
    Returns a list of (distance, chromosome, start, end, origValue, normValue, annotation). 
    offsetFrac can be 0.0 for start, 1.0 for end, 0.5 for middle of regions etc.
    """
    
    offsetFrac = max(min(1.0, offsetFrac), 0.0)
    
    if not models:
      models = self.getModels(structure)

    queryCoords = array([self.getPositionCoords(m, [position,], chromo)[0] for m in models]) # nModels, 3
    trackCoords = self.getDataTrackCoords(source, code, structure=structure, offsetFrac=offsetFrac) # nModels, nGenes, 3

    distanceList = []
    
    for chromo in trackCoords:    
      modelDists = []
      
      for m in models:
        coordData1 = queryCoords[m]
        coordData2 = trackCoords[chromo][m]
        
        deltas = coordData1 - coordData2
        dists = sqrt((deltas*deltas).sum(axis=1))
        modelDists.append(dists)
      
      meanDists = array(modelDists).mean(axis=0)
        
      regions = array(dataTrack[chromo]['regions'], int32).tolist()
      values = array(dataTrack[chromo]['values']).tolist()
      n = len(regions)
      
      if 'annotations' in dataTrack[chromo]:
        annos = array(dataTrack[chromo]['annotations']).tolist()
      else:
        annos = [None] * n
         
      for i in range(n):
        if radius and (meanDists[i] > radius):
          continue
        
        datum = (meanDists[i], chromo, regions[i][0], regions[i][1], values[i][0], values[i][1], annos[i])      
        distanceList.append(datum)
    
    return distanceList
    
    
  def calcDataTrackPairIteractions(self, code, source=None, label=None):
   
    if not label:
      label = code + '_pairs'
    
    data_track, source = self._getDataTrackSource(code, source)
    data_track = self.getRefDataTrackGroup(source, code)
    
    region_dict = {}
    value_dict = {}
    
    chromos = list(data_track.keys())
    
    for a, chr_a in enumerate(chromos):
      regions_a = array(data_track[chr_a]['regions'], int32)
      na = len(regions_a)
       
      for chr_b in chromos[a:]:
        regions_b = array(data_track[chr_b]['regions'], int32)
        nb = len(regions_b)

        pair_regions = empty((na*nb, 4))
        values = ones((na*nb, 2))
        
    
  def getPositionCoords(self, model, positions, chromo=None, structure=None, clipEnds=False):
    """Get the (interpolated) 3D coordinates for a list of sequential positions"""

    from cUtil.apiUtil import getInterpolatedCoords
     
    coordsGroup = self._getCoordsGroup(structure)

    if coordsGroup:
      particGroup = self._getParticleGroup(structure)
      
      
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
            points = posDict[chromo] = array(particGroup[chromo]['positions'], int32)
 
          else:
            #coordData.append(None)
            continue
          
          p = getInterpolatedCoords(points, coords, array([pos], int32), clipEnds)
          if len(p):
            coordDataAppend(p[0])
        
        coordData = array(coordData)
        
      else:
        coords = array(coordsGroup[chromo][model])
        points = array(particGroup[chromo]['positions'], int32)
        pos = array(positions, int32)
        
        coordData = getInterpolatedCoords(points, coords, pos, clipEnds)
      
      if len(coordData):
        return coordData
      
      
  def getPositionDistances(self, positionsA, positionsB, models=None,
                           chromoA=None, chromoB=None, structure=None):
    """Get the 3D model distances between pairs of sequence positions"""
    
    nModels = self.getNumModels(structure)
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
      coordData1 = self.getPositionCoords(m, positionsA, chromoA, structure)
      coordData2 = self.getPositionCoords(m, positionsB, chromoB, structure)

      if coordData1 is None:
        continue

      if coordData2 is None:
        continue
      
      if coordData1.shape != coordData2.shape:
        continue
     
      deltas = coordData1 - coordData2
      dists = sqrt((deltas*deltas).sum(axis=1))
      modelDists.append(dists)
    
    if modelDists:
      meanDists = array(modelDists).T.mean(axis=1)
    
      return meanDists
  
  
  def getCoordDistances(self, chromoA, chromoB, models=None, structure=None):
    """Get the 3D model distances between particle coords of a given pair of chromosomes"""
    
    coordsGroup = self._getCoordsGroup(structure)
    
    if chromoA not in coordsGroup:
      return
      
    if chromoB not in coordsGroup:
      return  
    
    nModels = self.getNumModels(structure)
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
    
  
  def getRandomCoords(self, model, chromosome, num, structure=None):
    """Get 3D coordinates for random sequence positions along a backbone path"""
    
    self._checkChromosome(chromosome)
    
    particGroup = self._getParticleGroup(structure)
    start = particGroup[chromosome]['positions'][0]
    end = particGroup[chromosome]['positions'][-1]
    
    positions = random.randint(start, end, num)
    
    return self.getPositionCoords(model, positions, chromosome, structure)
  
  
  def getBackboneSpacing(self, structure=None):
    """Get the maximum sequence separation between backbone particles"""
    
    particGroup = self._getParticleGroup(structure)
              
    for chromo in particGroup:
      positions = array(particGroup[chromo]['positions'])
      seps = positions[1:] - positions[:-1]
      return seps.max()
        
  
  def calcNormDataTrack(self, code, source=None, chromosomes=None, method='max', clipValues=None):
    """
    Normalise data track values, preserving any original data.
    Allowed methods are 'untity', 'max', 'orig', 'log', 'quantile'
    If set clipValues is the min and maximum allowed values to which values ouseide this range will be limited   
    """    
             
    dataGroup, source = self._getDataTrackSource(code, source)

    if code in dataGroup:
      dataLayer = dataGroup[code]
      refNuc = self.getRefDataTrackNuc(source, code, 'a')
      refLayer = refNuc._getDataTrackGroup(source)[code]
      
      if not chromosomes:
        chromosomes = self.getChromosomes()

      if not chromosomes:
        chromosomes = dataLayer.keys()
      
      normMax = []
      maxVals = []
      minVals = []
      chromos = []
      stds = []
      meds = []
      
        
       
      for chromo in chromosomes:
        if chromo not in refLayer:
          continue
        
        subGroup = refLayer[chromo]
        valueData  = array(subGroup['values']).T # origValue, normValue
        
        if not len(valueData):
          continue
        
        chromos.append(chromo)
        origValues = valueData[0]
        
        vmin = origValues.min()
        vmax = origValues.max()
        
        if clipValues:
          if vmin < clipValues[0]:
            vmin = clipValues[0]
          
          if vmax > clipValues[1]:
            vmax = clipValues[1]
           
        maxVals.append( vmax ) 
        minVals.append( vmin )
      
      maxVal = max(maxVals)
      minVal = min(minVals)
      deltaVal = maxVal-minVal
      
      if clipValues:
        if clipValues[0] is None:
          clipValues[0] = minVal

        if clipValues[1] is None:
          clipValues[1] = maxVal

      for chromo in chromos:
        subGroup = refLayer[chromo]
        valueData = array(subGroup['values']).T # origValue, normValue
        origValues = valueData[0]
       
        if clipValues:
          minVal, maxVal = clipValues
          origValues = clip(origValues, minVal, maxVal)
                
        if method == 'unity':
          values = ones(len(origValues), float)
          normMax.append(1.0)

        elif method == 'max':
          values = origValues - minVal
          values /= deltaVal
          normMax.append(1.0)

        elif method == 'sqrt':
          values = origValues - minVal
          values /= deltaVal
          values = sqrt(values)
          normMax.append(1.0)
         
        elif method == 'orig':
          values = origValues
          normMax.append(values.max())

        elif method == 'log':
          values = origValues - minVal
          values = log(1.0 + values)
          vMax = log(1.0 + (maxVal-minVal)) or 1.0
          values /= vMax
          normMax.append(1.0)

        elif method == 'inverse':
          values = 1.0 + origValues - minVal
          values = 1.0 / values
          values /= 1.0 / (minVal + 1.0)
          normMax.append(1.0)
         
        elif method == 'quantile':
          order = origValues.argsort()
          step = 1.0/len(origValues)
          refValues = arange(0.0, 1.0+step, step)
          values = refValues[order.argsort()]
          normMax.append(1.0)
  
        valueData[1] = values
        refNuc._setData('values', subGroup, float, valueData.T)
      
      dispVals = list(dataLayer.attrs['display'])
      dispVals[7] = max(normMax or [0.0])
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
        if ('regions' in group[chromo]) and len(group[chromo]['regions']):
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
   

  def getRefInteractionsGroup(self, code, mode='r'):    
    
    intGroup = self.interactions
    
    if code in intGroup:
      group = intGroup[code]
      for chrA in group:
        for chrB in group[chrA]:
          if ('regions' in group[chrA][chrB]) and  len(group[chrA][chrB]['regions']):
            return group
    
    tNuc = self.getExperimentRefNuc(mode)
    if tNuc:
      intGroup = tNuc.interactions
      if code in intGroup:
        return intGroup[code]
    
    gNuc = self.getGenomeRefNuc(mode)
    if gNuc:
      intGroup = gNuc.interactions
      if code in intGroup:
        return intGroup[code]
    
    # Empty, but no ref
    if code in intGroup:
      return intGroup[code]
    
    
  def getRefDataTrackNuc(self, source, code, mode='r'):
  
    dataGroup = self._getDataTrackGroup(source)
    if code in dataGroup:
      group = dataGroup[code]
      for chromo in group:
        if ('regions' in group[chromo]) and len(group[chromo]['regions']):
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
  
  
  def getInteractionsCodes(self):
  
    return sorted(self.interactions.keys())
  
  
  def getDataTrackCodes(self, source):
  
     dataGroup = self._getDataTrackGroup(source)
     
     return sorted(dataGroup.keys())
  
  
  def getInteractions(self, code, chromosomes=None):
  
    intGroup = self.interactions
    
    if code in intGroup:
      outerGroup = intGroup[code]
      refGroup = self.getRefInteractionsGroup(code)
  
      availChromo = set(refGroup.keys())
      for chromo in refGroup:
        availChromo.update(set(refGroup[chromo].keys()))
      
      if not chromosomes:
        chromosomes = availChromo
      
      dataDict = {}
      
      for chromoA in outerGroup:
        if chromoA not in chromosomes:
          continue
      
        groupA = outerGroup[chromoA]
        
        for chromoB in groupA:
          if chromoB not in chromosomes:
            continue
          
          groupB = groupA[chromoB]
          
          chromoPair = (chromoA, chromoB)
          
          regionData = array(groupB['regions']) # start1, end1, start2, end2
          valueData  = array(groupB['values']) # origValue, normValue
          n = len(regionData)          
          
          if 'annotations' in groupB:
            annotations = groupB['annotations']
          
          else:
            annotations = None
          
          dataDict[chromoPair] = (regionData, valueData, annotations)
               
      return dataDict

       
  def getDataTrack(self, source, code, chromosomes=None, model=None):
    """Get a list of superposed data values and annotations"""
  
    dataGroup = self._getDataTrackGroup(source)

    if code in dataGroup:
      dataLayer = dataGroup[code]
      refLayer  = self.getRefDataTrackGroup(source, code)
    
      stranded = dataLayer.attrs['stranded'] # 0, 1
      
      if not chromosomes:
        chromosomes = [chromo for chromo in refLayer]
        #chromosomes = self.getChromosomes()
      
      dataDict = {}
      
      for chromo in chromosomes:
        if chromo not in refLayer:
          continue
        
        regionData = array(refLayer[chromo]['regions'], int32) # start, end
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
    
    dataTrack = self._getDataTrackGroup(source)[code]
    refLayer = self.getRefDataTrackGroup(source, code)
   
    if chromosomes is None:
      chromosomes = [c for c in refLayer]
   
    for chromo in chromosomes:
      if chromo not in refLayer:
        continue
      
      if chromo not in self.chromosomes:
        continue
      
      regions = array(refLayer[chromo]['regions'], int32)   
      
      starts = regions[:,1]
      if not len(starts):
        continue
         
      start, end = self.getChromosomeLimits(chromo)
      n = int(ceil(end/1e6))
      
      starts /= 1000000
      
      indices = zeros(n, int32)
      prev = -1
      for i, j in enumerate(starts):
        if j != prev:
          indices[prev+1:j+1] = i
          prev = j

      if j < n:
        indices[j:] = n
      
      subGroup = self._getGroup(chromo, dataTrack)
      subGroup.attrs['indicesMb'] = indices
           
  
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
      minVal = dispAttrs[5]
      maxVal = dispAttrs[6]
      numBins = float(width)
      binSize = int32((end-start)/numBins)
      maxNorm = self.getDataTrackMaxValue(source, code)
      
      if binSize > 3 * DATA_TRACK_CACHE_BIN:
        key = (source, code, chromo)
 
        if key in self._dataTrackHistogramCache:
          cacheThrsh, cacheValues, cacheRegions = self._dataTrackHistogramCache[key]
 
          if cacheThrsh == minVal:
            values = cacheValues
            regions = cacheRegions
 
          else:
            values, regions = self._setDataTrackHistogramCache(refLayer, source, code,
                                                               chromo, minVal, maxNorm)
 
        else:
          values, regions = self._setDataTrackHistogramCache(refLayer, source, code,
                                                             chromo, minVal, maxNorm)
      
      else:      
        regions = array(refLayer[chromo]['regions'], int32)
        idx = regions[:,0].argsort()
        regions = regions[idx]
        values = array(refLayer[chromo]['values'], float)[idx]
        i, j = searchsorted(regions[:,0],[start, end])
          
        regions = regions[i:j+1]  # start, end
        values = values[i:j+1,1]  #  normValue  


      hist = regionBinValues(regions, values, binSize, int32(start), int32(end), maxNorm, scale, minVal)
         
      colors = array(color + [0, 0, 0, 255] + color, int32) # pos, zero, neg
      pixmap = addPixmapHistogram(pixmap, hist[:width], colors, asDensity)
    
    return pixmap
    
  
  def getChromosomeDataValues(self, code, chromosome, source=None, model=None):
    print('NucApi: Function getChromosomeDataValues() id depricated. Use getDataTrackValues() instead')
    
    return self.getDataTrackValues(code, chromosome, source, model)
    

  def getDataTrackValues(self, code, chromosome, source=None, model=None, minValue=None):
    """Get a list of chromsome data values in stored order, without position information"""
    
    dataGroup, source = self._getDataTrackSource(code, source)

    if dataGroup:
      dataLayer = self.getRefDataTrackGroup(source, code)

      if chromosome not in dataLayer:
        return 
     
      valueData  = array(dataLayer[chromosome]['values'], float) # origValue, normValue
      
      if (model is not None) and ('models' in dataLayer[chromosome]):
        sIndices = (dataLayer[chromosome]['models'] == model).nonzero()
        valueData = valueData[sIndices,:]
 
      if minValue is not None:
        sIndices = (abs(valueData[:,1]) >= minValue).nonzero()
        valueData = valueData[sIndices[0],:]
       
      return valueData
 
 
  def getDataTrackAnnotations(self, code, chromosome, source=None, model=None, minValue=None):
    """For a given chromosomes get a list of superposed regions with data values"""
   
    dataGroup, source = self._getDataTrackSource(code, source)

    if not dataGroup:
      return
      
    dataLayer = self.getRefDataTrackGroup(source, code)
    if chromosome not in dataLayer:
      return  
    
    chromoData = dataLayer[chromosome]    
    
    if 'annotations' not in chromoData:
      return
    
    anno_data = array(chromoData['annotations'])
    
    
    print anno_data.shape, '222'
    
    sIndices = None
      
    if (model is not None) and ('models' in dataLayer[chromosome]):
      idx = (dataLayer[chromosome]['models'] == model).nonzero()
      anno_data = anno_data[idx,:]
    
    if minValue is not None:
      value_data  = array(dataLayer[chromosome]['values']) # origValue, normValue
      
      if idx is not None:
        value_data = value_data[idx,:]
      
      idx2 = (abs(value_data[:,1]) >= minValue).nonzero()
      anno_data = anno_data[idx2[0],:] 
    
    return anno_data
 
 
  def getDataTrackRegions(self, code, chromosome, source=None, model=None, minValue=None):
    """For a given chromosomes get a list of superposed regions with data values"""
   
    dataGroup, source = self._getDataTrackSource(code, source)

    if not dataGroup:
      return
      
    dataLayer = self.getRefDataTrackGroup(source, code)
    if chromosome not in dataLayer:
      return  
    
    chromoData = dataLayer[chromosome]    
    regionData = array(chromoData['regions'], int32) # start, end
    sIndices = None
      
    if (model is not None) and ('models' in dataLayer[chromosome]):
      sIndices = (dataLayer[chromosome]['models'] == model).nonzero()
      regionData = regionData[sIndices,:]
    
    if minValue is not None:
      valueData  = array(dataLayer[chromosome]['values']) # origValue, normValue
      
      if sIndices is not None:
        valueData = valueData[sIndices,:]
      
      sIndices = (abs(valueData[:,1]) >= minValue).nonzero()
      regionData = regionData[sIndices[0],:] 
    
    return regionData
      
  
  def getStructureGroup(self, code=None, name=None, isActive=None):
    # TBD: Private
    """
    Fetch a structure (HDF group) with a given code or make if doesn't exist.
    Resets the structure name if specified.
    """   
    
    if code is None:
      code = self.getNextStructureCode()
    else:
      code = str(code)
    
    strucGroup = self._getGroup(code, self.structures)  
    
    if name:
      if ('name' not in strucGroup.attrs) or (name != strucGroup.attrs['name']):
        strucGroup.attrs['name'] = string_(str(name))
      
    if 'name' not in strucGroup.attrs:
      n = len(self.structures.keys())
      strucGroup.attrs['name'] = string_(str(n))
    
    if isActive is None:
      if 'isActive' not in strucGroup.attrs:
        if len(self.structures.keys()) == 1:
          isActive = 1
        else:
          isActive = 0
        
        strucGroup.attrs['isActive'] =  isActive
    
    elif isActive:
      strucGroup.attrs['isActive'] = 1
     
    else:
      strucGroup.attrs['isActive'] = 0
       
    if 'displayModels' not in strucGroup.attrs:
      strucGroup.attrs['displayModels'] =  array([0,], int)
   
    for subGroup in ('particles', 'restraints', 'transforms', 'coords'):
      self._getGroup(subGroup, strucGroup)
    
    calcGroup = self._getGroup('calculation', strucGroup)
    calcAttrs = calcGroup.attrs
 
    for name in STRUC_CALC_DEFAULTS:
      if name not in calcAttrs:
        calcAttrs[name] = STRUC_CALC_DEFAULTS[name]
   
    return strucGroup
  
  
  def _copyObject(self, parent, source, dest_name):
      
    stack = [(source, parent, dest_name)]
    
    dest = None
    
    while stack:
      s, p, name = stack.pop(0)
      
      if not name:
        name = s.name.split('/')[-1]
      
      if name in p:
        i = 1
        
        if name == '0':
          name = ''
        
        name2 = '%s%d' % (name, i)
        
        while name2 in p:
          i += 1
          name2 = '%s%d' % (name, i)
        
        name = name2
      
      if isinstance(s, Group):
        d = self._getGroup(name, p)
      
      else:
        data = array(s)
        d = self._setData(name, p, s.dtype, data, s.compression)
      
      if not dest:
        dest = d
      
      #print('Copied %s' % (d.name))
        
      for att in s.attrs:
        d.attrs[att] = s.attrs[att]
      
      if isinstance(s, Group):
        for child in s:
          stack.append((s[child], d, None))
    
    return dest
      

  def importCoords(self, filePath, format, structure=None, chromosomes=None):
     
    # Standard formats: {'JSON','TSV','HDF5','NDArray','Nuc3D'}
    from formats.Util import STRUCTURE_FORMATS
        
    if not self._filePathExists(filePath):
      msg = 'File %s does not exist' % filePath
      self._warning(msg)
      return
    
    
    if format == 'Nuc3D':
      nuc = Nucleus(filePath)

      for imp_struc in nuc.structures:
        if not nuc.getNumModels(imp_struc):
          continue
        
        obj = self._copyObject(self.structures, nuc.structures[imp_struc], self.getNextStructureCode())
        code = obj.name.split('/')[-1]
        
        parti_group = nuc._getStructureSubGroup(code, 'particles')
        posDict = {}
        for chromo in parti_group:
          pos = array(parti_group[chromo]['positions'], uint32)
          posDict[chromo] = pos
        
        self.addChromosomes(posDict)
     
      codes = [x for x in self.structures]
      
      if codes:
        self.setCurrentStructure(codes[-1])
          
    else:
   
      if not structure:
        structure = self.getNextStructureCode()
 
      if format.upper() in STRUCTURE_FORMATS:
        format = format.upper()
 
      coordsGroup = self._getCoordsGroup(structure)
      particGroup = self._getParticleGroup(structure)
 
      moduleName = 'formats.%s' % format
      module = __import__(moduleName, globals(), locals(), [moduleName], -1)
      importFunc = getattr(module, 'importCoords')
 
      posDict, coordsDict = importFunc(filePath)

      chromoIds = list(posDict.keys())
      nameDict = self._getChromoNameMapping(chromoIds)
 
      nPos = 0
 
      for key in chromoIds:
        chromo = nameDict[key]
 
        if chromosomes:
          if chromo not in chromosomes:
            continue
 
        pos = posDict[key]
        coords = coordsDict[key]
 
        nPos += len(pos)
        self.addChromosomes({chromo:pos}, structure=structure)
        self.setModelChromosomeCoords(coords, chromo, structure=structure)
    
      return len(posDict), len(coords), nPos
    
  
  def exportCoords(self, filePath, format, structure=None, chromosomes=None, models=None, scale=None, **kw):
    
    # Standard formats: {'PDB','JSON','TSV','XLSX','HDF5','NDArray'}
    
    if structure is None:
      structure = self.structure.name
    
    if not chromosomes:
      chromosomes = self.getDisplayedChromosomes(structure)
    
    moduleName = 'formats.%s' % format
    module = __import__(moduleName, globals(), locals(), [moduleName], -1)
    exportFunc = getattr(module, 'exportCoords')
    
    coordsGroup = self._getCoordsGroup(structure)
    particGroup = self._getParticleGroup(structure)
    
    nModels = self.getNumModels(structure)
    
    if models:
      moddels = tuple([m for m in models if m < nModels])
      nModels = len(models)
    
    nPos = 0
    posDict = {}
    coordsDict = {}
    
    for chromo in chromosomes:
      chromoStr = 'chr%s' % chromo
      positions = array(particGroup[chromo]['positions'])
      
      if not len(positions):
        continue
      
      nPos += len(positions)

      if models:
        models = tuple(models)
        coords = array(coordsGroup[chromo])[models,:]
 
      else:
        coords = array(coordsGroup[chromo])
 
      if scale:
        coords *= scale
      
      posDict[chromoStr] = positions
      coordsDict[chromoStr] = coords
     
    if format == 'PDB': # Requires extra for mappability, remarks etc.
      kw['nuc'] = self
      kw['structure'] = structure
      kw['scale'] = scale or 1.0
    
    exportFunc(filePath, posDict, coordsDict, **kw)
    
    return len(posDict), nModels, nPos
    
  
  def exportPdbFile(self, filePath, chromosomes=None, scale=0.15, extended=False, structure=None):
    """Save any 3D chromosome structure data in a pseudo PDB format"""
    
    msg = '* * * *  Function exportPdbFile() is deprecated. Use exportCoords(format="PDB") instead. * * * *'
    self._warning(msg)
    
    kw = {'extended':extended}
    self.exportCoords(filePath, 'PDB', structure, chromosomes, None, scale, **kw)

  
  def exportContacts(self, filePath, format='PFE', groupName=None, binSize=None,
                     chromosomes=None, cis=True, trans=True, model=None):
    """Export chromosomal contacts as formatted text file"""

    moduleName = 'formats.%s' % format
    module = __import__(moduleName, globals(), locals(), [moduleName], -1)
    
    if binSize:
      exportFunc = getattr(module, 'exportContactMatrix')  
    else:
      exportFunc = getattr(module, 'exportContacts')
    
    if not groupName:
      names = self.getContactGroups()
      
      if not names:
        msg = 'No contacts present to export'
        self._warning(msg)
        return
    
      groupName = names[0][0]
      
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    n = self.getNumContacts(groupName, chromosomes, cis, trans)
    
    if binSize:    
      matrixDict = {}
      startDict = {}
      m = len(chromosomes)
      for i in range(m):
        chrA = chromosomes[i]
        limA = self.getChromosomeLimits(chrA)
      
        for j in range(i,m):
          if (i == j) and not cis:
            continue
          
          if (i != j) and not trans:
            continue  
        
          chrB = chromosomes[j]
          limB = self.getChromosomeLimits(chrB)
          
          matrix = self.getContactMatrix(chrA, chrB, binSize, groupName, model=model)
          chromoPair = (chrA, chrB)
          
          matrixDict[chromoPair] = matrix
          startDict[chromoPair] = (limA, limB)
         
      exportFunc(filePath, matrixDict, startDict, binSize)
    
    else:
      contactDict = self.getContacts(groupName, chromosomes, cis, trans, model)
      exportFunc(filePath, contactDict)
    
    return n
  
  
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
    
    # TBD add more colour options
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
  
  
  def getContactImage(self, groupName, chromosomes=None, binSize=int(1e6), gamma=0.5, model=None):
  
    if not chromosomes:
      chromosomes = self.getChromosomes()
    
    matrices = []
    for chrA in chromosomes:
      row = []
      
      for chrB in chromosomes:
        matrix = array(self.getContactMatrix(chrA, chrB, binSize, groupName, model), float)
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
    
    
  def exportRmsdMatrixImage(self, fileName, format='PNG', models=None,
                            backboneOnly=False, structure=None):
    """Export a pairwise model RMSD matrix as an image file"""
    
    if not models:
      models = range(self.getNumModels(structure))
    
    fileName = self._checkImageFileName(fileName, format)
  
    matrix = self.calcModelRmsdMatrix(models, backboneOnly, structure=structure)
    image = self._getPixmapImage(matrix)
    image.save(fileName, format)
    
    return image
    

  def exportDistanceMatrixImage(self, fileName, chromosome, format='PNG', models=None,
                                binSize=int(1e6), useMinVal=True, structure=None):
    """Export a chromosome distance matrix as an image file"""
 
    if not models:
      models = range(self.getNumModels(structure))
    
    fileName = self._checkImageFileName(fileName, format)
    from itil.Image import pixmapToImage
  
    matrix = self.getDistanceMatrix(chromosome, models, binSize, useMinVal, structure=structure)
    maxVal = matrix.max()
    if maxVal:
      matrix /= maxVal
    
    image = self._getPixmapImage(matrix)
    image.save(fileName, format)

    return image
    
    
  def exportDataTrack(self, source, code, filePath, format='TSV', chromosomes=None):
    """Export a layer of chromosomal data that may be superposed on structures to an external file"""
    
    moduleName = 'formats.%s' % format
    module = __import__(moduleName, globals(), locals(), [moduleName], -1)
    
    exportFunc = getattr(module, 'exportDataTrack')
    
    dataGroup = self._getDataTrackGroup(source)
    dataLayer = dataGroup[code]
    stranded = dataLayer.attrs['stranded']
 
    # TBD: DataTrack.getDataDict(chromosomes)
    dataDict = self.getDataTrack(source, code, chromosomes)
    exportFunc(filePath, dataDict, code, stranded)
    
    numChromos = len(dataDict)
    numRegions = sum([len(dataDict[c][0]) for c in dataDict])
    
    return numChromos, numRegions
    
    
  def exportInteractions(self, code, filePath, format='TSV', chromosomes=None):
    
    moduleName = 'formats.%s' % format
    module = __import__(moduleName, globals(), locals(), [moduleName], -1)
    
    exportFunc = getattr(module, 'exportInteractions')

    dataDict = self.getInteractions(code, chromosomes)
    exportFunc(filePath, dataDict, code)
    
    usedChromos = set()
    for pair in dataDict:
      usedChromos.update(pair)
    
    numChromos = len(usedChromos)
    numRegions = sum([len(dataDict[c][0]) for c in dataDict])
    
    return numChromos, numRegions
    
    
  def _readChainMappingFile(self, filePath):
    
    chain_file_obj = open(filePath)
 
    region_map = {}
 
    for line in chain_file_obj:
      data = line.split()
 
      if not data:
        pos_a = None
        pos_b = None
 
      elif data[0] == 'chain':
        chromo = data[2][3:] # no "chr"
        pos_a = int(data[5])
        pos_b = int(data[10])
 
        if chromo not in region_map:
          region_map[chromo] = []
 
      elif len(data) == 3:
        delta = int(data[0])
        gap_a = int(data[1])
        gap_b = int(data[2])
        region_map[chromo].append((pos_a, pos_a+delta, pos_b, pos_b+delta))
 
        pos_a += gap_a + delta
        pos_b += gap_b + delta
 
      elif len(data) == 1:
        delta = int(data[0])
        region_map[chromo].append((pos_a, pos_a+delta, pos_b, pos_b+delta))

        pos_a = None
        pos_b = None
 
    for chromo in region_map:
      region_map[chromo] = array(sorted(region_map[chromo]), int32)

    return region_map
  
  
  def _convertContactChainMapping(self, contacts, map_regions_a, map_regions_b):
   
    from cUtil.apiUtil import pointRegionsIntersection
    
    contacts = array(contacts, int32)
    
    pos_a = contacts[0]
    pos_b = contacts[1]
    n_obs = contacts[2]

    map_a = map_regions_a[:,:2]
    map_b = map_regions_b[:,:2]
    n = len(pos_a)
 
    mask = zeros(n, int)
 
    idx = pointRegionsIntersection(pos_a, map_a, exclude=False)
    mask[idx] += 1
 
    idx = pointRegionsIntersection(pos_b, map_b, exclude=False)
    mask[idx] += 1
 
    idx = (mask == 2).nonzero() # ChrA and ChrB can be mapped over
 
    pos_a = pos_a[idx]
    pos_b = pos_b[idx]
    n_obs = n_obs[idx]

    new_contacts = zeros((3, len(pos_a)), int32)
    new_contacts[2] = n_obs
 
    idx_a = searchsorted(map_a[:,1], pos_a)
    idx_b = searchsorted(map_b[:,1], pos_b)
 
    for i, j in enumerate(idx_a):
      p = pos_a[i]
      a1, a2, b1, b2 = map_regions_a[j]

      new_contacts[0,i] =  b1 + (p-a1)

    for i, j in enumerate(idx_b):
      p = pos_b[i]
      a1, a2, b1, b2 = map_regions_b[j]

      new_contacts[1,i] = b1 + (p-a1)
    
    return new_contacts
    
    
  def _importSamFile(self, filePath, hdfGroup, binSize=None, cis=True, trans=True, chainMapFile=None):

    # TBD: Needs more checks to make sure it is a properly processed and paired BAM file
    # Rename to readSamFile()
     
    from util.Io import isFileBinary
    from cUtil.samread import readPairedSam
    from cUtil.apiUtil import binContacts, intMatrixToSparse
    import time, gc
    
    t0 = time.time()
    
    if isFileBinary(filePath):
      fileMode = 'rb'
    else:
      fileMode = 'r'

    dataDict, positionDict = readPairedSam(filePath, fileMode)
    
    chromoIds = list(positionDict.keys())
    nameDict = self._getChromoNameMapping(chromoIds)
    
    if chainMapFile:
      region_map = self._readChainMappingFile(chainMapFile)
    
    if binSize and binSize > 1:
      binSize = int32(binSize)
      
      for chrA in dataDict:
        sa, ea = positionDict[chrA]
        n = ceil((ea-sa)/float(binSize))
        if n <= 0:
          continue
 
        subGroup = self._getGroup(nameDict[chrA], hdfGroup)

        for chrB in dataDict[chrA]:
          sb, eb = positionDict[chrB]
          m = ceil((eb-sb)/float(binSize))
          if m <=0:
            continue
          
          contacts = dataDict[chrA][chrB].T
          
          if chainMapFile:
            contacts = self._convertContactChainMapping(contacts, region_map[chrA], region_map[chrB])
           
          binMatrix = zeros((n,m), int32)
          
          bin_start_a = int32(sa/binSize)*binSize
          bin_start_b = int32(sb/binSize)*binSize
          
          binContacts(contacts, binMatrix, bin_start_a, bin_start_b, binSize)
          
          contacts = intMatrixToSparse(binMatrix, binSize, binSize, bin_start_a, bin_start_b)
          self._setData(nameDict[chrB], subGroup, uint32, contacts, compression='gzip')
   
    else:
      for chrA in dataDict:
        subGroup = self._getGroup(nameDict[chrA], hdfGroup)
 
        for chrB in dataDict[chrA]:
          contacts = dataDict[chrA][chrB].T
          
          if chainMapFile:
            contacts = self._convertContactChainMapping(contacts, region_map[chrA], region_map[chrB])
         
          self._setData(nameDict[chrB], subGroup, uint32, contacts, compression='gzip')

    gc.collect()
    
    posDict = {}
    for chromo in positionDict:
      posDict[nameDict[chromo]] = positionDict[chromo]
    
    print('BAM read time: %.4f' % (time.time()-t0))
 
    return posDict
    
  
  def importContactArray(self, dataArray, chromoA, chromoB, posA, posB,
                         groupName=None, isSingleCell=True, updateChromos=True):
    
    if not groupName:
      groupName = self.getDefaultContactGroup()
    
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
        
  def importContactFile(self, *args, **kw):
   
    msg = 'Function Nucleus.importContactFile() depricated, use  Nucleus.importContacts() instead'
    self._warning(msg)
    
    return self.importContacts(*args, **kw)
          
  def _getChromoNameMapping(self, names):
  
    from util.Genome import getChromosomeNames
    
    nameDict = {}
    chromoIds = []
    
    for name in names:
      if len(name) < 4:
        nameDict[name] = name
      elif name.lower() == 'chromosome':
        nameDict[name] = 'chr'
      elif name[:3].lower() == 'chr':
        nameDict[name] = name[3:] or 'chr' # Never blank
      else:
        chromoIds.append(name)
    
    if chromoIds:
      nameDict.update(getChromosomeNames(chromoIds))
        
    return nameDict
    
           
  def importContacts(self, filePath, format='NCC', groupName=None,
                     binSize=None, isSingleCell=True, updateChromos=True, chainMapFile=None):
    """Import chromosomal contacts from a formatted text file or BAM/SAM file.
       Effectively defines the .nuc file, if needed."""
    
    from cUtil.apiUtil import intMatrixToSparse
    from formats.Util import CONTACT_FORMATS
       
    # Removing existing data?
    if not groupName:
      groupName = self.getDefaultContactGroup() 
    
    if not self._filePathExists(filePath):
      msg = 'File %s does not exist' % filePath
      self._warning(msg)
      return
    
    prevChromosomes = self.getChromosomes()
    
    if format == 'Nuc3D':
      nuc = Nucleus(filePath)
      
      posDict = {}
      for chromo in nuc.getChromosomes():
        posDict[chromo] = nuc.getChromosomeLimits(chromo)
      
      for name in nuc.origContacts:
        self._copyObject(self.origContacts, nuc.origContacts[name], name)

      for name in nuc.workContacts:
        self._copyObject(self.workContacts, nuc.workContacts[name], name)
    
    else:
 
      if format.upper() in CONTACT_FORMATS:
        format = format.upper()
 
      posDict = {}
 
      if binSize is not None:
        binSize = int(binSize)
 
      group = self._getContactGroup(groupName, self.origContacts)
      group.attrs['isSingleCell'] = 1 if isSingleCell else 0
      group.attrs['filePath'] = string_(filePath)
        
      if format == 'SAM':
        posDict = self._importSamFile(filePath, group, binSize, chainMapFile)

      else:
        moduleName = 'formats.%s' % format
        module = __import__(moduleName, globals(), locals(), [moduleName], -1)
 
        if isSingleCell or format in ('NCC', 'PFE'): # Sparse
          importFunc = getattr(module, 'importContacts')
 
          # TBD: No exception hiding now, but add try/except here in release version.
          dataDict = importFunc(filePath, binSize)
 
        else: # Matrix
          importFunc = getattr(module, 'importContactMatrix')
          matrixDict, startDict, binSize = importFunc(filePath)
 
          dataDict = {}
          for chromoPair in matrixDict:
            matrix = matrixDict[chromoPair]
            sa, sb = startDict[chromoPair]
 
            contacts = intMatrixToSparse(matrix, binSize, binSize, int32(sa), int32(sb))
            dataDict[chromoPair] = contacts.T
 
        if not dataDict:
          msg = 'No data in file to import'
          self._warning(msg)
          return
 
        chromoPairs = sorted(list(dataDict.keys()))
 
        chromoIds = set()
        for pair in chromoPairs:
          chromoIds.update(pair)
 
        nameDict = self._getChromoNameMapping(chromoIds)
 
        for chromoPair in chromoPairs:
          chromoA, chromoB = chromoPair
          chrA = nameDict[chromoA]
          chrB = nameDict[chromoB]
 
          contacts = array(dataDict[chromoPair]).T
          subGroup = self._getGroup(chrA, group)
          self._setData(chrB, subGroup, uint32, contacts, compression='gzip')
 
          if updateChromos:
            for i, chromo in enumerate(chromoPair):
              chrA = nameDict[chromo]
              pos = contacts[i] # 0 or 1
              pMin, pMax = pos.min(), pos.max()
 
              if chrA in posDict:
                pMin0, pMax0 = posDict[chrA]
                pMin = min(pMin, pMin0)
                pMax = max(pMax, pMax0)
 
              posDict[chrA] = (pMin, pMax)

      if isSingleCell:
        group.attrs['binSize'] = 1
      else:
        group.attrs['binSize'] = binSize or 1

    
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
    
    if not prevChromosomes:
      self.resetChromoColors()
      
    
  def importDataTrack(self, filePath, format, source, code, binSize=1e3, updateChromos=True, feature=None):
    """
    Import a layer of chromosomal data that may be superposed on structures
    from an external file
    """
    
    if not self._filePathExists(filePath):
      return

    if format == 'Nuc3D':
      nuc = Nucleus(filePath)
      
      chrLimits = {}
      for chromo in nuc.getChromosomes():
        chrLimits[chromo] = nuc.getChromosomeLimits(chromo)
      
      for name in nuc.externalData:
        if name not in self.externalData:
          self._copyObject(self.externalData, nuc.externalData[name], name)

      for name in nuc.innateData:
        if name not in self.innateData:
          self._copyObject(self.innateData, nuc.innateData[name], name)
    
      for name in nuc.derivedData:
        if name not in self.derivedData:
          self._copyObject(self.derivedData, nuc.derivedData[name], name)
      
      counts = None
      
    else: 
      moduleName = 'formats.%s' % format
      module = __import__(moduleName, globals(), locals(), [moduleName], -1)
      importFunc = getattr(module, 'importDataTrack')
 
      if binSize is not None:
        binSize = int(binSize)
 
      if format == 'GFF':
        dataDict, name = importFunc(filePath, None, feature)
      else:
        dataDict, name = importFunc(filePath, binSize)
 
      dataGroup = self._getDataTrackGroup(source)
 
      chrLimits = {}
      regionDict = {}
      valueDict = {}
      annoDict = {}
 
      chrIds = list(dataDict.keys())
      nameDict = self._getChromoNameMapping(chrIds)
      stranded = False
 
      max_val = 0
      min_val = None
 
      for chrId in chrIds:
        chromo = nameDict[chrId]
        regionData, valueData, annotations = dataDict[chrId]
        regionDict[chromo] = regionData
        valueDict[chromo] = valueData

        if annotations:
          annoDict[chromo] = annotations
 
        if not stranded and regionData.argmax(axis=1).min() == 0:
          stranded = True
 
        if updateChromos:
          chrLimits[chromo] = regionData.min(), regionData.max()
 
      self.setDataTrack(code, source, regionDict, valueDict, annoDict, stranded)

      counts = (len(dataDict), sum([len(dataDict[c][0]) for c in dataDict]))


    if updateChromos:
      addDict = {}
      chromos = set(self.getChromosomes())
 
      for chrA, limits in chrLimits.iteritems():
        if chrA not in chromos:
          addDict[chrA] = limits
 
      if addDict:
        self.addChromosomes(addDict)
 
    return counts

  
  def importInteractions(self, filePath, format, code, updateChromos=True):
    """
    Import chromosomal interactions that may be superposed on structures
    from an external file
    """
    
    if not self._filePathExists(filePath):
      return
      
    if format == 'Nuc3D':
      nuc = Nucleus(filePath)
      
      chrLimits = {}
      for chromo in nuc.getChromosomes():
        chrLimits[chromo] = nuc.getChromosomeLimits(chromo)
        
      for name in nuc.interactions:
        self._copyObject(self.interactions, nuc.interactions[name], name)
      
      counts = None
      
    else: 
 
      moduleName = 'formats.%s' % format
      module = __import__(moduleName, globals(), locals(), [moduleName], -1)
      importFunc = getattr(module, 'importInteractions')
 
      dataDict, name = importFunc(filePath)
 
      intGroup = self.interactions
 
      chrLimits = {}
      regionDict = {}
      valueDict = {}
      annoDict = {}
 
      chrIds = set()
      for pair in dataDict:
        chrIds.update(pair)
 
      nameDict = self._getChromoNameMapping(chrIds)
 
      for pair in dataDict:
        chrIdA, chrIdB = pair
        chrA = nameDict[chrIdA]
        chrB = nameDict[chrIdB]
        chromoPair = (chrA, chrB)
 
        regionData, valueData, annotations = dataDict[pair]
 
        regionDict[chromoPair] = regionData
        valueDict[chromoPair] = valueData
 
        if annotations:
          annoDict[chromoPair] = annotations
 
        if updateChromos:
          regionsA = regionData[:,:2]
          regionsB = regionData[:,2:]
          chrLimits[chrA] = regionsA.min(), regionsA.max()
          chrLimits[chrB] = regionsB.min(), regionsB.max()
 
      self.setInteractions(code, regionDict, valueDict, annoDict)
 
      numChromos = len(chrIds)
      numRegions = sum([len(dataDict[c][0]) for c in dataDict])
 
      counts = (numChromos, numRegions)

    if updateChromos:
      addDict = {}
      chromos = set(self.getChromosomes())
 
      for chrA, limits in chrLimits.iteritems():
        if chrA not in chromos:
          addDict[chrA] = limits
 
      if addDict:
        self.addChromosomes(addDict)
      
    return counts

  
  def calcStructure(self, contacts, numCpus=1, updateFunc=None,
                   trans=True, bgCalc=False, structure=None):
    # TBD : Structure.calculate               
    """Calculate structure using stored parameters"""
    
    chromosomes = self.getDisplayedChromosomes()
    
    calcGroup = self._getCalculationGroup(structure)
    attrs = calcGroup.attrs
  
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
      domDict = {}
      stages.append( (domDict, (bboneSpace, bboneSpace)) ) # Binned restraints, regular backbone
      
    elif bboneReg:
      domDict = {}
      for chromo in chromosomes:
        start, end = self.getChromosomeLimits(chromo)
        pos = arange(start, end, bboneSpace)
        domDict[chromo] = array([pos[:-1], pos[1:]]).T
      
      stages.append( (domDict, (bboneSpace, bboneSpace)) ) # No binning, regular spacer backbone
      
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
                               updateFunc, bgCalc, numCpus, None, structure)
                        
    return job
    
  
  def _addDuplicateChromosomes(chromosomes, ext='b'):
    
    posDict = {}
    backboneDict = {}
    chromos_b = []
    
    for chromo in chromosomes:
      chromo_b = chromo + ext
      
      if chromo_b in self.chromosomes:
        continue
      
      chromos_b.append(chromo_b)
      chromoGroup = self._getGroup(chromo, self.chromosomes)
      
      posDict[chromo_b] = chromoGroup['positions']
      backboneDict[chromo_b]= chromoGroup['backbone']
    
    self.addChromosomes(positionsDict, backboneDict=None, interpolateCoords=True, structure=None)
  
    for i, chromo in enumerate(chromosomes):
      color = self.getChromoColor(chromo)
      self.setChromoColor(chromos_b[i], color.rgba)

    return chromos_b
    
  
  def annealStructure(self, groupName, chromosomes=None, numModels=1, calcParams=None,
                      callback=None, bgCalc=False, numCpus=1, contactSubsets=None,
                      structure=0, ambigChromos=None):
                      
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
      chromosomes = self.getChromosomes(structure)
    
    if ambigChromos:
      ambigChromos = sorted(set(ambigChromos) & set(chromosomes))
      chromosomes += self._addDuplicateChromosomes(ambigChromos)
    else:
      ambigChromos = []
      
    models = list(range(numModels))
    contactDict  = self.getCachedContacts(groupName)
    
    if groupName in self.origContacts:
      isSingleCell = self.origContacts[groupName].attrs['isSingleCell']
    
    else:
      isSingleCell = self.workContacts[groupName].attrs['isSingleCell']
    
    if structure is None:
      structure = self.getNextStructureCode()
      print('Structure:%s' %  structure)
          
    name = '%s: %d chrs x%d' % (structure, len(chromosomes), numModels)
    strucGroup = self.getStructureGroup(structure, name)
    
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
            
      restrDict, posDict, ambigDict, bboneDict, binLstDict = calcRestraints(chromosomes, ambigChromos, contactDict, isSingleCell,
                                                               annealStage.bboneSep, annealStage.domainDict, 1.0, # scale not used
                                                               calcParams.distPowerLaw, calcParams.distLower, calcParams.distUpper,
                                                               calcParams.minNumObs, calcParams.maxPopDist, contactSubsets)  
      
      if not bgCalc:
        self._cacheRestrDict.append(restrDict)
      
      if i == 0:
        # For first stage setup starting coords for seq positions
        self.addChromosomes(posDict, bboneDict, interpolateCoords=True, structure=structure)

        coords = []
        if calcParams.randStart:
          self.setRandomCoords(models, chromosomes, maxStep=calcParams.randRad/100.0,
                               randWalk=calcParams.randWalk, randSeed=calcParams.randSeed,
                               structure=structure)
 
        for model in models:
          mCoords = self.getModelCoords(model, chromosomes, structure)
 
          if (mCoords is None) or not len(mCoords):
            self.setRandomCoords([model], chromosomes, randWalk=calcParams.randWalk,
                                 randSeed=calcParams.randSeed, structure=structure)
            mCoords = self.getModelCoords(model, chromosomes, structure)
 
          coords.append(mCoords)
 
        coords = array(coords)
      
      # add backbone restraints, get single concatenated restraint arrays      
      indices, dists, ambig = concatenateRestraints(restrDict, posDict, ambigDict, seqScale,
                                                    calcParams.bboneLower, calcParams.bboneUpper)
      
      stageMasses = ones(len(indices),  float) # self.massesFromBins(binLstDict)

      masses.append( stageMasses ) 
      
      #radii.append( self.radiiFromMasses(stageMasses) )
      radii.append( ones(len(indices),  float) )
      
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

    
    radii = concatenate(radii, axis=0)
    restrIndices = concatenate(restrIndices, axis=0)
    restrDists = concatenate(restrDists, axis=0)
    restrAmbig = concatenate(restrAmbig, axis=0)
    
    masses = concatenate(masses, axis=0)
    
    if (numCpus == 1) and not bgCalc:
      coords = self._runAnnealingLocal(structure, coords, posDicts, stageCounts, restrIndices, restrDists, restrAmbig,
                                       temps, repulsScales, timeSteps, dynSteps, callback, bgCalc, numCpus,
                                       masses=masses, radii=radii)
      job = None
      
    else:  
      # run the annealing via wrapper that handles temp schedule, repulsion and parallelisation
      job = self._runAnnealing(structure, coords, posDicts, stageCounts, restrIndices, restrDists, restrAmbig,
                               temps, repulsScales, timeSteps, dynSteps, callback, bgCalc, numCpus,
                               masses=masses, radii=radii)
 
      # Need to know which coords to set when background job completes
      job.structure = structure
      
    # store final restraints
    restraintGroup = self._getRestraintsGroup(structure)
    used = set()
     
    for chrA in restrDict:
      rGroup = self._getGroup(chrA, restraintGroup)
    
      for chrB in restrDict[chrA]:
        dataArray = restrDict[chrA][chrB]
        self._setData(chrB, rGroup, float32, dataArray)
        used.add((chrA, chrB))
    
    # remove old, unused groups
    for chrA in list(restraintGroup.keys()):
       for chrB in list(restraintGroup[chrA].keys()):
         if (chrA, chrB) not in used:
           self._delete(chrB, restraintGroup[chrA])
       
       
       if not [c for c in restraintGroup[chrA]]:
         self._delete(chrA, rGroup)
    
    # store seq positions for coords 
    self.addChromosomes(posDict, bboneDict, interpolateCoords=False, structure=structure)

    if not bgCalc:
      # store final coords
      
      if job:
       coords = job.getResult()
      
      self.setAllCoords(coords, chromosomes, structure)
 
      # superimpose ensemble
      if len(coords) > 1:
        self.modelAlign(chromosomes=chromosomes, structure=structure)
 
      if callback:
        callback()
      
      self.save()
        
    return job


  def nldrStructure(self, groupName, chromosomes=None, numModels=1, numCpus=1, contactSubsets=None,
                    structure=0, calcParams=None, ambigChromos=None):
                      
    """Calculate chromosome structures using distance restraints from
       contact data using simple annealing protocol"""
    
    from cUtil.apiUtil import calcRestraints, concatenateRestraints

    # Future: consider masses and fixed particles    
     
    # Some initialisation 
    
    if not calcParams:
      calcParams = StrucCalcParams() # i.e. use defaults
      calcParams.addAnnealStage()    # One default, binned stage  
    
    if chromosomes:
      chromosomes = sorted(chromosomes)
    else:
      chromosomes = self.getChromosomes(structure)
      
    models = list(range(numModels))
    contactDict  = self.getCachedContacts(groupName)
    
    if groupName in self.origContacts:
      isSingleCell = self.origContacts[groupName].attrs['isSingleCell']
    
    else:
      isSingleCell = self.workContacts[groupName].attrs['isSingleCell']
    
    if structure is None:
      structure = self.getNextStructureCode()
      print('Structure:%s' %  structure)

    if ambigChromos:
      ambigChromos = sorted(set(ambigChromos) & set(chromosomes))
      chromosomes += self._addDuplicateChromosomes(ambigChromos)
    else:
      ambigChromos = [] 
               
    name = '%s: %d chrs x%d' % (structure, len(chromosomes), numModels)
    strucGroup = self.getStructureGroup(structure, name)
    
    # Loop through annealing stages, collate restraints and calc parameters
    
    restrIndices = []  # Particle indices of restrained pairs
    restrDists = []    # Distrance restraint limits
    restrAmbig = []    # Restraint ambiguity data
    posDicts = []      # Particle seq positions for each chromosome
    stageCounts = []   # Num restraints in each stage
    
    annealStage = calcParams.annealStages[0]
    
    
    if annealStage.bboneSep is not None: # TBD: Think about what to do here in the long-run
      seqScale = annealStage.bboneSep[1]
    else:
      seqScale = calcParams.seqScale   
            
    restrDict, posDict, ambigDict, bboneDict, binLstDict = calcRestraints(chromosomes, ambigChromos, contactDict, isSingleCell,
                                                                          annealStage.bboneSep, annealStage.domainDict, 1.0, # scale not used
                                                                          calcParams.distPowerLaw, calcParams.distLower, calcParams.distUpper,
                                                                          calcParams.minNumObs, calcParams.maxPopDist, contactSubsets)

    self.addChromosomes(posDict, bboneDict, interpolateCoords=True, structure=structure)

    coords = []
    if calcParams.randStart:
      self.setRandomCoords(models, chromosomes, maxStep=calcParams.randRad/100.0,
                           randWalk=calcParams.randWalk, randSeed=calcParams.randSeed,
                           structure=structure)
 
    for model in models:
      mCoords = self.getModelCoords(model, chromosomes, structure)
 
      if (mCoords is None) or not len(mCoords):
        self.setRandomCoords([model], chromosomes, randWalk=calcParams.randWalk,
                             randSeed=calcParams.randSeed, structure=structure)
        mCoords = self.getModelCoords(model, chromosomes, structure)
 
      coords.append(mCoords) 
    coords = array(coords)    
    
    from sklearn.neighbors import NearestNeighbors
    from sklearn.manifold import Isomap, MDS
   
    n_neighbors=10
   
    # add backbone restraints, get single concatenated restraint arrays      
    indices, dists, ambig = concatenateRestraints(restrDict, posDict, ambigDict, seqScale,
                                                  calcParams.bboneLower, calcParams.bboneUpper)
    
    dists = dists.mean(axis=1)
    ensemble_new = []
    
    for model in models:
      n_coords = len(coords[model])
      dist_mat = full((n_coords, n_coords), 100.0)
      
      for k, (i, j) in enumerate(indices):
        dist_mat[i,j] = dists[k]
        dist_mat[j,i] = dists[k]
        
      """
      nn = NearestNeighbors(n_neighbors=n_neighbors, radius=2.0, algorithm='auto',
                            metric='precomputed', n_jobs=numCpus)
      nn.fit(dist_mat)
      im = Isomap(n_neighbors=n_neighbors, n_components=3,
                  eigen_solver='auto', tol=0,
                  max_iter=1000, path_method='D',
                  neighbors_algorithm='auto')
      im.fit(nn) # shape (n_samples, n_components)
      """
      #mds = MDS(metric=True, n_components=3,
      #          verbose=1,
      #          max_iter=300, n_jobs=numCpus,
      #          dissimilarity='precomputed')

      im = Isomap(n_neighbors=n_neighbors, n_components=3,
                  eigen_solver='auto', tol=0,
                  max_iter=1000, path_method='D',
                  neighbors_algorithm='auto')
      
      coords_new = im.fit_transform(dist_mat)

      ensemble_new.append(coords_new)
   
    # store final restraints
    restraintGroup = self._getRestraintsGroup(structure)
    used = set()
     
    for chrA in restrDict:
      rGroup = self._getGroup(chrA, restraintGroup)
    
      for chrB in restrDict[chrA]:
        dataArray = restrDict[chrA][chrB]
        self._setData(chrB, rGroup, float32, dataArray)
        used.add((chrA, chrB))
    
    # remove old, unused groups
    for chrA in list(restraintGroup.keys()):
       for chrB in list(restraintGroup[chrA].keys()):
         if (chrA, chrB) not in used:
           self._delete(chrB, restraintGroup[chrA])
       
       
       if not [c for c in restraintGroup[chrA]]:
         self._delete(chrA, rGroup)
    
    # store seq positions for coords 
    self.addChromosomes(posDict, bboneDict, interpolateCoords=False, structure=structure)

    self.setAllCoords(array(ensemble_new), chromosomes, structure)
 
    # superimpose ensemble
    if len(ensemble_new) > 1:
      self.modelAlign(chromosomes=chromosomes, structure=structure)
 
    self.save()     
      

  def setupStrucCalcJobFile(self, file_path, groupName, chromosomes=None, numModels=1, calcParams=None,
                           callback=None, bgCalc=False, contactSubsets=None, structure=0, ambigChromos=None):
                      
    """Calculate chromosome structures using distance restraints from
       contact data using simple annealing protocol"""
    
    from cUtil.apiUtil import calcRestraints, concatenateRestraints
    
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
      chromosomes = self.getChromosomes(structure)
    
    if ambigChromos:
      ambigChromos = sorted(set(ambigChromos) & set(chromosomes))
      chromsomes += self._addDuplicateChromosomes(ambigChromos)
    else:
      ambigChromos = []
      
    models = list(range(numModels))
    contactDict  = self.getCachedContacts(groupName)
    
    if groupName in self.origContacts:
      isSingleCell = self.origContacts[groupName].attrs['isSingleCell']
    
    else:
      isSingleCell = self.workContacts[groupName].attrs['isSingleCell']
    
    if structure is None:
      structure = self.getNextStructureCode()
      print('Structure:%s' %  structure)
          
    name = '%s: %d chrs x%d' % (structure, len(chromosomes), numModels)
    strucGroup = self.getStructureGroup(structure, name)
    
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
            
      restrDict, posDict, ambigDict, bboneDict, binLstDict = calcRestraints(chromosomes, ambigChromos, contactDict, isSingleCell,
                                                               annealStage.bboneSep, annealStage.domainDict, 1.0, # scale not used
                                                               calcParams.distPowerLaw, calcParams.distLower, calcParams.distUpper,
                                                               calcParams.minNumObs, calcParams.maxPopDist, contactSubsets)
      
      if not bgCalc:
        self._cacheRestrDict.append(restrDict)
      
      if i == 0:
        # For first stage setup starting coords for seq positions
        self.addChromosomes(posDict, bboneDict, interpolateCoords=True, structure=structure)

        coords = []
        if calcParams.randStart:
          self.setRandomCoords(models, chromosomes, maxStep=calcParams.randRad/100.0,
                               randWalk=calcParams.randWalk, randSeed=calcParams.randSeed,
                               structure=structure)
 
        for model in models:
          mCoords = self.getModelCoords(model, chromosomes, structure)
 
          if (mCoords is None) or not len(mCoords):
            self.setRandomCoords([model], chromosomes, randWalk=calcParams.randWalk,
                                 randSeed=calcParams.randSeed, structure=structure)
            mCoords = self.getModelCoords(model, chromosomes, structure)
 
          coords.append(mCoords)
 
        coords = array(coords)
      
      # add backbone restraints, get single concatenated restraint arrays      
      indices, dists, ambig = concatenateRestraints(restrDict, posDict, ambigDict, seqScale,
                                                    calcParams.bboneLower, calcParams.bboneUpper)
      
      stageMasses = ones(len(indices),  float) # self.massesFromBins(binLstDict)

      masses.append( stageMasses ) 
      
      #radii.append( self.radiiFromMasses(stageMasses) )
      radii.append( ones(len(indices),  float) )
      
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
    
    radii = concatenate(radii, axis=0)
    masses = concatenate(masses, axis=0)
    
    stageCounts = array(stageCounts)
    temps = array(temps)
    repulsScales = array(repulsScales)
    timeSteps = array(timeSteps)
    dynSteps = array(dynSteps)
 
    from numpy import savez

    kwArgs = {'coords':coords,
              'stageCounts':stageCounts, 'restrIndices':restrIndices,
              'restrDists':restrDists, 'restrAmbig':restrAmbig,
              'temps':temps, 'repulsScales':repulsScales,
              'timeSteps':timeSteps, 'dynSteps':dynSteps}
              #'masses':masses , 'radii':radii}
 
    for i, pos_dict in enumerate(posDicts):
      for chromo in pos_dict:
        key = 'pos_%d_%s' % (i, chromo)
        kwArgs[key] = pos_dict[chromo]
 
    # store restraints used
    restraintGroup = self._getRestraintsGroup(structure)
    used = set()
     
    for chrA in restrDict:
      rGroup = self._getGroup(chrA, restraintGroup)
    
      for chrB in restrDict[chrA]:
        dataArray = restrDict[chrA][chrB]
        self._setData(chrB, rGroup, float32, dataArray)
        used.add((chrA, chrB))
    
    # remove old, unused groups
    for chrA in list(restraintGroup.keys()):
       for chrB in list(restraintGroup[chrA].keys()):
         if (chrA, chrB) not in used:
           self._delete(chrB, restraintGroup[chrA])
       
       if not [c for c in restraintGroup[chrA]]:
         self._delete(chrA, rGroup)
    
    # store seq positions for coords 
    self.addChromosomes(posDict, bboneDict, interpolateCoords=False, structure=structure)
    
    savez(file_path, **kwArgs)
     
    self.save()  # Save restraints used etc.
        

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
      radii.append(mass**(1/3.0))
      # radii.append(2.25)

    radii = array(radii, float)
    #print(radii)

    return radii
        
      
  def _strucUpdateCallback(self, callback, structure, newCoords, posDict, chromosomes, stage, step):

    # update restraints
    if step == 0 and (stage < len(self._cacheRestrDict)):
      restrDict = self._cacheRestrDict[stage]
      restraintGroup = self._getRestraintsGroup(structure)
      for chrA in restrDict:
        rGroup = self._getGroup(chrA, restraintGroup)
 
        for chrB in restrDict[chrA]:
          dataArray = restrDict[chrA][chrB]
          self._setData(chrB, rGroup, float32, dataArray)
    
    # update chromsome positions
    if step == 0:
      self.addChromosomes(posDict, interpolateCoords=False, structure=structure)
    
    # update coordinates
    self.setAllCoords(newCoords, chromosomes, structure)
    callback()
    
  
  def _runAnnealing(self, structure, coords, posDicts, stageCounts, indices, dists, ambig, temps,
                    repulsions, timeSteps, dynSteps, callback=None, bgCalc=False,
                    numCpus=1, minDist=1.5, repDist=1.5, printInterval=100,
                    masses=None, radii=None):
    """
    Wrapper to calulate structures via simulated annealing using concatenated restraint arrays.
    Controls the temperature and repulsive schedule. Handles parallelisation.
    Handles the display callback.
    """
    
    from solve.SimAnnealJob import simAnnealjob
      
    nModels = len(coords)
       
    engine = self._getParallelEngine()
    
    sameArgs = {'posDicts':posDicts, 'stageCounts':stageCounts,
                'indices':indices, 'dists':dists, 'ambig':ambig,
                'temps':temps, 'repulsions':repulsions, 'timeSteps':timeSteps,
                'dynSteps':dynSteps, 'minDist':minDist,
                'repDist':repDist, 'masses':masses, 'radii':radii}
    
    printIntervals = zeros(nModels, int)
    printIntervals[0] = printInterval  # First model prints progress     
    
    # An input arg to specify func for graphical update
    _callback = lambda v, w, x, y, z: self._strucUpdateCallback(callback, structure, v, w, x, y, z)
    callbackArg = ('callback', _callback) if callback and not bgCalc else None
        
    diffArgs = {'coords':coords,
                'printInterval':printIntervals}
    
    job = engine.run(simAnnealjob, sameArgs, diffArgs, numCpus,
                     wait=not bgCalc, combineFunc=array,
                     callbackArg=callbackArg)

    return job 
  
  
  def _runAnnealingLocal(self, structure, coords, posDicts, stageCounts, indices, dists, ambig, temps,
                         repulsions, timeSteps, dynSteps, callback=None, bgCalc=False,
                         minDist=1.5, repDist=1.5, printInterval=100,
                         masses=None, radii=None):
    
    from solve.SimAnnealJob import simAnnealjob

    new_coords = []
    
    for model_coords in coords: 
      model_coords = simAnnealjob(model_coords, posDicts, stageCounts, indices, dists, ambig, temps, repulsions,
                                  timeSteps, dynSteps, minDist, repDist, printInterval, callback=None, masses=masses, radii=radii)
    
      new_coords.append(model_coords)
 
    new_coords = array(new_coords)

    return new_coords  
    
    
  def calcBootstrapRmsds(self, groupName, chromosomes=None, rmsdWeightScale=10.0):
    """
    Calculate structure model RMSDS within and between bootstrap partition bundles and overall
    """
    #TBD: Confidence interval on coorindate mean .e.g width of 95% tail limits
    
    from util.Structure import superimposeCoordArray, superimposeCoordPair
    
    contactGroup = self.origContacts[groupName]
    numPartitions, numResamples, strucCodes = self._getCrossValidParams(contactGroup)

    # Cache coords
    coordDict = {}
    for strucCode in strucCodes:
      if not chromosomes:
        chromosomes = self.structures[strucCode]['coords'].keys()
    
      nModels = self.getNumModels(strucCode)
      coords = []
      
      for i in range(nModels):
        coords.append(self.getModelCoords(i, chromosomes, strucCode))
      
      coordDict[strucCode] = array(coords)

    rmsdsResample = []
    rmsdsBundle = []
    rmsdsPartition = []
    for a, strucCodeA in enumerate(strucCodes[:-1]):
      coordsA = coordDict[strucCodeA]
      s1 = a // numPartitions
      p1 = a % numPartitions
      
      for b, strucCodeB in enumerate(strucCodes[a:], a):
        coordsB = coordDict[strucCodeB]
        s2 = b // numPartitions
        p2 = b % numPartitions
        
        rmsds = []
        for i in range(len(coordsA)):
          for j in range(len(coordsB)):
            result = superimposeCoordPair(coordsA[i], coordsB[j], rmsdWeightScale)
            coords2, rot, rmsd, atomRmsds = result
            rmsds.append(rmsd)
        
        if a == b:  # Within same bundle
          rmsdsBundle.extend(rmsds)    
            
        elif s1 == s2:  
          if p1 != p2:# Between partitions of the same sampling
            rmsdsPartition.extend(rmsds) 
          
        else:
          rmsdsResample.extend(rmsds) # Compare across resamplings
    
    return array(rmsdsBundle), array(rmsdsPartition), array(rmsdsResample)
  
  
  def setBootstrapOverviewStructure(self, groupName, chromosomes=None):
  
    from cUtil.apiUtil import interpolateChromoModelCoords
    contactGroup = self.origContacts[groupName]
    numPartitions, numResamples, strucCodes = self._getCrossValidParams(contactGroup)
  
    # Cache coords, get unified seq positions
    coordDict = {}
    prevPosDict = {}
    posDict = {}
    for strucCode in strucCodes:
      if not chromosomes:
        chromosomes = self.structures[strucCode]['coords'].keys()
    
      prevPosDict[strucCode] = {}
      nModels = self.getNumModels(strucCode)
      particGroup = self._getParticleGroup(strucCode)
      
      for chromo in chromosomes:
        chromoPartGroup = self._getGroup(chromo, particGroup)
        pos = array(chromoPartGroup['positions'], int32)
        prevPosDict[strucCode][chromo] = pos
        
        coords = []
        for i in range(nModels):
          coords.append((strucCode, self.getModelCoords(i, [chromo,], strucCode)))
              
        if chromo in posDict:
          posDict[chromo].update(pos)
          coordDict[chromo].extend(coords)
        else:
          posDict[chromo] = set(pos)
          coordDict[chromo] = coords
    
    for chromo in posDict:
      pos = posDict[chromo]
      posDict[chromo] = array(sorted(pos), int32)
    
    # Make combined structure  
    strucCode = self.getNextStructureCode()
    strucGroup = self.getStructureGroup(strucCode, 'BootstrapCombined')
    particGroup = self._getParticleGroup(strucCode)
    
    for chromo in posDict:
      chromoPartGroup = self._getGroup(chromo, particGroup)
      positions = posDict[chromo]
      n = len(positions)
      
      backbone = zeros(n, uint32)
      mapability = ones(n, float32)
       
      self._setData('backbone', chromoPartGroup, uint32, backbone)
      self._setData('positions', chromoPartGroup, float32, positions)
      self._setData('mapability', chromoPartGroup, float32, mapability)
      
      pDictA = {chromo: positions}
      chromoCoords = []
      
      for m, (strucCodePrev, coords) in enumerate(coordDict[chromo]):
        pDictB = {chromo: prevPosDict[strucCodePrev][chromo]}
        coords = interpolateChromoModelCoords(pDictA, pDictB, coords)
        chromoCoords.append(coords)
      
      chromoCoords = array(chromoCoords)
      
      self.setModelChromosomeCoords(chromoCoords, chromo, None, strucCode)

  
  def _getCrossValidParams(self, contactGroup):
  
    if 'cvPartitions' not in contactGroup.attrs:
      err = Exception('Contact group "%s" does not contain bootstrap calculation information')
      raise(err % contactGroup.name)
    
    numPartitions = contactGroup.attrs['cvPartitions']
    numResamples = contactGroup.attrs['cvResamples']
    strucCodes =  contactGroup.attrs['cvStructures']
    strucCodes = [str(c) for c in strucCodes]
    
    for code in strucCodes:
      if code not in self.structures:
        err = Exception('Bootstrap calculation structure %s missing')
        raise(err % code)
    
    return numPartitions, numResamples, strucCodes
      
      
  def calcBootstrapCrossValid(self, groupName):
    # Look at distance distributions of contacts left-out from bootstrap calculation
    # e.g. random k(10)-fold cross validation

    contactDict = self.getCachedContacts(groupName)
    contactGroup = self.origContacts[groupName]
    numPartitions, numResamples, strucCodes = self._getCrossValidParams(contactGroup)
    
    distsExclude = []
    distsInclude = []
 
    for i in range(numResamples):
      
      for j in range(numPartitions):
        k = (i * numPartitions) + j
        structCode = strucCodes[k]
        coordsGroup = self._getCoordsGroup(structCode)
        numModels = self.getNumModels(structCode)
        
        for chrA in contactDict:
          if chrA not in coordsGroup:
            continue
        
          subGroup = contactGroup[chrA]
        
          for chrB in contactDict[chrA]:
            if chrB not in coordsGroup:
              continue
              
            if 'cvRandSeeds' not in subGroup[chrB].attrs:
              continue # This pair not used in calc 
          
            seed( subGroup[chrB].attrs['cvRandSeeds'][i] ) # Seed used to get random indices at calc time
            contacts = contactDict[chrA][chrB]
            n = len(contacts[0])
            partitions = empty(n, int32)
            indices = list(range(n))
            shuffle(indices)
          
            for a in range(n):
              partitions[indices[a]] = a % numPartitions
          
            contacts = contactDict[chrA][chrB]
            indices = (partitions == j).nonzero() # Excluded indices
            
            posA = contacts[0, indices].flatten() # Excluded contacts
            posB = contacts[1, indices].flatten()
            
            modelDists = []
            for m in range(numModels):
              coordData1 = self.getPositionCoords(m, posA, chrA, structCode)
              coordData2 = self.getPositionCoords(m, posB, chrB, structCode)
 
              deltas = coordData1 - coordData2
              mDists = sqrt((deltas*deltas).sum(axis=1))
              modelDists.append(mDists)
 
            meanDists = array(modelDists).T.mean(axis=1)
            
            distsExclude.append(meanDists)
            
            indices = (partitions != j).nonzero() # Included indices
            
            posA = contacts[0, indices].flatten() # Included contacts
            posB = contacts[1, indices].flatten()
            
            modelDists = []
            for m in range(numModels):
              coordData1 = self.getPositionCoords(m, posA, chrA, structCode)
              coordData2 = self.getPositionCoords(m, posB, chrB, structCode)
 
              deltas = coordData1 - coordData2
              mDists = sqrt((deltas*deltas).sum(axis=1))
              modelDists.append(mDists)
 
            meanDists = array(modelDists).T.mean(axis=1)
            
            distsInclude.append(meanDists)
    
    return array(distsExclude).ravel(), array(distsInclude).ravel()
            
  
  def annealStructureBootstrap(self, groupName, chromosomes=None, 
                               numModels=10, numPartitions=10, numResamples=10,
                               calcParams=None, numCpus=None):
     
    if not calcParams:
      calcParams = StrucCalcParams() # i.e. use defaults
      calcParams.addAnnealStage()    # One default, binned stage    
     
    if chromosomes:
      chromosomes = sorted(chromosomes)
    else:
      chromosomes = self.getChromosomes()
      
      
    # Clear all previous structures
    strucCodes = [c for c in self.structures]
    for code in strucCodes:
      self._delete(code, self.root)    
    
    # Make new structure slots
    strucCodes = []
    strucGroups = []
    for i in range(numResamples):
      for j in range(numPartitions):
        k = (i * numPartitions) + j
        code = self.getNextStructureCode()
        name = '%d: s%02d p%02d %d chrs x%d' % (k+1, i+1, j+1, len(chromosomes), numModels)
        strucGroup = self.getStructureGroup(code, name)
        strucCodes.append(code)

    # Normal BG annealing for each partition+resample, separate jobs
    # Must have same final fixed bins with spacers (not native resolution due to differing contacts)
    
    contactGroup = self.origContacts[groupName]
    randSeeds = {} # So that the partitioning can be regnerated without storing the indices
    for chrA in contactGroup:
      randSeeds[chrA] = {}
      for chrB in contactGroup[chrA]:
        randSeeds[chrA][chrB] = []
   
    structureStack = list(strucCodes)
    for i in range(numResamples):
      if groupName in self._contactsCache:
        del self._contactsCache[groupName]
      contactDict = self.getCachedContacts(groupName) # Clean copy
      
      # Mark partitions in contacts
      for chrA in contactDict:
        subGroup = contactGroup[chrA]
      
        for chrB in contactDict[chrA]:
          tSeed = int(time.time())
          seed(tSeed)
          randSeeds[chrA][chrB].append(tSeed)
          
          contacts = contactDict[chrA][chrB]
          n = len(contacts[0])
          partitions = empty((1,n), int32)
          indices = list(range(n))
          shuffle(indices)
          
          for a in range(n):
            partitions[0,indices[a]] = a % numPartitions
          
          contacts = concatenate([contacts, partitions], axis=0).astype(uint32)
          
          # Store in cache, not permanently           
          self._contactsCache[groupName][chrA][chrB] = contacts
          
      
      for j in range(numPartitions):
        strucCode = structureStack.pop(0)
        parts = [k for k in range(numPartitions) if k != j]
        
        job = self.annealStructure(groupName, chromosomes, numModels, calcParams,
                                   callback=None, bgCalc=False, numCpus=numCpus,
                                   contactSubsets=parts, structure=strucCode)
    
    # Store parameters on contact group for later cross-validation
    contactGroup.attrs['cvPartitions'] = numPartitions
    contactGroup.attrs['cvResamples']  = numResamples
    contactGroup.attrs['cvStructures'] = array([int(x) for x in strucCodes])
    
    for chrA in contactGroup:
      for chrB in contactGroup[chrA]:
        seeds = randSeeds[chrA][chrB]
 
        if seeds:
           contactGroup[chrA][chrB].attrs['cvRandSeeds'] = array(seeds)
                              
    return job
  
  
  def calcContactVoidRegions(self, groupName=None, trackName='void', close_cis=int(1e6), source=DERIVED):
  
    from cUtil import dataLayer
    from numpy import histogram
    
    #from matplotlib import pyplot
    #from matplotlib.colors import LinearSegmentedColormap
    
    if not groupName:
      groupName = self.getDefaultContactGroup()    
      
    contDict = self.getContacts(groupName, None, cis=True, trans=False)
    binSize = self.getContactsBinSize(groupName)
    
    #cmap = LinearSegmentedColormap.from_list(name='YB', colors =['#000000', '#0080FF', '#FFFF00'], N=51)
 
    regionDict = {}
    valueDict = {}
    
    for key in contDict:
      chrA, chrB = key
      
      start, end = self.getChromosomeLimits(chrA)

      num_bins = int((end-start)/binSize)
      
      if not num_bins:
        continue
    
      obs = self.getContactMatrix(chrA, chrA, binSize, groupName).astype(float)
      obs -= diag(diag(obs))
      obs_sum = obs.sum()
      
      if not obs_sum:
        continue
      
      obs /= obs_sum
      n = len(obs)
      
      data = contDict[key]
      
      hist, edges = histogram(data[:2].ravel(), bins=num_bins, range=(start, end+binSize))            
      
      if len(hist) < 2:
        continue
      
      selection = zeros(len(hist), int)
      
      med = median(hist)

      idx = (hist < 0.5 * med).nonzero()
      selection[idx] = 1
      
      idx = (hist > 2.0 * med).nonzero()
      selection[idx] = 1
      
      seq_seps = abs(data[0]-data[1])
      
      idx = (seq_seps < close_cis).nonzero()[0]
      close_cis_data = data[:,idx]
      
      values = close_cis_data[:2].ravel() # All contact ends, i.e. flatten to diagonal
      c_hist, edges = histogram(values, bins=num_bins, range=(start, end+binSize))
      
      med = median(c_hist)
      idx = (c_hist < 0.5 * med).nonzero()
      selection[idx] = 1
 
      counts = zeros(n)
      sig = zeros(n)
 
      for i in range(n):
        for j in range(i,n):
          d = j-i
          sig[d] += obs[i,j]
          counts[d] += 1.0
 
      for c, j in enumerate(counts):
        if c:
          sig[j] /= c
 
      sig /= sig.sum()
 
      exp = zeros((n,n), float)
      for i in range(n):
        exp[i,:i+1] = sig[:i+1][::-1]
        exp[i,i:] = sig[:n-i]
 
      vals = obs.sum(axis=0)
      vals /= vals.sum()
 
      exp *= outer(vals, vals)
      exp /= exp.sum()
 
      idx = (exp * obs).nonzero()
      z = ((exp * obs) == 0.0).nonzero()
 
      h = array(obs)
      h[idx] /= exp[idx]
      
      g = h.sum(axis=0)
      g -= median(g)
      std = 1.4826 * median(abs(g))
      g /= std or 1.0
      
      #print(chrA, len(g), len(hist))
      
      #selection2 = zeros(len(h), int)
      idx2 = (abs(g) > 4.0).nonzero()
      idx2 = clip(0, len(selection)-1, idx2[0])
      selection[idx2] = 1
      
      idx = selection.nonzero()
      
      """
      
      h = clip(h, 0.0, 2.0*h.std())
    
      fig, axarr = pyplot.subplots(2, 2) 
      fig.suptitle('Hi-C clean-up & void zones Chr %s' % chrA)
 
      axarr[0,0].matshow(log10(obs*1e8+1.0), cmap=cmap)
      axarr[0,0].set_title('Original data')
      
      axarr[0,1].matshow(h, cmap=cmap)
      axarr[0,1].set_title('Abnormal num contacts at seq sep')
      
      
      h2 = zeros(h.shape)
      
      h2[idx] = 0.5
      h2[:,idx] = 0.5
      
      h2[idx2] = 1.0
      h2[:,idx2] = 1.0
      
      axarr[1,0].matshow(h2, cmap=cmap)
      axarr[1,0].set_title('Void zones')
      
      obs2 = array(obs)
      obs2[idx] = 0.0
      obs2[:,idx] = 0.0
      
      obs2[idx2] = 0.0
      obs2[:,idx2] = 0.0
      axarr[1,1].matshow(log10(obs2*1e8+1.0), cmap=cmap)
      axarr[1,1].set_title('Cleaned data')
      """
      
      # TBD merge adjasent regions to save space
      
      starts = edges[:-1]
      ends = edges[1:]-1
      
      starts = starts[idx]
      ends = ends[idx]
      
      if len(starts) < 1:
        continue
      
      valueDict[chrA] = ones(len(ends), float)
      
      regionDict[chrA] = array(zip(starts, ends), int32)
    
      #pyplot.show()
    
        
    return self.setDataTrack(trackName, source, regionDict, valueDict, color=(1.0, 0.0, 1.0))
 
  
  def calcRegionViolDataTrack(self, group_name=None, track_name=None, structure=None,
                              win_size=1, threshold=4.0, dist_power_law=-0.33):
    
    from cUtil.apiUtil import calcRestraints
    
    if not group_name:
      group_name = self.getDefaultContactGroup()
 
    if not track_name:
      track_name = 'region_viols_' + group_name
    
    particle_size = self.getBackboneSpacing(structure)
    num_models = self.getNumModels()
    models = range(num_models)
    
    partic_group = self._getParticleGroup(structure)
    contact_dict = self.getCachedContacts(group_name)
    chromosomes = contact_dict.keys()
    
    chromo_counts = {}
    chromo_values = {}
    
    restrDict, posDict, ambigDict, bboneDict, binLstDict = calcRestraints(chromosomes, [], contact_dict, True,
                                                                          [particle_size, particle_size], {}, 1.0, # scale not used
                                                                          dist_power_law, 0.8, 1.2, 1)
    
    for chr_a in contact_dict:
      if chr_a not in posDict:
        continue
      
      chr_pos_a = array(partic_group[chr_a]['positions'])
      res_pos_a = posDict[chr_a]
      
      s_a = (chr_pos_a[0] - res_pos_a[0])/particle_size   # Num leading extra particles in restraints 
      coords_a = self.getChromoCoords(chr_a, structure)
      n_a = len(chr_pos_a)
      
      if coords_a is None:
        continue
      
      if chr_a not in chromo_counts:
        chromo_counts[chr_a] = zeros(n_a, float)
        chromo_values[chr_a] = zeros(n_a, float)
      
      
      for chr_b in contact_dict[chr_a]:
        if chr_b not in posDict:
          continue
          
        chr_pos_b = array(partic_group[chr_b]['positions'])
        res_pos_b = posDict[chr_b]
        
        s_b = (chr_pos_b[0] - res_pos_b[0])/particle_size   # Num leading extra particles in restraints 
        coords_b = self.getChromoCoords(chr_b, structure)
        n_b = len(chr_pos_b)
    
        if coords_b is None:
          continue
 
        if chr_b not in chromo_counts:
          chromo_counts[chr_b] = zeros(n_b, float)
          chromo_values[chr_b] = zeros(n_b, float)
 

        restraints = restrDict[chr_a][chr_b] # [i, j, weight, target, lower, upper] 
        
        idx_a = array(restraints[0], int) - s_a
        idx_b = array(restraints[1], int) - s_b
        target = restraints[3]
        
        idx = (idx_a >= 0).nonzero()
        idx_a = idx_a[idx]
        idx_b = idx_b[idx]
        target = target[idx]
        
        idx = (idx_b >= 0).nonzero()
        idx_a = idx_a[idx]
        idx_b = idx_b[idx]
        target = target[idx]
        
        idx = (idx_a < n_a).nonzero()
        idx_a = idx_a[idx]
        idx_b = idx_b[idx]
        target = target[idx]
        
        idx = (idx_b < n_b).nonzero()
        idx_a = idx_a[idx]
        idx_b = idx_b[idx]
        target = target[idx]
        
        if not len(idx_a):
          continue
        
        bead_coords_a = coords_a[:,idx_a]
        bead_coords_b = coords_b[:,idx_b]
        
        dists = []

        for m in models:
 
          deltas = bead_coords_a[m] - bead_coords_b[m]
          model_dists = sqrt((deltas*deltas).sum(axis=1))
          restraint_deltas = model_dists - target
          
          region_means = []
 
          for i in range(win_size):
            j = win_size - i
            region_means.append(restraint_deltas[i:-j])

          region_means = mean(region_means, axis=0)
 
          dists.append(region_means)
 
        dists = mean(dists, axis=0)
        
        i = win_size/2
        j = win_size-i
        idx_a = idx_a[i:-j]
        idx_b = idx_b[i:-j]
        
        idx = (dists > threshold).nonzero()[0]
        
        if len(idx):
          dists = dists[idx]
          idx_a = idx_a[idx]
          idx_b = idx_b[idx]
          
          chromo_values[chr_a][idx_a] = array([dists, chromo_values[chr_a][idx_a]]).max(axis=0)
          chromo_values[chr_b][idx_b] = array([dists, chromo_values[chr_b][idx_b]]).max(axis=0)
    
    region_dict = {}
    values_dict = {}
    
    from matplotlib import pyplot as plt
    
    for chromo in chromo_counts:
      chr_pos = array(partic_group[chromo]['positions'])

      idx = chromo_values[chromo].nonzero()[0]
      
      if not len(idx):
        continue
      
      values = chromo_values[chromo][idx]
      starts = chr_pos[idx] + 1
      ends   = chr_pos[idx] + particle_size
      
      region_dict[chromo] = array([starts, ends], int32).T
      values_dict[chromo] = array([values, values/values.max()]).T

    if values_dict and track_name:
      self.setDataTrack(track_name, EXTERNAL, region_dict, values_dict, color='#8000FF', threshold=0.0)      
    
    return region_dict, values_dict 
    
     
  def calcDataTrackSpatialDensity(self, track_name, source=None, out_track=None,
                                  structure=None, min_exp_value=1e-9, power=2,
                                  null_samples=25, min_seq_sep=21000):
    
    from cUtil import dataLayer, apiUtil
    from scipy.stats import rankdata
    
    dataGroup, source = self._getDataTrackSource(track_name, source)
    
    if not out_track:
      out_track = track_name + '_density'

    refLayer = self.getRefDataTrackGroup(source, track_name)
    if not refLayer:
      return
    
    coordsGroup  = self._getCoordsGroup(structure)
    particGroup  = self._getParticleGroup(structure)
    chromosomes  = [x for x in refLayer if x in particGroup and x in coordsGroup]
    particle_sep = self.getBackboneSpacing(structure)
    color = self.getDataTrackColor(source, track_name).rgb

    all_coords    = []
    all_seq_pos   = []
    track_coords  = []
    track_seq_pos = []
    track_values  = []
    chromo_ranges = []
    
    track_coords_n  = [[] for i in range(null_samples)]
    track_seq_pos_n = [[] for i in range(null_samples)]
    track_values_n  = [[] for i in range(null_samples)]
 
    a = 0
    max_point = 0
    
    for chromo in chromosomes:
 
      coords = self.getChromoCoords(chromo)[0]
      pos = self.getChromoSeqPos(chromo)
      all_coords.append(coords)
      all_seq_pos.append(pos)
      
      b = a + len(coords)
      chromo_ranges.append((a,b))
      a = b
      
      pos_0 = pos[0] - 0.5 * particle_sep
      pos_1 = pos[-1] - 2.5 *  particle_sep
      p_max = pos[-1]
      
      exp_values  = self.getDataTrackValues(track_name, chromo, minValue=min_exp_value)[:,1]
      exp_regions = self.getDataTrackRegions(track_name, chromo, minValue=min_exp_value)
    
      middles = exp_regions.mean(axis=1).astype(int)
      exp_regions[:,0] = middles
      exp_regions[:,1] = middles + 1
     
      binned = dataLayer.regionBinValues(exp_regions, exp_values, int32(particle_sep), pos_0, pos_1)
      idx_a = binned.nonzero()
      
      exp_coords = coords[idx_a]
      exp_points = pos[idx_a].astype(int32)

      track_seq_pos.append(exp_points.astype(int) + max_point)
      track_coords.append(exp_coords)
      track_values.append(binned[idx_a])
      
      for i in range(null_samples):
        j = particle_sep * randint(1, len(coords)-1)
        null_regions = array((exp_regions + j) % p_max, int32)
        null_values  = array(exp_values)
      
        binned = dataLayer.regionBinValues(null_regions, null_values, int32(particle_sep), int32(pos_0), int32(pos_1))

        idx_n = binned.nonzero()
        null_coords = coords[idx_n]
        null_points = pos[idx_n].astype(int)

        track_seq_pos_n[i].append(null_points + max_point)
        track_coords_n[i].append(null_coords) 
        track_values_n[i].append(binned[idx_n]) 
      
      if len(exp_points):
        max_point += exp_points.max() + particle_sep * 5
      
    all_seq_pos = concatenate(all_seq_pos, axis=0).astype(int64)
    all_coords = concatenate(all_coords, axis=0)
    
    track_seq_pos = concatenate(track_seq_pos, axis=0).astype(int64)
    track_coords = concatenate(track_coords, axis=0)
    track_values = concatenate(track_values, axis=0)
    track_values = rankdata(track_values, method='max').astype(float)
    track_values /= track_values.max() # Quantile in [0..1]
    
    for i in range(null_samples):
      track_seq_pos_n[i] = concatenate(track_seq_pos_n[i], axis=0).astype(int64)
      track_coords_n[i] = concatenate(track_coords_n[i], axis=0)
      track_values_n[i] = concatenate(track_values_n[i], axis=0)
      track_values_n[i] = rankdata(track_values_n[i], method='max').astype(float)
      track_values_n[i] /= track_values_n[i].max()
        
    track_seq_pos_n = array(track_seq_pos_n)
    track_coords_n = array(track_coords_n)
    track_values_n = array(track_values_n) 
    
    inv_dist_sums = apiUtil.getInvDistSums(all_coords, track_coords,
                                           all_seq_pos, track_seq_pos, min_seq_sep,
                                           track_values, power/2)
    n = float(len(all_coords))
    field = inv_dist_sums/n

    null_field = []
    null_neighbour_counts = []
    
    m = len(track_coords_n)
    for i in range(m):
      inv_dist_sums = apiUtil.getInvDistSums(all_coords, track_coords_n[i],
                                             all_seq_pos, track_seq_pos_n[i], min_seq_sep,
                                             track_values_n[i], power/2)
      null_field.append(inv_dist_sums/n)

    null_field_mean = mean(null_field, axis=0)
    null_field = concatenate(null_field)
                                                    
    density_bias = log2(field/null_field_mean)    
    v_max = density_bias.max()
    
    #from matplotlib import pyplot as plt
    #from numpy import histogram
    
    #density_bias_quant = rankdata(density_bias, method='max')
    #density_bias_quant /= float(density_bias_quant.max())
    
    region_dict = {}
    values_dict = {}
           
    for c, chromo in enumerate(chromosomes):
      a, b = chromo_ranges[c]
      
      starts = all_seq_pos[a:b] - particle_sep/2
      ends = starts + (particle_sep - 1)
      values = density_bias[a:b]
      
      idx = (values > 0).nonzero()
      starts = starts[idx]
      ends   = ends[idx]
      values = values[idx]
      
      values_norm = values/v_max
      
      #if chromo in ('1','3','7','12'):
      #  hist, edges = histogram(values_norm, bins=25, normed=True)
      #  
      #  print chromo, values_norm.shape
      #  
      #  plt.plot(edges[:-1], hist, label=chromo)
 
      
      region_dict[chromo] = array([starts, ends], int32).T
      values_dict[chromo]  = array([values, values_norm]).T      
  
    #plt.legend()
    #plt.show()
    
    if values_dict:
      self.setDataTrack(out_track, EXTERNAL, region_dict, values_dict, color=color, threshold=0.0)      
    
    return region_dict, values_dict

 
  
  def calcDataTrackMinDist(self, code, source=None, out_track=None, structure=None, nuc=None):
    
    from scipy.spatial import cKDTree
    
    if not out_track:
      out_track = 'dist_' + code
    
    if not nuc:
      nuc = self
    
    dataGroup, source = nuc._getDataTrackSource(code, source)
    
    coordsGroup = self._getCoordsGroup(structure)
    particGroup = self._getParticleGroup(structure)
    models = self.getModels(structure)    

    refLayer = nuc.getRefDataTrackGroup(source, code)
    if not refLayer:
      return
      
    chromosomes = [x for x in refLayer if x in particGroup]
    value_dict = {}
    region_dict = {}
    min_dists = {}
    
    start_coords = [[] for x in models]
    end_coords = [[] for x in models]
    chromo_ranges = {}

    for chromo in chromosomes:
      if chromo not in coordsGroup:
        continue

      if chromo not in particGroup:
        continue
 
      regions = array(refLayer[chromo]['regions'], int32) # start, end
      if not len(regions):
        continue
      
      starts = regions[:,0]
      ends = regions[:,1]
      n = len(starts)
         
      if not n:
        continue

      for m in models:
        start_coords[m].append(self.getPositionCoords(m, starts, chromo, structure, clipEnds=True))
        end_coords[m].append(self.getPositionCoords(m, ends, chromo, structure, clipEnds=True))
     
    for chromo in chromosomes:
       min_dists[chromo] = []
        
    for m in models:
      coords_a = concatenate(start_coords[m], axis=0)
      coords_b = concatenate(end_coords[m], axis=0)
    
      kt_a = cKDTree(coords_a, 10)
      kt_b = cKDTree(coords_b, 10)
  
      for chromo in chromosomes:
      
        coords = self.getChromoCoords(chromo, structure, model=m)
        
        if coords is None:
          continue
        
        if len(coords) < 2:
          continue
        
        dists_a, idx = kt_a.query(coords, k=1)
        dists_b, idx = kt_b.query(coords, k=1)
        
        min_dists[chromo].append(dists_a)
        min_dists[chromo].append(dists_b)
  
    for chromo in min_dists:
      if not len(min_dists[chromo]):
        continue
    
      dists = vstack(min_dists[chromo])
      dists = dists.mean(axis=0)
      
      pos = array(particGroup[chromo]['positions'])
      region_dict[chromo] = array([pos, pos+1], int32).T
      value_dict[chromo]  = array([dists, dists/dists.max()]).T    
  
    if value_dict:
      self.setDataTrack(out_track, DERIVED, region_dict, value_dict)      
    
    return region_dict, value_dict
  
  
  def calcContactDirectionality(self, group_name=None, track_name='directionality', bin_size=50000, width=5e6):
  
    if not group_name:
      group_name = self.getDefaultContactGroup()
   
    if not bin_size:
      bin_size = self.getContactsBinSize(group_name)
        
    contact_cache_dict = self.getCachedContacts(group_name)
    
    chromos = sorted(contact_cache_dict.keys())
       
    value_dict  = dict([(c,[]) for c in chromos])
    region_dict = dict([(c,[]) for c in chromos])
 
    m = int(width/bin_size)
    
    for chromo in chromos:
 
      if chromo not in contact_cache_dict[chromo]:
        continue
      
      start, end = self.getChromosomeLimits(chromo)
      
      obs = self.getContactMatrix(chromo, chromo, bin_size, group_name).astype(float)
      n = len(obs)
      rev = zeros(n, float)
      fwd = zeros(n, float)
      
      for i in range(1,n-1):
        a = max(0, i-m)
        b = min(n, i+m)
        rev[i] = obs[i,a:i].sum()
        fwd[i] = obs[i:i+1:m].sum()
      
      nz = (rev *fwd).nonzero()
      values = zeros(n)

      a = fwd[nz]
      b = rev[nz]
      
      values[nz] = a * log(a/b) + b * log(b/a)
      values -= values.min()
      
      pos = start + bin_size * arange(n)
      values = values[nz]
      #vales = clip(values, 0.0, 1.0)
      pos = pos[nz]
      
      d_max = values.max()
      
      region_dict[chromo] = array([pos, pos+(bin_size-1)], int32).T
      value_dict[chromo]  = array([values, values/d_max]).T    
    
    if value_dict:
      self.setDataTrack(track_name, DERIVED, region_dict, value_dict)      
    
    return region_dict, value_dict
    
    
    
   
  def calcContactDistances(self, groupName=None, structure=None, min_num=1, trackName='cont_dists'):
  
    if not groupName:
      groupName = self.getDefaultContactGroup()
  
    chromos = self.getChromosomes()
    contDict = self.getCachedContacts(groupName)
    models = range(self.getNumModels(structure))
    particGroup = self._getParticleGroup(structure)
    
    chromo_idx = dict([(c,i) for i, c in enumerate(chromos)])
    n_chromo = len(chromos)
    n_cols = n_chromo
    n_rows = n_chromo
    data = []
    labels = []
    d_max = 0.0
    value_dict = dict([(c,[])  for c in chromos])
    region_dict = dict([(c,[])  for c in chromos])
 
    for chr_a in chromos:
      if chr_a not in contDict:
        continue
 
      if len(particGroup[chr_a]['positions']) < 3:
        continue
 
      for chr_b in chromos:
        if chr_b not in contDict[chr_a]:
          continue
 
        if len(particGroup[chr_b]['positions']) < 3:
          continue

        contacts = contDict[chr_a][chr_b]
 
        if contacts.shape[1] < min_num:
          continue
 
        chrPosA = [(chr_a, pos) for pos in contacts[0]]
        chrPosB = [(chr_b, pos) for pos in contacts[1]]

        dists = self.getPositionDistances(chrPosA, chrPosB, models, structure=structure)
 
        if dists is None:
          continue
        
        d_max = max(d_max, dists.max())
        nz = dists.nonzero()[0] # not in same bead
        
        region_dict[chr_a].append(contacts[0,nz])
        region_dict[chr_b].append(contacts[1,nz])
        
        value_dict[chr_a].append(dists[nz])
        value_dict[chr_b].append(dists[nz])
   
    for chromo in list(value_dict.keys()):
      if not region_dict[chromo]:
        del region_dict[chromo]
        del value_dict[chromo]
        continue
    
      pos    = concatenate(region_dict[chromo])
      values = concatenate(value_dict[chromo])
      idx    = pos.argsort()
           
      pos    = pos[idx]
      values = values[idx]
      
      region_dict[chromo] = array([pos, pos+1], int32).T
      value_dict[chromo]  = array([values, values/d_max]).T    
    
    if value_dict and trackName:
      self.setDataTrack(trackName, DERIVED, region_dict, value_dict)      
    
    return region_dict, value_dict
    
  
  def calcRadGyration(self, chromosomes, window=11, structure=None, label='rad_of_gyration'):
    
    from cUtil.apiUtil import calcChainRadGyration
   
    particGroup = self._getParticleGroup(structure)
    valueDict = {}
    regionDict = {}
    models = list(range(self.getNumModels(structure)))
    
    for chromo in chromosomes:
      if chromo not in particGroup:
        continue
      
      rgList = []
      
      for model in models:
        coords = self.getModelCoords(model, [chromo,], structure)
        
        if coords is None:
          continue
        
        nCoords = len(coords)
        if nCoords < 2:
          continue

        rgs = calcChainRadGyration(coords, window)
        
        # Clip big distortions (structure calculation failures)
        med = median(rgs)
        std = 1.4826 * median(abs(rgs-med))
        rgs = clip(rgs, 0.0, med+5*std)
        
        rgList.append(rgs)
      
      if not rgList:
        continue
      
      rgs = array(rgList, float32).mean(axis=0)
      
      rgs_norm = rgs - rgs.min()
      rgs_norm /= rgs_norm.max() or 1.0
      
      positions = array(particGroup[chromo]['positions'])
      sep = int(positions[-1]-positions[0])/float(len(positions))
      
      regionDict[chromo] = vstack([positions-sep, positions+sep]).T
      valueDict[chromo] = vstack([rgs, rgs_norm]).T
      
    if valueDict:
      self.setDataTrack(label, DERIVED, regionDict, valueDict)    
    
    
  def calcDihedralAngles(self, chromosomes, structure=None, label='dihedrals'):
    
    from cUtil.apiUtil import getDihedralAngles
   
    particGroup = self._getParticleGroup(structure)
    valueDict = {}
    regionDict = {}
    angleDict = {}
    
    for chromo in chromosomes:
      if chromo not in particGroup:
        continue
      
      coords = self.getChromoCoords(chromo, structure=structure)
      
      if (coords is None) or (len(coords[0]) < 4):
        continue
      
      angles = getDihedralAngles(coords).mean(axis=0)
      angleDict[chromo] = array(angles)
      
      if label:
        angles_norm = abs(angles)/PI
        angles = angles % TAU
        positions = array(particGroup[chromo]['positions'])
        regionDict[chromo] = vstack([positions[1:-3], positions[2:-2]]).T
        valueDict[chromo] = vstack([angles, angles_norm]).T
        
    if valueDict and label:
      self.setDataTrack(label, DERIVED, regionDict, valueDict)    
    
    return angleDict
    
    
  def calcChromoIntermingling(self, chromosomes, radius=2.0, nIntersect=8, models=None,
                              structure=None, label='intermingled'):
  
    from cUtil.apiUtil import calcSurfaceBuried
   
    particGroup = self._getParticleGroup(structure)
    
    if not models:
      models = list(range(self.getNumModels(structure)))
    
    nBuried = 0
    nTotal = 0
    
    buriedPoints = {}
    for chromo in chromosomes:
      n = len(particGroup[chromo]['positions'])
      buriedPoints[chromo] = zeros(n, int)
    
    for model in models:
      chrRanges = []
      chromos = []
      allCoords = []
 
      n = 0
      for chromo in chromosomes:
        if chromo not in particGroup:
          continue
 
        coords = self.getModelCoords(model, [chromo,], structure)
        
        if coords is None:
          continue
        
        nCoords = len(coords)
 
        if nCoords < 2:
          continue
         
        chrRanges.append( (n, n + nCoords) )
        chromos.append(chromo)
 
        n += nCoords
        allCoords.append(coords)
      
      buried = zeros(n, int)
      allIdx = list(range(n))
      
      for i, chromo in enumerate(chromos):
        a, b = chrRanges[i]
        origIdx = array(allIdx[:a] + allIdx[b:]) # Chop out query chromo
        
        surfCoords = allCoords[i]
        testCoords = vstack(allCoords[:i] + allCoords[i+1:])
                        
        buriedIdx = calcSurfaceBuried(surfCoords, testCoords, radius, 200, nIntersect=nIntersect)
        
        origIdx = origIdx[buriedIdx] # Original indices of buried indices
        buried[origIdx] = 1
      
      for i, chromo in enumerate(chromos):
        a, b = chrRanges[i]
        buriedPoints[chromo] += buried[a:b]
      
      nBuried += len(buried.nonzero()[0])
      nTotal += n
    
    chromoCounts = {}
    
    if label:
      valueDict = {}
      regionDict = {}
      m = len(models) / 2
 
      for chromo in chromos:
 
        idx = (buriedPoints[chromo] > m).nonzero()
        positions = array(particGroup[chromo]['positions'])
        positions = positions[idx]
        p = len(positions)
 
        regionDict[chromo] = vstack([positions, positions+1.0]).T
        valueDict[chromo] = ones((p, 2), float)
        
        chromoCounts[chromo] = p/float(len(buriedPoints[chromo]))
 
      if valueDict:
        self.setDataTrack(label, DERIVED, regionDict, valueDict)
    
    else:
      m = len(models) / 2
      for chromo in chromos:
        idx = (buriedPoints[chromo] > m).nonzero()
        chromoCounts[chromo] = len(idx[0])/float(len(buriedPoints[chromo]))
        
    return float(nBuried)/nTotal, chromoCounts
  
  
  def getModelDepths(self, models, chromosomes, radius=5.0, nVoxels=100, structure=None,
                     separateChromos=False, transInterface=False, filterTrack=None, filterNuc=None):
    
    from cUtil.apiUtil import calcCoordDepth, pointRegionsIntersection
    
    particGroup = self._getParticleGroup(structure)
    chromo_ranges = {}
    
    if filterNuc is None:
      filterNuc = self
    
    if not hasattr(models, '__iter__'):
      models = [models]
       
    model_depths = []
    
    for model in models:
    
      n = 0
      allCoords = []
      chromo_sizes = []
      for chromo in chromosomes:
        if chromo not in particGroup:
          continue
 
        coords = self.getModelCoords(model, [chromo,], structure)
        
        if filterTrack:
          filter_regions = filterNuc.getDataTrackRegions(filterTrack, chromo)
          pos = array(particGroup[chromo]['positions'], int32)
          
          print model, chromo, filterTrack
          
          idx = pointRegionsIntersection(pos, filter_regions)
          coords = coords[idx]
          
        nCoords = len(coords)
 
        if nCoords < 2:
          continue
        
        chromo_sizes.append(nCoords)
        chromo_ranges[chromo] = (n, n + nCoords)
 
        n += nCoords
        allCoords.append(coords)
 
      if transInterface:
        allCoords = vstack(allCoords)
        allDepths = calcCoordDepth(allCoords, radius, nVoxels, chromo_sizes, chromoDepth=False)
 
      elif separateChromos:
        allCoords = vstack(allCoords)
        allDepths = calcCoordDepth(allCoords, radius, nVoxels, chromo_sizes, chromoDepth=True)

      else:
        allCoords = vstack(allCoords)
        allDepths = calcCoordDepth(allCoords, radius, nVoxels)
      
      model_depths.append(array(allDepths))
    
    model_depths = array(model_depths)
  
    return model_depths, chromo_ranges
   
   
  def calcDepths(self, models, chromosomes, radius=5.0, nVoxels=100, structure=None,
                 separateChromos=False, transInterface=False, label=None, filterTrack=None, filterNuc=None):
    """
    (Re)calculate the depth below the surface for particles in
    a 3D chromosomes model. Option to treat chromsomes separetely.
    """
   
    from cUtil.apiUtil import pointRegionsIntersection
    
    if filterNuc is None:
      filterNuc = self
    
    valueDict = {}
    regionDict = {}
    particGroup = self._getParticleGroup(structure)
    
    model_depths, chromo_ranges = self.getModelDepths(models, chromosomes, radius, nVoxels, structure,
                                                      separateChromos, transInterface, filterTrack, filterNuc)
   
    mean_depths = model_depths.mean(axis=0)
    
    dMax = mean_depths.max()
    dMin = mean_depths.min()
    
    for chromo in chromo_ranges:
      start, end = chromo_ranges[chromo]
      depths = mean_depths[start:end]
      
      positions = array(particGroup[chromo]['positions'], int32)

      if filterTrack:
        filter_regions = filterNuc.getDataTrackRegions(filterTrack, chromo)        
        idx = pointRegionsIntersection(positions, filter_regions)
        positions = positions[idx]
      
      valueArray = vstack([depths, (depths-dMin)/(dMax-dMin)]).T
      regionArray = vstack([positions, positions+1.0]).T
    
      regionDict[chromo] = regionArray
      valueDict[chromo] = valueArray
      
    if valueDict:
      if transInterface and not label:
        label = 'trans_dist'
  
      elif separateChromos and not label:
        label = 'chromoDepth'
 
      elif not label:
        label ='depth'
        
      self.setDataTrack(label, DERIVED, regionDict, valueDict)
   
   
  def calcDensity(self, model, chromosomes, radius=100.0, structure=None):
    """(Re)calculate the spatial density of model particles"""
    
    from cUtil.apiUtil import calcCoordDensity, getInvDistSums

    particGroup = self._getParticleGroup(structure)
    valueDict = {}
    regionDict = {}
    allCoords = []
    allPos = []
    chrRanges = []
    chromos = []
    
    n = 0
    for chromo in chromosomes:
      if chromo not in particGroup:
        continue
    
      coords = self.getModelCoords(model, [chromo,], structure)
      
      if coords is None:
        continue
      
      nCoords = len(coords)
      chrRanges.append((n, n + nCoords))
      chromos.append(chromo)
      
      n += nCoords
      allCoords.append(coords)
      allPos.append( array(particGroup[chromo]['positions'], int) ) 
    
    allCoords = vstack(allCoords)
    allPos = hstack(allPos)
    #allDensities = calcCoordDensity(allCoords, radius)
    
    allDensities = getInvDistSums(allCoords, allCoords, allPos, allPos, 0)
    
    std = allDensities.std()
    mean = allDensities.mean()
    
    minDensity = mean - 3*std
    maxDensity = mean + 3*std
    
    allDensities = clip(allDensities, minDensity, maxDensity)
    
    for i, chromo in enumerate(chromos):
      start, end = chrRanges[i]
      
      densities = allDensities[start:end]
      
      densNorm = (densities - minDensity) / (maxDensity - minDensity)
      positions = array(particGroup[chromo]['positions'])
      
      valueArray = vstack([densities, densNorm]).T
      regionArray = vstack([positions, positions+1.0]).T
    
      regionDict[chromo] = regionArray
      valueDict[chromo] = valueArray
      
    if valueDict:
      self.setDataTrack('density', DERIVED, regionDict, valueDict)
  
  
  def calcCentreDists(self, chromosomes, structure=None, label='centre_dist'):

    particGroup = self._getParticleGroup(structure)
    nuc_coords = []
    chromo_idx_ranges = []
    chromos = []
    
    n = 0
    for chromo in chromosomes:
      if chromo not in particGroup:
        continue
    
      coords = self.getChromoCoords(chromo, structure=structure)
      
      if coords is None:
        continue
      
      nCoords = len(coords[0])
      chromo_idx_ranges.append((n, n + nCoords))
      chromos.append(chromo)
      
      n += nCoords
      nuc_coords.append(coords)
    
    nuc_coords = concatenate(nuc_coords, axis=1)
    all_dists = []
    m = len(nuc_coords)
    
    for i in range(m):
      model_coords = nuc_coords[i]
      centre = model_coords.mean(axis=0)
      
      deltas = model_coords - centre
      dists = sqrt((deltas * deltas).sum(axis=1))
      all_dists.append(dists)
     
    all_dists = vstack(all_dists)
    all_dists = all_dists.mean(axis=0)
      
    d_max = all_dists.max()
    
    valueDict = {}
    regionDict = {}
    for i, chromo in enumerate(chromos):
      start, end = chromo_idx_ranges[i]
      
      dists = all_dists[start:end]
      dists_norm = dists / d_max
      
      positions = array(particGroup[chromo]['positions'])
      
      valueArray = vstack([dists, dists_norm]).T
      regionArray = vstack([positions, positions+1.0]).T
    
      regionDict[chromo] = regionArray
      valueDict[chromo] = valueArray
      
    if valueDict and label:
      self.setDataTrack(label, DERIVED, regionDict, valueDict)
    
    return regionDict, valueDict
    
    
  
  def calcRestraintViolations(self, chromosomes, models=None, cis=True,
                              trans=True, upperOnly=False, reportAll=True,
                              usePositions=True, structure=None):
    """Calculate 3D model distances of directly restrained particles and compare to restraint bounds"""
    
    if not models:
      models = range(self.getNumModels(structure))
    
    restrDict = self.getRestraints(chromosomes, cis, trans, usePositions=usePositions, structure=structure)
    violDict = {}
    
    for key in restrDict: # pos0, pos1, weight, target, lower, upper
      chromoA, chromoB = key
      restraints = restrDict[key].T
      if not len(restraints[0]):
        continue
      
      posA = array(restraints[0], uint32)
      posB = array(restraints[1], uint32)
       
      if usePositions:
        dists = self.getPositionDistances(posA, posB, models, chromoA, chromoB, structure)
        
      else:
        dists = self.getCoordDistances(chromoA, chromoB, models, structure=structure)
      
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
    
    
  def calcRestraintDistances(self, chromosomes, models=None, cis=True, trans=True, structure=None):
    """Calculate 3D model distances of directly restrained particles"""
    
    if not models:
      models = range(self.getNumModels(structure))

    restrDict = self.getRestraints(chromosomes, cis, trans, usePositions=True, structure=structure)
    distDict = {}
    
    for key in restrDict: # pos0, pos1, weight, target, lower, upper
      chromoA, chromoB = key
      restraints = restrDict[key]
      posA = array(restraints[0], uint32)
      posB = array(restraints[1], uint32)
      
      chrPosA = [(chromoA, pos) for pos in posA]
      chrPosB = [(chromoB, pos) for pos in posB]
    
      dists = self.getPositionDistances(chrPosA, chrPosB, models, structure=structure)
      distDict[key] = dists
      
      
    return distDict
  
  
  def calcBackboneDistances(self, structure=None, model=None, trackName='backbone_dist'):

    coordsGroup = self._getCoordsGroup(structure)
    particGroup = self._getParticleGroup(structure) 
       
    if model:
      models = [model]
    
    else:  
      models = range(self.getNumModels(structure))
    
    dist_dict = {}
      
    for chromo in coordsGroup:
      coords = array(coordsGroup[chromo])
      model_dists = []

      for model in models:
        deltas = coords[model,1:] - coords[model,:-1]
        dists = sqrt((deltas*deltas).sum(axis=1))
        model_dists.append(dists)
        
      model_dists = vstack(model_dists).mean(axis=0)
      dist_dict[chromo] = model_dists.astype(float32)

      
    if trackName:
      region_dict = {}
      value_dict = {}
    
      for chromo in dist_dict:
      
        pos = array(particGroup[chromo]['positions'], int32)
        vals = dist_dict[chromo]
        
        norm_vals = vals - vals.min()
        norm_vals /= norm_vals.max()
                
        value_dict[chromo]  = vstack([vals, norm_vals]).T
        region_dict[chromo] = vstack([pos[:-1], pos[1:]-1]).T

      self.setDataTrack(trackName, DERIVED, region_dict, value_dict)
    
    return dist_dict
    
      
  def calcModelRmsds(self, models=None, chromosomes=None, backbone=None, weightThreshold=1.0, structure=None):
    """Calculate RMSD of structural models as a measure of coordinate precision.
       Backbone True:only, False:exclude, None:everything
       """
       
    from util.Structure import superimposeCoordArray
 
    # Single value for bundle
   
    if not models:
      models = range(self.getNumModels(structure))
    
    if not models:
      return [], []
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    ensemble = [self.getModelCoords(m, chromosomes, structure, backbone) for m in models]
    ensemble, rotations, rmsds, atomRmsds = superimposeCoordArray(ensemble, weightThreshold)
 
    return rmsds, atomRmsds


  def calcRmsdDataTrack(self, chromosomes=None, weightThreshold=1.0, structure=None, trackName='RMSD'):
    """Calculate RMSD of structural models along the chromosomal sequence
    """
       
    from util.Structure import superimposeCoordArray

    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
   
    particGroup = self._getParticleGroup(structure)
    regionDict = {}
    valueDict = {}
    
    for chromo in chromosomes:
    
      positions = array(particGroup[chromo]['positions'])
      
      if (positions is None) or not len(positions):
        continue
      
      ensemble = self.getChromoCoords(chromo, structure=structure)
      
      if (ensemble is None) or not len(ensemble):
        continue
      
      ensemble, rotations, rmsds, atomRmsds = superimposeCoordArray(ensemble, weightThreshold)
 
      atomRmsdsNorm = atomRmsds / (atomRmsds.max() or 1.0)

      valueArray = vstack([atomRmsds, atomRmsdsNorm]).T
      regionArray = vstack([positions, positions+1.0]).T
 
      regionDict[chromo] = regionArray
      valueDict[chromo] = valueArray

    self.setDataTrack(trackName, DERIVED, regionDict, valueDict)
    

  def calcModelRmsdMatrix(self, models=None, backbone=None, weightThreshold=1.0, structure=None):
    """Calculate a pairwise model RMSD matrix
       Backbone True:only, False:exclude, None:everything"""
    
    from util.Structure import superimposeCoordPair
    
    if not models:
      models = range(self.getNumModels(structure))
    
    chromosomes = self.getChromosomes(structure)
    ensemble = [self.getModelCoords(m, chromosomes, structure, backbone) for m in models]
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
    
  
  def clusterStructures(self, structures=None, weightThreshold=1.0):
    """
    Hierarchically cluster structures
    """
    
    from scipy.spatial import distance
    from scipy.cluster import hierarchy
    from util.Structure import superimposeCoordPair
    
    if not structures:
      structures = [s for s in self.structures]
    
    # Does any coord interpolation
    ensemble = self.mergeStructures(structures, target=None, nModels=None, storeCoords=False)
    
    n = len(ensemble)
    matrix = zeros((n,n), float)
    
    for i in range(n-1):
      coords1 = ensemble[i]
    
      for j in range(i+1, n):
        coords2 = ensemble[j]
        coords3, rot, rmsd, atomRmsds = superimposeCoordPair(coords1, coords2, weightThreshold)
        
        matrix[i,j] = rmsd
        matrix[j,i] = rmsd
 
    matrix = distance.squareform(matrix)
    
    linkage = hierarchy.average(matrix) # (idxA, idxB, dist, nChildren)
     
    #from matplotlib import pyplot    
    #pyplot.figure()
    #pyplot.title('Structure hierarchy')        
    #graph = hierarchy.dendrogram(linkage)
    
    return linkage
    
    
  def calcModelClusters(self, k=4, backbone=None, structure=None):
    """Group structural models into a given mumber of clusters"""
    
    from util.Cluster import kMedioids, hierarchicalCluster
    
    rmsdMatrix = matrixOrig = self.calcModelRmsdMatrix(backbone=backbone, structure=structure)
    
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
        regions = array(dataLayerA[chromo]['regions'], int32) # start, end
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
          
  
  def modelCentre(self, models=None, chromosomes=None, structure=None):
    """Centre chromosome model coordinates at the orgin"""
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    if not models:
      models = range(self.getNumModels(structure))
    
    centers = []
    
    for model in models:
      coords = self.getModelCoords(model, chromosomes, structure)
      center = coords.sum(axis=0)/float(len(coords))
      coords -= center
      
      self.setModelCoords(coords, model, chromosomes, structure)
      centers.append(center)
      
    return centers
  
  def mergeStructures(self, structures=None, target=None, nModels=None, storeCoords=True):
    """
    Make a structure by combining existing chromosome models into a larger ensemble,
    interpolationg coords into a unified set of positions.
    """

    from cUtil.apiUtil import interpolateChromoModelCoords
    
    if not structures:
      structures = list(self.structures.keys())
    
    chromos = set()
    structs = []
    
    # Collate chromosomes and valid structures
    
    for s in structures:
      if not self.getNumModels(s):
        continue
    
      particGroup = self._getParticleGroup(s)
      avail = set(particGroup.keys())
      
      if not avail:
        continue
      
      if chromos:
        chromos &= avail
        
      else:
        chromos.update(avail)
        
      structs.append(s)

    if chromos:
      chromos = sorted(chromos)
      
    else:
      msg = 'No common set of chromosomes for selected structures'
      self._warning(msg)
      return
    
    # Get particle positions for interpolation
    
    posDictMerge = {}
    posDict = {}
    
    for s in structs:
      posDict[s] = {}
      particGroup = self._getParticleGroup(s)
      
      for chromo in chromos:
        positions = array(particGroup[chromo]['positions'], int32)
        
        if len(positions) > 2:
          if chromo not in posDictMerge:
            posDictMerge[chromo] = set()
        
          posDictMerge[chromo].update(positions) # combined positions for interpolation
          posDict[s][chromo] = positions
    
    for chromo in posDictMerge:
       posDictMerge[chromo] = array(sorted(posDictMerge[chromo]), int32)
    
    
    # Create new ensemble of models by interpolation
      
    ensemble = []
    for s in structs:
      n = self.getNumModels(s)
      
      if nModels:
        n = min(nModels, n)
      
      for m in range(n):
        coords = self.getModelCoords(m, chromos, s)
        
        if len(coords) > 2:
          coords = interpolateChromoModelCoords(posDictMerge, posDict[s], coords)
          ensemble.append(coords)
   
    ensemble = array(ensemble)
    
 
    # Set merged positions and coordinates
    
    if storeCoords:
      if not target:
        target = self.getNextStructureCode()
 
      self.addChromosomes(posDictMerge, structure=target)
      self.setAllCoords(ensemble, chromos, structure=target)
    
    return ensemble
   

  def structureAlign(self, structures=None, chromosomes=None, backbone=None, weightThreshold=10.0):
    """
    Superpose all models of stated structures using RMSD weighted iterative SVD
    Backbone True:only, False:exclude, None:everything
    """

    from util.Structure import superimposeCoordArray
    from cUtil.apiUtil import interpolateChromoModelCoords
    
    if not structures:
      structures = list(self.structures.keys())
       
    nModels = None
    structs = []
    for s in structures:
      m = self.getNumModels(s)
      
      if not m:
        continue
      
      if nModels:
        if m < nModels:
          nModels = m
        
      else:
        nModels = m
        
      particGroup = self._getParticleGroup(s)
      
      if not chromosomes:
        chromosomes = list(particGroup.keys())
        structs.append(s)
        
      else:
        for chromo in chromosomes:
          if chromo not in particGroup:
            break
             
        else:
          structs.append(s)
      
    if len(structs) < 2:
      msg = 'No structures with matching chromosomes (and at least one model) to align'
      self._warning(msg)
      return
    
    # Get particle positions for interpolation to universal set
        
    posDictMerge = {}
    for chromo in chromosomes:
      posDictMerge[chromo] = set()
    
    posDict = {}
    for s in structures:
      posDict[s] = {}
      particGroup = self._getParticleGroup(s)
      
      for chromo in chromosomes:
        positions = array(particGroup[chromo]['positions'], int32)
        posDictMerge[chromo].update(positions) # combined positions for interpolation
        posDict[s][chromo] = positions
    
    for chromo in posDictMerge:
       posDictMerge[chromo] = array(sorted(posDictMerge[chromo]), int32)
    
    # Interpolate coords and build ensemble
    
    ensemble = []
    for s in structures:
      for m in range(nModels):
        coords = self.getModelCoords(m, chromosomes, s, backbone)
        coords = interpolateChromoModelCoords(posDictMerge, posDict[s], coords)
        ensemble.append(coords)
    
    # SVD alignment
    
    ensemble, rotations, rmsds, atomRmsds = superimposeCoordArray(array(ensemble), weightThreshold)
    
    
    # Apply alignment transforms to original coords

    i = 0
    for s in structures:
      for m in range(nModels):
        coords = self.getModelCoords(m, chromosomes, s)
        coords -= coords.mean(axis=0)
        coords = dot(coords, rotations[i])
        self.setModelCoords(coords, m,  chromosomes, s)
        i += 1 
  
  
  def modelAlign(self, models=None, chromosomes=None, backbone=None, weightThreshold=10.0, structure=None):
    """
    Superpose the chromosome models of one structure using RMSD weighted iterative SVD
    Backbone True:only, False:exclude, None:everything
    """

    from util.Structure import superimposeCoordArray
      
    if not models:
      models = range(self.getNumModels(structure))
 
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    ensemble = [self.getModelCoords(m, chromosomes, structure, backbone) for m in models]
    
    ensemble, rotations, rmsds, atomRmsds = superimposeCoordArray(ensemble, weightThreshold)
    
    n_models = len(models)
    
    """
    for c in ['9',]:
      for s1 in self.structures:
        pg1 = self._getParticleGroup(s1)
 
        for s2 in self.structures:
          pg2 = self._getParticleGroup(s2)
 
          pos1 = array(pg1[c]['positions'])
          pos2 = array(pg2[c]['positions'])
 
          a2 = 0
          b2 = len(pos2)
          while pos1[0]-pos2[a2] > 1e5:
            pos2 = pos2[1:]
            a2 += 1
 
          a1 = 0
          b1 = len(pos1)
          while pos2[0]-pos1[a1] > 1e5:
            pos1 = pos1[1:]
            a1 += 1
 
          while (b2-a2) > (b1-a1):
            b2 -= 1

          while (b1-a1) > (b2-a2):
            b1 -= 1
 
          pair_rmsds = []
 
          for m in models:
            ensemble = [self.getModelCoords(m, [c], s1)[a1:b1], self.getModelCoords(m+1 % n_models, [c], s2)[a2:b2]]
            ensemble, rotations, rmsds, atomRmsds = superimposeCoordArray(ensemble, weightThreshold)
            pair_rmsds.append(rmsds[0])
 
          pair_rmsds = array(pair_rmsds)
          print c, s1, s2, pair_rmsds.mean()
    
    """
     
    for i, modelCoords in enumerate(ensemble):
      self.setModelCoords(modelCoords, models[i],  chromosomes, structure)


  def modelAlignAxes(self, models=None, chromosomes=None, structure=None):
    """Align largest orthogonal model dimensions with coordinate axes"""

    from util.Cluster import principleComponentAnalysis
    
    if not models:
      models = range(self.getNumModels(structure))
      
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
      
    coordsGroup = self._getCoordsGroup(structure)   
    
    if coordsGroup:
      for model in models:
        modelCoords = self.getModelCoords(model, chromosomes, structure)
        eigenMat, energy = principleComponentAnalysis(modelCoords)
        modelCoords = dot(modelCoords, eigenMat)
        self.setModelCoords(modelCoords, model, chromosomes, structure)


  def modelScale(self, factor, models=None, chromosomes=None, structure=None):
    """Scale 3D chromosome model coordinates by specified factor"""
  
    matrix = factor * eye(3)
    
    self.modelTransform(matrix, models, chromosomes, structure)
    
    
  def modelRotate(self, angle, axis=(0,0,1), models=None, chromosomes=None, structure=None):
    """Rootate chromosome model coordinates by specified angle around an axis"""
    
    from util.Structure import getRotationMatrix
     
    matrix = getRotationMatrix(axis, angle)
    
    self.modelTransform(matrix, models, chromosomes, structure)
    
    
  def modelTranslate(self, vector, models=None, chromosomes=None, structure=None):
  
    vector = array(vector) 
     
    for model in models:
      modelCoords = self.getModelCoords(model, chromosomes, structure)
      modelCoords += vector
      self.setModelCoords(modelCoords, model, chromosomes, structure)
     
     
  def modelMirror(self, models=None, chromosomes=None, structure=None):
    """Make 3D chromosome model coordinates a mirror image"""

    matrix = eye(3)
    matrix[0,0] = -1
    
    self.modelTransform(matrix, models, chromosomes, structure)
    
    
  def modelTransform(self, matrix, models=None, chromosomes=None, structure=None):
    """Apply a 3D matrixs tranformation to the coordinates of chromosomes in structural models"""

    if not models:
      models = range(self.getNumModels(structure))
      
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
   
    matrix = array(matrix)

    if matrix.shape == (3,3):
      for model in models:
        modelCoords = self.getModelCoords(model, chromosomes, structure)
        modelCoords = dot(modelCoords, matrix)
        self.setModelCoords(modelCoords, model, chromosomes, structure)

    else:
      matShape = 'x'.join(['%d' % x for x in matrix.shape])
      raise Exception('Transformation matrix must be 3x3 not %s' % matShape)
  
  
  def setSpiralCoords(self, radius=10.0, chromosomes=None, structure=None):
    """Set stucture to a spriral, given chromosomal positions"""
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
      
    particGroup = self._getParticleGroup(structure)
    nPoints = sum([len(particGroup[c]['positions']) for c in chromosomes])
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
      self.setAllCoords(modelCoords, chromosomes, structure)
    
  
  def setRandomCoords(self, models=None, chromosomes=None, centre=(0.0,0.0,0.0),
                      randWalk=True, randSeed=None, maxStep=1.0, structure=None):
    """Set the coordinates of chromsome structural models to random positions"""
    
    uniform = random.uniform
    
    if randSeed:
      seed(randSeed)
    else:
      seed()
    
    if not models:
      models = range(self.getNumModels(structure))
      
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    
    centre = array(centre)
    particGroup = self._getParticleGroup(structure)
    nPoints = sum([len(particGroup[c]['positions']) for c in chromosomes])
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
    self.setAllCoords(allCoords, chromosomes, structure)
  
  
  def setRandomContacts(self, groupName, model, numRestraints, chromosomes=None,
                        numNeighbours=10, replace=True, structure=None):
    """Make restraints based upon random close points, 
       e.g. given a synthetic structure"""
    
    from cUtil.apiUtil import getClosestPoints
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    particGroup = self._getParticleGroup(structure)
    coords = self.getModelCoords(model, chromosomes, structure)
    
    start = 0
    chromos = []
    chromoStarts = []
    posDict = {}
    
    coordsGroup = self._getCoordsGroup(structure)

    for chromo in coordsGroup:
      if chromo in chromosomes:
        chromos.append(chromo)
        chromoStarts.append(start)
        pos = array(particGroup[chromo]['positions'])
        posDict[chromo] = pos
        start += len(pos)
    
    n = 0
    contactDict = {}
    append = {}
    numCoords = len(coords)
    done = set()
    doneAdd = done.add
    numObs = 1
    #numRestraints = min(numCoords/2, numRestraints)
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
      
      chrA = chromos[iChr]
      chrB = chromos[jChr]
 
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
                        rad1=0.4, rad2=0.2, contact=0.06, scale=45.0, structure=None):
                    
    """Set the coordinates of chromsome structural models to a synthetic test structure"""
    
    from util.Structure import getCoiledCoilCoords, HEX_GRID
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    particGroup = self._getParticleGroup(structure)
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
      numPoints = len(particGroup[chromo]['positions'])
    
      x0, y0 = hexGrid[i]
      baseCoord = [x0, y0, 0.0]
      coords = getCoiledCoilCoords(baseCoord, numPoints, dAngle1,
                                   dAngle2, rad1, rad2, contact)
      coords *= scale
      
      self.setModelChromosomeCoords(coords, chromo, model, structure)
 
    self.setRandomContacts('singleCell', model, numRestraints, chromosomes, structure=structure)
  
  
  def setTestCoordsHilbert(self, chromosomes=None, centre=(0.0,0.0,0.0), scale=5.0, structure=None):
    """Set the coordinates of chromsome structural models to a synthetic test structure"""
    
    from util.Structure import indexToHilbert3D
    
    particGroup = self._getParticleGroup(structure)
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    model = 0
    chromoSizes = [len(particGroup[c]['positions']) for c in chromosomes]
    nPoints = sum(chromoSizes)
    
    coords = array([indexToHilbert3D(i) for i in range(nPoints)]) * scale
    coords -= coords.mean(axis=0) - array(centre) 
 
    self.setAllCoords(coords, chromosomes, structure)
    self.setRandomContacts('singleCell', model, chromosomes, structure=structure)
  
  
  def setContactMatrix(self, groupName, chromoA, chromoB, matrix, binSize,
                       startA=0, startB=0, isSingleCell=False):
    """Sets the contacts for a specified chromosome pair from a matrix
       with specified binSize. Option to state the starting point (in bp)
       of chromosome for each matrix axis. 
       Only affects the working contacts; originals are preserved.
    """ 
  
    from cUtil.apiUtil import intMatrixToSparse
  
    group = self._getContactGroup(groupName, self.origContacts)
    group.attrs['isSingleCell'] = 1 if isSingleCell else 0
    group.attrs['binSize'] = binSize or 1
   
    n, m = matrix.shape
    matrix = matrix.astype(int32)
    binSize = int32(binSize)
    
    if (chromoA in self.origContacts) and (chromoB in self.origContacts[chromoA]):
      
      contacts = intMatrixToSparse(matrix, binSize, binSize, int32(startA), int32(startB))
      contacts = array(contacts, uint32)
 
      subGroup = self._getGroup(chromoA, group)   
      ds = self._setData(chromoB, subGroup, uint32, contacts, compression='gzip')
 
      if groupName in self._contactsCache:
        if chromoA in self._contactsCache[groupName][chromoA]:
          self._contactsCache[groupName][chromoA][chromoB] = contacts
    
    else:
      contacts = intMatrixToSparse(matrix.T, binSize, binSize, int32(startB), int32(startA))
      contacts = array(contacts, uint32)
      
      subGroup = self._getGroup(chromoB, group)   
      ds = self._setData(chromoA, subGroup, uint32, contacts, compression='gzip')    
           
      if groupName in self._contactsCache:
        if chromoB in self._contactsCache[groupName][chromoB]:
          self._contactsCache[groupName][chromoB][chromoA] = contacts
  
  
  def setContacts(self, groupName, contactDict, replace=True,
                  selected=None, color=None, isSingleCell=True,
                  transposed=True):
    """Sets working contacts, never originals; they are only loaded
       By default any previous contacts will be replaced
    """
    
    if groupName in self.origContacts:
      group = self._getContactGroup(groupName, self.workContacts)
    
    else: # New
      group = self._getContactGroup(groupName, self.origContacts)
      group.attrs['isSingleCell'] = 1 if isSingleCell else 0
      
      if isSingleCell:
        binSize = 1
      else:
        binSize = 1
 
        for key in contactDict:
          if len(contactDict[key]):
            pos = array(contactDict[key])[:,0]
            sz = (pos[1:]-pos[:-1]).min()
            if sz:
              binSize = sz
              break

      group.attrs['binSize'] = binSize
    
    if color is not None:
      group.attrs['color'] = Color(color).rgba
        
    if selected is not None:
      group.attrs['selected'] = selected
      
    for chromoKey in contactDict:
      chromoA, chromoB = chromoKey
      
      if transposed:
        contacts = array(contactDict[chromoKey]).T
      else:
        contacts = array(contactDict[chromoKey])
        
      subGroup = self._getGroup(chromoA, group)
      
      if (chromoB in subGroup) and not replace:
        prev = array(subGroup[chromoB])
        prevDim = prev.shape[1]
        contDim = contacts.shape[1] 
        
        if prevDim < contDim:
          msg  = 'Attempt to add contacts with model indices to a dataset '
          msg += 'without model indices. Model indices will be removed.'
          self._warning(msg)
          contacts = contacts[:3]
          
        elif contDim < prevDim:
          msg =  'Attempt to add contacts without model indices to a dataset '
          msg += 'which requires them'
          raise Exception(msg)
          
        contacts = concatenate([prev, contacts], axis=1)
        
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


  def removeStructure(self, structure):
  
    self._delete(str(structure), self.structures)
  
  
  def removeModels(self, structure, models):
  
    coordsGroup = self._getCoordsGroup(structure)
    
    n = self.getNumModels(structure)
    
    models = [i for i in set(models) if 0 <= i < n]
    
    if len(models) < n:
      idx = [i for i in range(n) if i not in models]
    
      for chromo in coordsGroup:
        coords = array(coordsGroup[chromo])        
        coords = coords[idx]
        self._setData(chromo, coordsGroup, float, coords)
    
    else:
      self._delete(structure, self.structures)

  
  def removeContacts(self, groupName):

    if groupName in self.workContacts:
      self._delete(groupName, self.workContacts)
    
    if groupName in self.origContacts:
      self._delete(groupName, self.origContacts)
   
  
  def removeViolatedContacts(self, groupName, chromosomes=None, threshold=5.0,
                             models=None, structure=None):
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
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
      dists = self.getPositionDistances(chrPosA, chrPosB, models, structure=structure)
      
      if dists is None:
        continue
      
      meanDist = dists.mean()
      
      indices = (dists < threshold).nonzero()
      
      if not group:
        group = self._getContactGroup(groupName, self.workContacts)

      subGroup = self._getGroup(chromoA, group)
      
      if len(indices):
        n +=  len(posA) - indices[0].shape[0]
        data = data[:,indices[0]]
        self._setData(chromoB, subGroup, uint32, data, compression='gzip')
      
      else:
        self._delete(chromoB, subGroup)
      
    if groupName in self._contactsCache:
      del self._contactsCache[groupName]
    
    return n  
  
  
  def resolveAmbiguousContacts(self, groupName, chromosomes=None, ideal_dist=2.0,
                               models=None, structure=None, min_contrib=0.9):
    
    from collections import defaultdict

    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    n = 0
    if models and len(models) == 1:
      model = models[0]
    else:
      model = None
    
    contDict = self.getContacts(groupName, chromosomes, cis=True, trans=True, model=model)
    group = None
     
    for chromoKey in contDict:
      data = contDict[chromoKey]
      
      if len(data) != 4:
        continue
 
      resolved = []
      chromoA, chromoB = chromoKey
      posA = data[0]
      posB = data[1]
      
      ambig_groups = defaultdict(list)
      for i, ambig_group in enumerate(data[3]):
        ambig_groups[ambig_group].append(i)
      
      chrPosA = [(chromoA, pos) for pos in posA]
      chrPosB = [(chromoB, pos) for pos in posB]
      dists = self.getPositionDistances(chrPosA, chrPosB, models, structure=structure)
      
      if dists is None:
        continue
        
      for ambig_group, group_idx in ambig_groups.iteritems():
        group_dists = dists[group_idx]           
       
        if len(group_idx) == 2:
          group_dists /= ideal_dist
          group_dists = 1.0/(group_dists*group_dists)
          contrib = group_dists / (group_dists.sum() or 1.0)
          i_largest = contrib.argsort()[-1]
          
          if contrib[i_largest] >= min_contrib:
             resolved.append(group_idx[i_largest]) 
 
        elif len(group_idx) == 1:
          if group_dists[0] < 2.5 * ideal_dist:
            resolved.append(group_idx[0])       
      
      if not group:
        group = self._getContactGroup(groupName, self.workContacts)

      subGroup = self._getGroup(chromoA, group)
      
      if len(resolved):
        n += len(resolved)
        data = data[:,resolved]
        self._setData(chromoB, subGroup, uint32, data, compression='gzip')
      
      else:
        self._delete(chromoB, subGroup)
      
    if groupName in self._contactsCache:
      del self._contactsCache[groupName]
    
    return n  
  
  
  def mergeContacts(self, name_a, name_b, name=None, remove_frac=None):
  
    from cUtil.apiUtil import binContacts, intMatrixToSparse
    
    cont_dict_a = self.getCachedContacts(name_a)
    cont_dict_b = self.getCachedContacts(name_b)
    
    group_a = self.getContactGroup(name_a)
    group_b = self.getContactGroup(name_b)
   
    bs_a = group_a.attrs['binSize'] or 1
    bs_b = group_b.attrs['binSize'] or 1
    
    if bs_b > bs_a:
      bs_a, bs_b = bs_b, bs_a
      group_a, group_b = group_b, group_a
      cont_dict_a, cont_dict_b = cont_dict_b, cont_dict_a
    
    f = bs_a/float(bs_b)
      
    if abs(f-int(f)) > 1e-12:
      msg = 'Can only merge contacts with unequal bin sizes if one size is an integer multiple of the other'
      raise Exception(msg)  
 
    if not name:
      name = '%s+%s' % (group_a,group_b)
    
    merge_dict = {}
    
    chromos1 = set(cont_dict_a.keys()) | set(cont_dict_b.keys())
    
    for c1 in chromos1:
      s1, e1 = self.getChromosomeLimits(c1)
      merge_dict[c1] = {}
      chromos2 = set(cont_dict_a[c1].keys()) | set(cont_dict_b[c1].keys())
      
      for c2 in chromos2:
        s2, e2 = self.getChromosomeLimits(c2)
        merge_dict[c1][c2] = []
        
        if (c1 in cont_dict_a) and (c2 in cont_dict_a[c1]):
          conts = cont_dict_a[c1][c2]
          
          if remove_frac:
            idx = arange(conts.shape[1])
            random.shuffle(idx)
            idx = idx[:int(remove_frac*len(idx))]
            conts = conts[:,idx]
          
          merge_dict[c1][c2].append(conts)
        
        if (c1 in cont_dict_b) and (c2 in cont_dict_b[c1]):
          conts = cont_dict_b[c1][c2]

          if remove_frac:
            idx = arange(conts.shape[1])
            random.shuffle(idx)
            idx = idx[:int(remove_frac*len(idx))]
            conts = conts[:,idx]
        
          merge_dict[c1][c2].append(conts)
        
        if bs_a == 1:
          data = concatenate(merge_dict[c1][c2], axis=1)
          idx = data[0].argsort()
          data = data[:,idx].T.tolist()
          
          prev_row = data[0]
          uniq_data = [prev_row]
          for row in data[1:]:
            if row[:2] != uniq_data[-1][:2]:
              uniq_data.append(row)
          
          data = array(uniq_data).T
          
        else:
          s_a = int(s1/bs_a)
          s_b = int(s2/bs_a)
          e_a = int(e1/bs_a)
          e_b = int(e2/bs_a)
         
          n = e_a - s_a + 1
          m = e_b - s_b + 1
          
          matrix = zeros((n,m), int32)
          
          for conts in data:
            binContacts(conts, matrix, int32(s1), int32(s2), int32(bs_a))
  
          data = intMatrixToSparse(matrix, int32(bs_a), int32(bs_a), int32(s1), int32(s2))
  
        merge_dict[c1][c2] = data
    
    ca = array(Color(group_a.attrs['color']).rgba)
    cb = array(Color(group_b.attrs['color']).rgba)
    
    group = self._getContactGroup(name, self.workContacts)
    group.attrs['isSingleCell'] = int(group_a.attrs['isSingleCell'] and group_b.attrs['isSingleCell'])
    group.attrs['binSize'] = bs_a
    group.attrs['color'] = (ca+cb)/2.0
    group.attrs['selected'] = 1
    
    for c1 in merge_dict:
      sub_group = self._getGroup(c1, group)
      
      for c2 in merge_dict[c1]:
        self._setData(c2, sub_group, uint32, merge_dict[c1][c2], compression='gzip')
 
        
       
  def splitIntraDataTrackContacts(self, groupName, dataSource, dataCode, targetName=None):
    """
    Sub divide a contact map into those contacts that occur WITHIN
    a specified data track's regions

    """
    
    # TBD: test this function
    # TBD: Add get contact data overlap func to GUI
    
    if not targetName:
      targetName = groupName + '_' + dataCode
    
    chromosomes = self.getChromosomes()
    dataDict = self.getDataTrack(dataSource, dataCode, chromosomes)
    # dataDict[chromo] = (regionData, valueData, strands, annotations)
    contDict = self.getContacts(groupName, chromosomes, cis=True, trans=True)    
    group = None
    
    for chromoKey in contDict:
      chromoA, chromoB = chromoKey
      
      regionsA = sorted(dataDict[chromoA][0])
      regionsB = sorted(dataDict[chromoB][0])
      nr = len(regionsA)
      
      contacts = contDict[chromoKey]
      posA = sorted(contacts[0])
      posB = sorted(contacts[1])
      
      selected = zeros(len(contacts), int)

      j = 0
      for i, p in enumerate(posA):
        p1, p2 = regionsA[j]
        
        while p > p2:
          j += 1
          
          if j == nr:
            break
          
          p1, p2 = regionsA[j]
        
        if p1 < p < p2:
          selected[i] += 1
       
        if j == nr:
          break

      j = 0
      for i, p in enumerate(posB):
        p1, p2 = regionsB[j]
        
        while p > p2:
          j += 1
          if j == nr:
            break
            
          p1, p2 = regionsB[j]
        
        if p1 < p < p2:
          selected[i] += 1
        
        if j == nr:
          break
      
      indices = (selected > 1).nonzero()[0] # Both posA and posB were within a data region
      contacts = contacts[:,indices]
 
      if not group:
        group = self._getContactGroup(targetName, self.workContacts)
 
      subGroup = self._getGroup(chromoA, group)
      self._setData(chromoB, subGroup, uint32, contacts, compression='gzip')
     
    
  def getViolatedContacts(self, groupName, targetGroupName, chromosomes=None,
                          threshold=4.0, models=None, structure=None):
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
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
      dists = self.getPositionDistances(chrPosA, chrPosB, models, structure=structure)
      meanDist = dists.mean()
      
      indices = (dists > meanDist*threshold).nonzero()
      
      if len(indices):
        n +=  indices[0].shape[0]
        data = data[:,indices[0]]
 
        if not group:
          group = self._getContactGroup(targetGroupName, self.workContacts)
 
        subGroup = self._getGroup(chromoA, group)
        self._setData(chromoB, subGroup, uint32, data, compression='gzip')
    
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
    
    
  def setModelChromosomeCoords(self, coords, chromosome, model=None, structure=None):
    """Set all the 3D coordinates of one chromosome.
       Num models comes from size of coords array if model not set"""
   
    particGroup = self._getParticleGroup(structure)
    allChromos = self.getChromosomes(structure)
    nPoints = len(particGroup[chromosome]['positions'])
    
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
       
    coordsGroup = self._getCoordsGroup(structure)

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
      
      
  def setAllCoords(self, coords, chromosomes=None, structure=None):
    """Set all the 3D model coordinates of specified chromosomes.
       Coords must be in chromosome order"""

    particGroup = self._getParticleGroup(structure)
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    elif not isinstance(chromosomes, (tuple, list, ndarray)):
      chromosomes = [chromosomes,]
    
    
    chromoSizes = [len(particGroup[c]['positions']) for c in chromosomes]
    nPoints = sum(chromoSizes)
    
    if coords.ndim == 2:
      coords = array([coords,])
    
    m, n, dims = coords.shape

    if n != nPoints:
      msg = 'Model coordinates must be an array of numModels x %d' % (nPoints,)
      raise(Exception(msg))
    
    coordsGroup = self._getCoordsGroup(structure)
        
    j = 0
    for i, chromo in enumerate(chromosomes):
      span = chromoSizes[i]
      chromoCoords = coords[:,j:j+span] # input models
      self._setData(chromo, coordsGroup, float, chromoCoords)
      
      j += span
    
         
  def setModelCoords(self, coords, model, chromosomes=None, structure=None):
    """Set all the 3D coordinates for one model of selected chromosomes.
       Corods must be in chromosome order"""
    
    particGroup = self._getParticleGroup(structure)

    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    elif not isinstance(chromosomes, (tuple, list, ndarray)):
      chromosomes = [chromosomes,]
    
    else:
      chromosomes = sorted(chromosomes)     
    
    if coords.ndim != 2:
      msg = 'Model coordinates must be an array of rank 2'
      raise(Exception(msg))
      
    
    chromoSizes = [len(particGroup[c]['positions']) for c in chromosomes]
    nPoints = sum(chromoSizes)
     
    n, dims = coords.shape
    
    if n != nPoints:
      msg = 'Model coordinates must be an array of length %d' % (nPoints,)
      raise(Exception(msg))
    
    coordsGroup = self._getCoordsGroup(structure) 
    
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

    
  def addChromosomes(self, positionsDict, backboneDict=None, interpolateCoords=True, structure=None):
    """Add new chromosome/DNA region identities"""
    
    
    from cUtil.apiUtil import interpolateChromoModelCoords

    coordsGroup = self._getCoordsGroup(structure)
    particGroup = self._getParticleGroup(structure)
    models = list(range(self.getNumModels(structure)))
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
      
      chromoPartGroup = self._getGroup(chromo, particGroup)
      coordArray = None
      
      if interpolateCoords and ('positions' in chromoPartGroup):
        prevPosDict[chromo] = array(chromoPartGroup['positions'], int)
       
      self._setData('backbone', chromoPartGroup, uint32, backbone)
      self._setData('positions', chromoPartGroup, float32, positions)
      self._setData('mapability', chromoPartGroup, float32, mapability)
      
      chromoGroup = self._getGroup(chromo, self.chromosomes)
      chromoGroup.attrs['limits'] = (positions[0], positions[-1])
      chromoGroup.attrs['region_of_interest'] = (positions[0], positions[-1], 0)
     
    nModels = self.getNumModels(structure)
    
    for chromo in prevPosDict:
      pDictA = {chromo: array(positionsDict[chromo], int32)}
      pDictB = {chromo: array(prevPosDict[chromo], int32)} 
      chromoCoords = []
      
      for model in range(nModels):
        coords = array(coordsGroup[chromo][model])
        coords = interpolateChromoModelCoords(pDictA, pDictB, coords)
        chromoCoords.append(coords)
      
      if chromoCoords:
        coords = array(chromoCoords)
        self.setModelChromosomeCoords(coords, chromo, structure=structure)
    
    names = set(self.chromosomes)
    names.update(positionsDict.keys())
    names = list(names)
    names.sort(key=naturalKey())
    self._chromosomes = names
        
    #if not nModels:
    #  self.setSpiralCoords(structure=structure)
      
    
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
    
    for structure in self.structures:
      strucGroup = self.structures[structure]
       
      if 'coords' in strucGroup:
        coordGroup = strucGroup['coords']
 
        for chromo in names:
          if chromo in coordGroup:
            self._delete(chromo, coordGroup)
    
      restraints = strucGroup['restraints']
      
      chromos = [x for x in restraints]
      for chromoA in chromos:
        if chromoA in names:
          self._delete(chromoA, restraints)
 
      for chromoA in restraints:
        chromos = [x for x in restraints[chromoA]]
 
        for chromoB in chromos:
          if chromoB in names:
            self._delete(chromoB, restraints[chromoA])
      
    for mainGroup in (self.workContacts, self.origContacts):
      for groupName in mainGroup:
        group = mainGroup[groupName]
        chrsA = [x for x in group]
        
        for chrA in chrsA:
          if chrA in names:
            self._delete(chrA, group)
          
          chrsB = [x for x in group.get(chrA, [])]
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
    
  
  def setRestraints(self, chromosomes=None, groupName=None, bboneSep=int(25e3),
                    binned=True, scale=1.0, exponent=-0.33, lower=0.8, upper=1.2,
                    minCount=1, maxPopDist=5.0, model=None, structure=None, ambigChromos=None):
    """Setup backbone and contact restraints accoring to specified parameters"""
    
    
    from cUtil.apiUtil import calcRestraints
    if not groupName:
      groupName = self.getDefaultContactGroup()
    
    if not structure:
      strucGroup = self.structures['0']
    else:
      strucGroup = self.structures[str(structure)]
    
    if chromosomes:
      chromosomes = sorted(chromosomes)
    else:
      chromosomes = self.getChromosomes(structure)
    
    ambigChromos = []
    contactDict  = self.getCachedContacts(groupName)
    chromosomes = [c for c in chromosomes if c in contactDict]
    
    if groupName in self.origContacts:
      isSingleCell = self.origContacts[groupName].attrs['isSingleCell']    
    else:
      isSingleCell = self.workContacts[groupName].attrs['isSingleCell']

    restrDict, posDict, ambigDict, bboneDict, binLstDict = calcRestraints(chromosomes, ambigChromos, contactDict, isSingleCell,
                                                                          [bboneSep, bboneSep], {}, 1.0, # scale not used
                                                                          exponent, lower, upper,
                                                                          minCount, maxPopDist, None)  
      
    if 'restraints' in strucGroup:
      del strucGroup['restraints']

    restGroup = self._getGroup('restraints', strucGroup)

    for chrA in restrDict:
      rGroup = self._getGroup(chrA, restGroup)
    
      for chrB in restrDict[chrA]:
        dataArray = restrDict[chrA][chrB]
        self._setData(chrB, rGroup, float32, dataArray)

    self.addChromosomes(posDict, bboneDict) 
 
  
  def setBackboneSpacing(self, bbSep=int(5e4), offset=None, structure=None):
    """Set the chromosomal sequence locations of backbone model particles"""
    
    chromos = self.getChromosomes(structure)
    limits = [self.getChromosomeLimits(c) for c in chromos]
    
    if offset is not None:
      for i, lim in enumerate(limits):
        lim = list(lim)
        lim[0] = offset
        limits[i] = lim
    
    self.setChromosomes(chromos, limits, bbSep)
    

  def setRefMapability(self, chromosome, values, structure=None):
    """Set the sequence mapability for organism reference"""
    
    particGroup = self._getParticleGroup(structure)
    dataArray = array(values, float)
    
    if dataArray.shape != particGroup['positions'].shape:
      s1 = ','.join(['%d' % x for x in dataArray.shape])
      s2 = ','.join(['%d' % x for x in particGroup['positions'].shape])
      raise Exception('Mapability data position size mismatch %s != %s' % (s1, s2))
    
    self._setData('mapability', particGroup, float, dataArray)
  
  
  def setMapability(self, chromosomes=None, structure=None):
    """Set the sequence mapability for particle positions"""
    
    if not self.genomeRef:
      return
    
    if not chromosomes:
      chromosomes = self.getChromosomes(structure)
    
    particGroup = self._getParticleGroup(structure)
    
    for chromosome in particGroup:
      if chromosome not in chromosomes:
        continue
      
      # Depends on genome build reference of file specified
      posMaps = self.genomeRef.getMapability(chromosome)
      posMaps = posMaps.tolist().reverse()
 
      chromoGroup = particGroup[chromosome]
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
  
  
  def setGlobalTransform(self, matrix, structure=None):
    """Set affine transformation associated with all models"""
       
    tformGroup = self._getTransformGroup(structure)
    
    self._setData('global', tformGroup, float, array(matrix))
    
     
  def setModelTransform(self, model, matrix, structure=None):
    """Set affine transformation associated with structure.
       This doesn't change the stored coordinates.
       This is for display and comparison only."""
    
    tformGroup = self._getTransformGroup(structure)
    
    if 'models' in tformGroup:
      transforms = array(tformGroup['models'])
      
      while len(transforms) <= model:
        append(transforms, eye(3), axis=0)
    
    else:
      numModels = self.getNumModels(structure)
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
  
  
  def getDataTrackMaxValue(self, source, code, recalc=False):
    """Fetch the maximum original and normalised data track values
       using cached values where possible"""
    

    dataGroup = self._getDataTrackGroup(source)
    
    if code in dataGroup:
      refLayer = self.getRefDataTrackGroup(source, code)
      dataLayer = dataGroup[code]
      maxNorm = dataLayer.attrs['display'][7]
    
      if recalc or not maxNorm:
        maxVals = []
 
        for chromo in refLayer:
          maxVals.append( array(refLayer[chromo]['values']).max(axis=0) )
 
        maxVals = array(maxVals).max(axis=0)
        maxOrig, maxNorm = maxVals
 
        dispVals = list(dataLayer.attrs['display'])
        dispVals[7] = maxNorm
        dataLayer.attrs['display'] = dispVals
      
      return maxNorm
      
  
  def renameDataTrack(self, old_code, new_code, source=None):
  
     source_group, source = self._getDataTrackSource(old_code, source)
    
     if not source_group:
       return
     
     old_data_track = source_group[old_code]
     new_data_track = self._getGroup(new_code, source_group)
     
     new_data_track.attrs['stranded'] = old_data_track.attrs['stranded']
     new_data_track.attrs['options']  = old_data_track.attrs['options'] 
     new_data_track.attrs['display']  = old_data_track.attrs['display'] 
     new_data_track.attrs['details']  = old_data_track.attrs['details'] 
           
     for chromo in old_data_track:
       old_chromo_group = old_data_track[chromo]
       new_chromo_group = self._getGroup(chromo, new_data_track)
       
       regions = array(old_chromo_group['regions'])
       values = array(old_chromo_group['values'])
       
       self._setData('regions', new_chromo_group, uint32, regions)
       self._setData('values', new_chromo_group, float32, values)
     
       if 'annotations' in old_chromo_group:
         self._setData('annotations', new_chromo_group, VL_str, array(old_chromo_group['annotations']))
       
       if 'models' in old_chromo_group:
         self._setData('models', new_chromo_group, uint32, array(old_chromo_group['models']))
     
     
     self._delete(old_code, source_group)

    
  def setDataTrack(self, code, source, regionDict, valueDict, annoDict=None, stranded=None,
                   modelDict=None, color=(1.0, 0.0, 0.0), scale=1.0, threshold=0.0,
                   showText=True, shape=0, details=None):
    """Set data values for an existing data layer, or make a new one should none exist.
       Regions, values and annotations set using chromsome keyed dictionaries."""
    # TBD : nuc.newDataTrack()
    
    
    dataGroup = self._getDataTrackGroup(source)
    annoDict = annoDict or {}
    
    if color:
      r, g, b = Color(color).rgb
    else:
      r, g, b = (1.0, 0.0, 0.0)
     
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
    
    if 'details' in dataLayer.attrs:
      if details:
        dataLayer.attrs['details'] = details
      
    else:
      dataLayer.attrs['details'] = details or code
    
    max_val_orig = []
    max_val_norm = []
    min_val_orig = []
    min_val_norm = []

    for chromo in regionDict:
      valueArray  = valueDict[chromo]
      
      if not len(valueArray):
        continue
      
      if valueArray.ndim == 1:
        orig_max = norm_max = valueArray.max()
        orig_min = norm_min = valueArray.min()
      
      else:
        orig_max, norm_max = valueArray.max(axis=0)
        orig_min, norm_min = valueArray.min(axis=0)
      
      max_val_orig.append(orig_max)
      max_val_norm.append(norm_max)
      min_val_orig.append(orig_min)
      min_val_norm.append(norm_min)
    
    if max_val_orig:
      max_val_orig = max(max_val_orig)
      max_val_norm = max(max_val_norm)
      min_val_orig = min(min_val_orig)
      min_val_norm = min(min_val_norm)
    
    for chromo in regionDict:
      regionArray = regionDict[chromo]
      n = len(regionArray)
      
      valueArray  = valueDict[chromo]
      if len(valueArray) != n:
        msg = 'Data track size mismatch: values (%d) must match number of regions (%d)'
        raise Exception(msg % (len(valueArray), n))
      
      if valueArray.ndim == 1:
        valueArray = array([valueArray, valueArray]).T
      
      if (max_val_orig == 0) and (min_val_orig == 0):
        if (max_val_norm == 0) and (min_val_norm == 0):
          valueArray = ones(valueArray.shape)
        
        else:
          valueArray[:,0] = valueArray[:,1]

      if min_val_norm < 0:
        valueArray[:,1] -= min_val_norm
        min_val_norm = 0.0
      
      if max_val_norm > 1.0:
        valueArray[:,1] /= max_val_norm
        max_val_norm = 1.0
      
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
   
    #self._setDataTrackIndicesMb(source, code)      
  
  
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
    
    #self._setDataTrackIndicesMb(source, code, [chromo])
    
    
    return chromoGroup
    

  def addDataTrack(self, code, regionDict, valueDict, source=EXTERNAL, annoDict=None, stranded=None, details=None):
    """Add a new data layer, with values for superimposing om structral models"""
    
    group = self._getDataTrackGroup(source)
    
    if code in group:
      msg = 'Attempt to add data track "%s" which already exists.' % code
      raise(Exception(msg))
      
    else:
      self.setDataTrack(code, source, regionDict, valueDict, annoDict, stranded, details=details)
    

  def removeDataTrack(self, source, code):
    """Remove a data layer"""
  
    group = self._getDataTrackGroup(source)
  
    if code in group:
      self._delete(code, group)
  
  
  def setInteractions(self, code, pairDict, valueDict=None, annoDict=None,
                      color=(0.0, 1.0, 1.0), style=0, width=1.0, showText=True):
     
    intGroup = self.interactions
    annoDict = annoDict or {}
    
    if color:
      r, g, b = Color(color).rgb
    else:
      r, g, b = (0.0, 1.0, 1.0)
     
    if code in intGroup:
      outerGroup = intGroup[code]
      
    else:
      outerGroup = self._getGroup(code, intGroup)
      showText = 1 if showText else 0
      
      outerGroup.attrs['options'] = (0, showText, int('0b1010101010101010', 2), 0)  # selected, showText, style, blank
      outerGroup.attrs['display'] = (r, g, b, 1.0, width, 0.0, 1.0) # r, g, b, a, width, minval, maxVal

    for chromoPair in pairDict:
      chromoA, chromoB = chromoPair
        
      pairArray = pairDict[chromoPair]
      n = len(pairArray)
      
      if not n:
        continue
      
      groupA = self._getGroup(chromoA, outerGroup)
      groupB = self._getGroup(chromoB, groupA)
      
      if pairArray.shape[1] != 4:
        msg = 'Interactions regions array must be of shape (n_interactions,4)'
        raise(Exception(msg))
      
      if valueDict:
        valueArray  = valueDict[chromoPair]
        if len(valueArray) != n:
          msg = 'Interactions data size mismatch: values (%d) must match number of pairs (%d)'
          raise Exception(msg % (len(valueArray), n))
 
        if valueArray.ndim == 1:
          valueArray = array([valueArray, valueArray]).T
      
      else:
        valueArray = ones((n, 2), float)         
 
      self._setData('regions', groupB, uint32, pairArray)
      self._setData('values', groupB, float32, valueArray)
 
      if annoDict and chromoA in annoDict:
        annoArray = annoDict[chromoPair]
 
        if len(annoArray):
          if len(annoArray) != n:
            msg = 'Data track size mismatch: annotations (%d) must match number of pairs (%d)'
            raise Exception(msg % (len(annoArray), n))
 
          self._setData('annotations', groupB, VL_str, annoArray)
 
  
  def addInteractions(self, code, regionDict, valueDict=None, annoDict=None):
  
    group = self.interactions
    
    if code in group:
      msg = 'Attempt to add interaction data "%s" which already exists.' % code
      raise(Exception(msg))
      
    else:
      self.setInteractions(code, regionDict, valueDict, annoDict)


  def removeInteractions(self, code): 
    
    group = self.interactions
    
    if code in group:
      self._delete(code, group)
  
  
  def _getImageGroup(self, code):  
    
    if code in self.images:
      return self.images[code]
    
    else:
      imgGroup = self._getGroup(code, self.images)
      self._setAttr(imgGroup, 'origin', zeros(3, float))
      self._setAttr(imgGroup, 'gridSize',ones(3, float))
    
      return imgGroup
  
  
  def getGridImages(self):
  
    codes = [c for c in self.images if 'grid' in self.images[c]]
    
    return codes
    
  
  def getCoordImages(self):
  
    codes = [c for c in self.images if 'coords' in self.images[c]]
    
    return codes
    
  
  def getImageGrid(self, code):
    """Get a full volumetric grid"""

    if code in self.images:
      imageGroup = self.images[code]
      
      if 'grid' in imageGroup:
        gridData = array(imageGroup['grid']) # Values[i,j,k]
        origin = imageGroup.attrs['origin']
        gridSize = imageGroup.attrs['gridSize']
      
        return origin, gridSize, gridData


  def setImageGrid(self, code, data, origin=(0.0,0.0,0.0), gridSize=(1.0,1.0,1.0)):
    """Set values for an existing full volumetric grid data"""
    
    data = array(data)
    if data.ndim != 3:
      raise Exception('Image voxel data must be three dimensional')
    
    imgGroup = self._getImageGroup(code)
    gridData = self._setData('grid', imgGroup, float, data, 'gzip')
    
    self._setAttr(imgGroup, 'origin', array(origin, float))
    self._setAttr(imgGroup, 'gridSize', array(gridSize, float))


  def addImageGrid(self, code, data, origin=(0.0,0.0,0.0), gridSize=(1.0,1.0,1.0)):
    """Add a new full volumetric grid data"""
    
    if code in self.images:
      raise Exception('Image with code "%s" already exists' % (code))
    
    self.setImageGrid(code, data, origin, gridSize)
    

  def removeImageGrid(self, code):
    """Remove a full volumetric grid data"""
    
    if (code in self.images) and ('grid' in self.images[code]):
      self._delete('grid', self.images[code])
      
      if 'coords' not in self.images[code]:
        self._delete(code, self.images)


  def getImageCoords(self, code, intensities=None, shapes=False):
    """
    Get coordinates of volumetric signal intensities (generally for sparse data).
    Optionally include fitted (not voxel grid) intensities.
    intensities is int or list of ints to specify intensity type
    """

    if code not in self.images:
      return
    
    imageGroup = self.images[code]  

    if 'coords' not in imageGroup:
      return
      
    data = array(imageGroup['coords'], float)
    
    # Add any intensities from a list of indices
    if intensities:
      if isinstance(intensities, int):
        intensities = [intensities]
        
      if 'intensities' in imageGroup:
        intGroup = imageGroup['intensities']
        n, m = intGroup.shape
 
        if n:
          idx = (i for i in intensities if i < m)
          idata = array(intGroup)[:,idx]
          data = concatenate([data, idata], axis=1)
      
      else: # Default intensities 1.0
        idata = ones([len(data), len(intensities)], float)
        data = concatenate([data, idata], axis=1)
      
    if shapes:
      if 'shapes' in imageGroup:
        shapes = array(imageGroup['shapes'], float)
      else: # Default shape (1.0, 1.0, 1.0)
        shapes = ones([len(data), 3], float)
    
      data = concatenate([data, shapes], axis=1)
  
  
    return data
        
  
  def getImageLabelCoords(self, code):
  
    if code not in self.images:
      return
    
    imageGroup = self.images[code]  

    if 'coords' not in imageGroup:
      return
      
    coords = array(imageGroup['coords'], float)
 
    if 'labels' in imageGroup:
      labels = list(imageGroup['labels'])
    else:
      labels = [''] * len(coords)
    
    coords, labels  
  
  
  def centreImageCoords(self):
  
    all_coords = None
    coord_cache = {}
    
    for code in self.images:
      data = array(self.images[code]['coords'])
      coord_cache[code] = data
      
      if all_coords is None:
        all_coords = array(data)
      else:
        all_coords = concatenate([all_coords, data], axis=0)
    
    if all_coords is not None:
 
      centre = all_coords.mean(axis=0)
 
      for code in self.images:
        imageGroup = self._getImageGroup(code)
        coords = coord_cache[code] - centre
        self._setData('coords', imageGroup, float, coords, 'gzip')
     
    
  def setImageCoords(self, code, coords, intensities=None, labels=None, shapes=None):

    imageGroup = self._getImageGroup(code)
    self._setData('coords', imageGroup, float, array(coords), 'gzip')
    
    n = len(coords)
    
    if intensities:
      if len(intensities) != n:
        msg = 'Number of intensities (%d) must match number of coords (%d)'
        raise Exception(msg % (len(intensities), n))
    
      self._setData('intensities', imageGroup, float, array(intensities), 'gzip')
    
    elif ('intensities' in imageGroup) and (len(imageGroup['intensities']) != n):
      self._delete('intensities', imageGroup)
    
    if labels:
      if len() != n:
        msg = 'Number of labels (%d) must match number of coords (%d)'
        raise Exception(msg % (len(labels), n))
    
      self._setData('labels', imageGroup, VL_str, labels,)

    elif ('labels' in imageGroup) and (len(imageGroup['labels']) != n):
      self._delete('labels', imageGroup)

    if shapes:
      if len() != n:
        msg = 'Number of shapes (%d) must match number of coords (%d)'
        raise Exception(msg % (len(shapes), n))
    
      self._setData('shapes', imageGroup, float, array(shapes), 'gzip')

    elif ('shapes' in imageGroup) and (len(imageGroup['shapes']) != n):
      self._delete('shapes', imageGroup)
    
    
  def addImageCoords(self, code, coords, intensities=None, labels=None, shapes=None):
   
    if code in self.images:
      raise Exception('Image with code "%s" already exists' % (code))
    
    self.setImageCoords(code, coords, intensities, labels, shapes)
    
    
  def removeImageCoords(self, code):
 
    if (code in self.images) and ('coords' in self.images[code]):
      self._delete('coords', self.images[code])
      
      if 'grid' not in self.images[code]:
        self._delete(code, self.images)
  
  
  def makeObsVsExpContactMap(self, group_name, default_bin_size=int(5e3), new_group_name=None):
    
    # Cis only at present
    
    bin_size = 2 * self.getContactsBinSize(group_name)
    
    if bin_size < 2:
      bin_size = default_bin_size
    
    if not new_group_name:
      new_group_name = group_name + '_ObsVsExp'
    
    contact_cache_dict = self.getCachedContacts(group_name)
    
    chromosomes = sorted(contact_cache_dict.keys())
    
    for chromo in chromosomes:
 
      if chromo not in contact_cache_dict[chromo]:
        continue
      
      start, end = self.getChromosomeLimits(chromo)
      
      obs = self.getContactMatrix(chromo, chromo, bin_size, group_name).astype(float)
       
      obs -= diag(diag(obs))
      obs /= obs.sum()
      n = len(obs)
      
      counts = zeros(n)
      sig = zeros(n)
      
      for i in range(n):
        for j in range(i,n):
          d = j-i
          sig[d] += obs[i,j]
          counts[d] += 1.0
      
      for c, j in enumerate(counts):
        if c:
          sig[j] /= c
      
      sig /= sig.sum()  
      
      exp = zeros((n,n), float)
      for i in range(n):
        exp[i,:i+1] = sig[:i+1][::-1]
        exp[i,i:] = sig[:n-i]
      
      vals = obs.sum(axis=0)
      vals /= vals.sum()
      
      exp *= outer(vals, vals)
      
      exp -= diag(diag(exp))
      exp /= exp.sum()     
      
      idx = (exp * obs).nonzero()
     
      matrix = obs
      matrix[idx] /= exp[idx]
      matrix = clip(matrix, 0.2, 5.0)
      #matrix[idx] = log(matrix[idx])
      #matrix /= 1000 * matrix.max()
    
      self.setContactMatrix(new_group_name, chromo, chromo, matrix, bin_size,
                            startA=start, startB=start, isSingleCell=False)
  
   
  def plotContactCovariance(self, groupName, chromosomes, binSize=int(5e5)):
    """Calculate covariance matrix for a specified set of
       intrachromosomal contacts
    """
    
    #from matplotlib import pyplot
    from util.Image import pixmapToImage   
    from numpy import outer
    
    cacheDict = self.getCachedContacts(groupName)
    
    for chromo in chromosomes:
    
      if chromo not in cacheDict:
        continue
 
      if chromo not in cacheDict[chromo]:
        continue
      
      obs = self.getContactMatrix(chromo, chromo, binSize,groupName).astype(float)
       
      obs -= diag(diag(obs))
      obs /= obs.sum()
      n = len(obs)
      
      counts = zeros(n)
      sig = zeros(n)
      
      for i in range(n):
        for j in range(i,n):
          d = j-i
          sig[d] += obs[i,j]
          counts[d] += 1.0
      
      for c, j in enumerate(counts):
        if c:
          sig[j] /= c
      
      sig /= sig.sum()  
      
      #pyplot.plot(sig)
      #pyplot.show()
      
      exp = zeros((n,n), float)
      for i in range(n):
        exp[i,:i+1] = sig[:i+1][::-1]
        exp[i,i:] = sig[:n-i]
      
      vals = obs.sum(axis=0)
      vals /= vals.sum()
      
      exp *= outer(vals, vals)
      
      exp -= diag(diag(exp))
      exp /= exp.sum()     
      
      idx = (exp * obs).nonzero()
     
      h = obs
      h[idx] /= exp[idx]
      h[idx] = log(h[idx])
      
      h = cov(clip(h, -2.0, 2.0))
    
      #z = ((exp * obs) == 0.0).nonzero()
      #h[z] = 0.0
      
      pos = (h > 0.0).nonzero()
      neg = (h < 0.0).nonzero()

      r = ones(obs.shape, float)
      g = ones(obs.shape, float)
      
      r[pos] = 0.0
      g[neg] = 0.0
      
      b = g
                   
      pixmap  = dstack([r,g,b])
      pixmap -= pixmap.min()
      pixmap /= pixmap.max()
      pixmap **= 0.25
      
      img = pixmapToImage(255*pixmap, 'RGB')
      img.show()

   
  def calcContactDataTrack(self, groupName=None, cis=True, trans=True, binSize=int(1e6),
                           trackName=None, includeEmpty=False, countInBin=False): 
    """Create a data track layer from contact information,
       per chromosome
    """
    
    from cUtil import dataLayer
    from collections import defaultdict
    
    if not groupName:
      groupName = self.getDefaultContactGroup()
      
    contDict = self.getContacts(groupName, None, cis, trans)
    regionDict = {}
    valueDict = {}
   
    for key in contDict:
      chrA, chrB = key
      valueDict[chrA] = []
      regionDict[chrA] = []
      valueDict[chrB] = []
      regionDict[chrB] = []
    
    
    pos_dict = defaultdict(list)
     
    for key in contDict:
      chrA, chrB = key
      
      data = contDict[key]
      
      if chrA == chrB:
        pos = data[:2].ravel()
        if not len(pos):
          continue
                  
        if binSize > 1:
          if not countInBin:
            # same-bin double counts
            bins_1 = array(data[0]/binSize, int)
            bins_2 = array(data[1]/binSize, int)
             
            idx = (bins_1 != bins_2).nonzero()
            pos = pos[idx]
          
          if len(pos):
            pos_dict[chrA].append(pos)
        
        else:
          pos_dict[chrA].append(pos)
        
      else:
        posA = data[0]
        posB = data[1]
        if not len(posA):
          continue
        
        pos_dict[chrA].append(posA)
        pos_dict[chrB].append(posB)
    
    for chrA in pos_dict:
      pos = concatenate(pos_dict[chrA], axis=0)
      regions = array([pos, pos+1], int32).T
      values = ones(len(regions), float)
       
      if binSize > 1:
        hist = dataLayer.regionBinValues(regions, values, int32(binSize))
        bins = arange(0.0, binSize*(1+len(hist)), binSize)
 
        regionDict[chrA] = vstack([bins[:-1], bins[1:]-1]).T
        valueDict[chrA] = vstack([hist, hist/(hist.max() or 1.0)]).T
      
      else:
        regionDict[chrA] = regions
        valueDict[chrA] = vstack([values, values])
         
      if not includeEmpty:
        valid = valueDict[chrA][:,0].nonzero()
        regionDict[chrA] = regionDict[chrA][valid]
        valueDict[chrA]  = valueDict[chrA][valid]
 
     
    if cis == trans:
      if not trackName:
        trackName = groupName
      color = (1.0, 1.0, 0.0)
    elif cis:
      if not trackName:
        trackName = groupName + '[cis]'
      color = (1.0, 0.0, 1.0)
    else:
      if not trackName:
        trackName = groupName + '[trans]'
      color = (0.0, 1.0, 1.0)
        
    return self.setDataTrack(trackName, DERIVED, regionDict, valueDict, color=color)
    
   
  def calcTransContactGscore(self, groupName=None, binSize=int(10e6)):
    """Calculate the G-test score (c.f. relative entropy) for
       single cell trans contacts."""
    
    # Cold get exp counts from cis/all contacts within bin
    
    if not groupName:
      groupName = self.getDefaultContactGroup()
    
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
      contacts = contactDictTrans[pairKey].T
      
      for posA, posB, count in contacts:
        binA = int(posA/binSize)
        binB = int(posB/binSize)
        
        key = frozenset([(chrA, binA), (chrB, binB)])
        obsCounts[key] += 1
    
    nCont = float(sum(obsCounts.values()))
    pExp = 1.0/float(len(obsCounts) or 1.0) # isotropic
    
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


class Color:
  """
  Class to store a single colour. Data is stored as floating point RGBA
  but may be extracted in various useful forms.
  """

  def __init__(self, *args):
   
    self._values = array([0.0,0.0,0.0,1.0])
    
    if args:
      if len(args) == 1:
        val = args[0]
        
        if type(val) is type(''):
          if len(val) == 9:
            self.setRgbaHex(val)
          else:
            self.setRgbHex(val)
        
        elif isinstance(val[0], (int, integer)):
          if len(val) == 4:
            self.setRgba32(*val)
          else:
            self.setRgb24(*val)
            
        elif len(val) == 4:
          self.setRgba(*val)
          
        else:
          self.setRgb(*val)
        
      elif len(args) == 3:
        if isinstance(args[0], (int, integer)):
          self.setRgb24(*val)
        else:
          self.setRgb(*val)
        
        
      elif len(args) == 4:
        if isinstance(args[0], (int, integer)):
          self.setRgba(*args)
        else:
          self.setRgba32(*args)
      
      else:
        raise Exception('Colour specification not understood')
  
  
  def _check8bitInt(self, *args):
  
    for v in args:
      if not isinstance(v, (int, integer)):
        valStr = ', '.join([str(v) for v in args])
        raise Exception('Value in color (%s) not of integer type' % valStr)
      
      if not (0 <= v <= 255):
        valStr = ', '.join([str(v) for v in args])
        raise Exception('Int value in color (%s) not in the range [0 .. 255]' % valStr)
      
      
  def _check01Float(self, *args):
  
    for v in args:
      if not isinstance(v, (float, floating)):
        valStr = ', '.join([str(v) for v in args])
        raise Exception('Value in color (%s) not of integer type' % valStr)
      
      if not (0.0 <= v <= 1.0):
        valStr = ', '.join([str(v) for v in args])
        raise Exception('Float value in color (%s) not in the range [0.0 .. 1.0]' % valStr)
      
      
  def _checkHexString(self, val, alpha):
    
    if val[0] != '#':
      raise Exception('Color hex string specification "%s" does not start with "#"' % val)

    if alpha:
      if len(val) != 9:
        raise Exception('Color hex string specification "%s" not length 9 "#RRGGBBAA"' % val)
      
      elif len(val) != 7:
        raise Exception('Color hex string specification "%s" not length 7 "#RRGGBB"' % val)
      
    for v in val[1:]:
      if v not in hexdigits:
        raise Exception('Color hex string specification "%s" must only use characters "%s"' % (val, hexdigits))

  
  def getRgb(self):
    """
    Get red, green, blue color ndarray as floats in interval [0.0, 1.0]
    """
    
    return tuple(self._values[:3])
  
  
  def setRgb(self, r, g, b):
    """
    Set red, green, blue with floats in interval [0.0, 1.0]
    """
    
    self._check01Float(r,g,b)    
    self._values[:3] = [r, g, b]
  
  
  def getHsv(self):
    """
    Get hue, saturation, value color ndarray as floats in interval [0.0, 1.0]
    """
    
    r, g, b = self._values[:3]
    return tuple(colorsys.rgb_to_hsv(r, g, b))
    
  
  def setHsv(self, h, s, v):
    """
    Set hue, saturation, value with floats in interval [0.0, 1.0]
    """
  
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    self.setRgb(r, g, b)
    
  
  def getRgba(self):
    """
    Get red, green, blue, alpha color ndarray as floats in interval [0.0, 1.0]
    """
  
    return tuple(self._values)
  
  
  def setRgba(self, r, g, b, a):
    """
    Set red, gree, blue, alpha with floats in interval [0.0, 1.0]
    """

    self._check01Float(r,g,b, a)    
    self._values = array([r, g, b, a])
  
  
  def getRgb24(self):
    """
    Get red, green, blue color ndarray as ints in interval [0, 255]
    """
  
    return tuple(array(self._values[:3]*255, int))
  
  
  def setRgb24(self, r, g, b):
    """
    Set red, green, blue with ints in interval [0, 255]
    """
    
    self._check8bitInt(r, g, b)
    self.setRgb(r/255.0, g/255.0, b/255.0)
  
  
  def getRgba32(self):
    """
    Get red, green, blue, alpha color ndarray as ints in interval [0, 255]
    """
  
    return tuple(array(self._values*255, int))
  
  
  def setRgba32(self, r, g, b, a):
    """
    Set red, green, blue, alpha with ints in interval [0, 255]
    """
  
    self._check8bitInt(r, g, b, a)
    self.setRgba(r/255.0, g/255.0, b/255.0, a/255.0)
    
  
  def getRgbaHex(self):
    """
    Get red, green, blue, alpha 32-bit hexadecimal string representation #RRGGBBAA
    """
  
    return '#%02X%02X%02X%02X' % self.getRgba32()
    
    
  def setRgbaHex(self, hexStr):
    """
    Set red, green, blue, alpha with 32-bit hexadecimal string representation #RRGGBBAA
    """

    self._checkHexString(hexStr, True)
    
    r = int(hexStr[1:3], 16)
    g = int(hexStr[3:5], 16)
    b = int(hexStr[5:7], 16)
    a = int(hexStr[7:9], 16)
    
    self.setRgba32(r, g, b, a)
  
  
  def getRgbHex(self):
    """
    Get red, green, blue 24-bit hexadecimal string representation #RRGGBB
    """
   
    return '#%02X%02X%02X' % self.getRgb24()
  
  
  def setRgbHex(self, hexStr):
    """
    Set red, green, blue with 24-bit hexadecimal string representation #RRGGBB
    """
    
    self._checkHexString(hexStr, False)
    
    r = int(hexStr[1:3], 16)
    g = int(hexStr[3:5], 16)
    b = int(hexStr[5:7], 16)
    
    self.setRgb24(r, g, b)
    
  
  def inverseGrey(self):
    """
    Returns a new Color object which is an inverse monochrome
    that provides good contrast, e.g. for text on a background
    """

    r, g, b, a = self.rgba32
 
    m = (11*r + 16*g + 5*b)/32
 
    if (m > 192) or (m < 64):
      m = 255-m
    elif m<128:
      m += 128
    elif m<192:
      m -= 128
  
    return Color((m,m,m,a))

  rgb = property(getRgb, setRgb)
  hsv = property(getHsv, setHsv)
  rgba = property(getRgba, setRgba)
  rgb24 = property(getRgb24, setRgb24)
  rgba32 = property(getRgba32, setRgba32)
  rgbHex = property(getRgbHex, setRgbHex)
  rgbaHex  = property(getRgbaHex, setRgbaHex)
 

class Chromosome:
  
  @argTypeCheck(nuc=Nucleus, name=(str, isWord), limits=isArrayN2)
  def __init__(self, nuc, name, limits=None):
    
    if name in nuc._chromoDict:
      raise NucApiException('Chromosome name "%s" already in use' % name)
    
    self.nuc = nuc
    self.name = name
    
    if not limits:
      limits = (0,1)
    
    nuc.setChromosomes([name], [limits])
  
  
  def getLimits(self):
  
    return self.nuc.getChromosomeLimits(self.name)
  
  
  def getIsShown(self):
  
    return self.nuc.getChromoDisplayed(self.name)
  
  
  def setIsShown(self, isShown):
  
    self.nuc.setChromoDisplayParams(self.name, isShown)
    
  
  def getColor(self):
  
    return self.nuc.getChromoColor(self.name)
  
  
  def setColor(self, color):
  
    self.nuc.setChromoColor(self.name, color)
  
  
  # colorMode = property
  # displayMode = property
  # useLabels = property
  
  # def delete()
    
  color = property(getColor, setColor)
  limits = property(getLimits)
  
class Structure:

  def __init__(self):
    pass

class ContactGroup:

  def __init__(self):
    pass

class DataTrack:

  def __init__(self):
    pass

class Image:

  def __init__(self):
    pass


if __name__ == '__main__':
  
  pass
  
  # TBD import and run testing
  
  
