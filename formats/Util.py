import os, gzip
from numpy import int32, float64, float32, zeros, arange, array
from math import ceil
from cUtil.apiUtil import binContacts, intMatrixToSparse
from cUtil.dataLayer import regionBinValues
from h5py import File
from numpy import load
from NucApi import Nucleus

STRUCTURE_FORMATS = {'HDF5':['*.hdf5','*.hdf','*.h5'],
                     'JSON':['*.json'],
                     'NDArray':['*.npz'],
                     'PDB':['*.pdb'],
                     'XLSX':['*.xlsx'],
                     'TSV':['*.tsv','*.n3d'],
                     'Nuc3D':['*.nuc']}

CONTACT_FORMATS = {'NCC':['*.ncc','*.ncc.gz'],
                   'SAM':['*.bam','*.sam','*.sam.gz'],
                   'PFE':['*.pfe','*.fend_pairs', '*.txt'],
                   'HDF5':['*.hdf5','*.hdf','*.h5'],
                   'JSON':['*.json'],
                   'NDArray':['*.npz'],
                   'TSV':['*.tsv','*.mat'],
                   'Nuc3D':['*.nuc']}
                   
DATA_TRACK_FORMATS = {'BED':['*.bed'],
                      'WIG':['*.wig'],
                      'bedGraph':['*.bedgraph'],
                      'broadPeak':['*.broadPeak'],
                      'TSV':['*.tsv'],
                      'HDF5':['*.hdf5','*.hdf','*.h5'],
                      'JSON':['*.json'],
                      'NDArray':['*.npz'],
                      'Nuc3D':['*.nuc']}

INTERACTIONS_FORMATS = {'HDF5':['*.hdf5','*.hdf','*.h5'],
                       'JSON':['*.json'],
                       'NDArray':['*.npz'],
                       'TSV':['*.tsv'],
                       'Nuc3D':['*.nuc']}

CONTACTS = 'contacts'
STRUCTURE = 'structure'
DATA_TRACK = 'dataTrack'
INTERACTIONS = 'interactions'

def addCompressedExts(format, exts):
  
  if format in ('HDF5', 'NDArray', 'SAM', 'XLSX', 'Nuc3D'):
    return exts
  
  exts2 = []
  
  for ext in exts:
    exts2.append(ext)
    exts2.append(ext + '.gz')
  
  return exts2


for fmtDict in (STRUCTURE_FORMATS, CONTACT_FORMATS, DATA_TRACK_FORMATS, INTERACTIONS_FORMATS):
  for f in fmtDict:
    fmtDict[f] = addCompressedExts(f, fmtDict[f])
 
  
def splitExtension(filePath):

  if filePath.endswith('.gz'):
    a, b = os.path.splitext(filePath[:-3])
    return a, b+'.gz'
    
  else:
    return os.path.splitext(filePath)
  
  
def getFileObj(filePath, mode='r'):

  if filePath.endswith('.gz'):
    return gzip.open(filePath, mode)
  
  else:
    return open(filePath, mode)


def getFileType(filePath):
  
  if filePath.endswith('.gz'):
    filePath = filePath[:-3]
  
  fileRoot, fileExt = os.path.splitext(filePath)
  
  if fileExt in ('.hdf5', '.hdf', '.h5', '.nuc', '.temp'):
    try:
      nuc = Nucleus(filePath, 'r')
      return 'nuc', 'nuc'
      
    except Exception:
      try:
        root = File(filePath, 'r')
        format = 'HDF5'
 
        if 'chromosomes' in root and 'id' in root.attrs:
          return 'nuc', 'nuc'
 
        elif DATA_TRACK in root:
          return DATA_TRACK, format
 
        elif STRUCTURE in root:
          return STRUCTURE, format

        elif CONTACTS in root or 'contactMatrix' in root:
          return CONTACTS, format

        elif INTERACTIONS in root:
          return INTERACTIONS, format
 
      except Exception:
        pass

  
  elif fileExt in ('.npz'):
    try:
      dataDict = load(filePath, 'r')
      format = 'NDArray'
      
      for key in dataDict:
        if key.startswith('dataTrack/'):
          return DATA_TRACK, format
        
        elif key.startswith('coords/'):
          return STRUCTURE, format

        elif key.startswith('contacts/'):
          return CONTACTS, format

        elif key.startswith('interactions/'):
          return INTERACTIONS, format
          
        elif key.startswith('contactMatrix/'):
          return CONTACTS, format
        
    except Exception:
      pass

  elif fileExt in ('.json'):
    
    fileObj = open(filePath)
    data = fileObj.read(1024)
    fileObj.close()
    format = 'JSON'
    
    if DATA_TRACK in data:
      return DATA_TRACK, format
    
    elif STRUCTURE in data:
      return STRUCTURE, format

    elif INTERACTIONS in data:
      return INTERACTIONS, format

    elif CONTACTS in data or 'contactMatrix' in data:
      return CONTACTS, format

  elif fileExt in ('.pfe','.fend_pairs','.bam','.sam'):
    
    if fileExt in ('.bam','.sam'):
      format = 'BAM'
    
    else:
      format = 'PFE'
    
    return CONTACTS, format
  
  elif fileExt in ('.bed', '.wig', '.bedgraph', '.broadPeak'):
    patt = '*' + fileExt
    
    for format, fileExts in DATA_TRACK_FORMATS:
      if patt in fileExts:
        return DATA_TRACK, format
    
  return None, None
      

def binSparseContacts(data, binSize):

  posA = data[:,0]
  posB = data[:,1]
  
  sa = posA.min()
  sb = posB.min()
  ea = posA.max()
  eb = posB.max()
  
  n = ceil((ea-sa)/float(binSize))
  if n <= 0:
    return
 
  m = ceil((eb-sb)/float(binSize))
  if m <=0:
    return
 
  binMatrix = zeros((n,m), int32)
  binContacts(data, binMatrix, int32(sa), int32(sb), int32(binSize))
  data = intMatrixToSparse(binMatrix, binSize, binSize, int32(sa), int32(sb)).T
  
  return data


def binTrackData(regions, values, binSize):

  binSize = int32(binSize)
  regions = array(regions, int32)
  
  histA = regionBinValues(regions, values[:,0].astype(float64), binSize)
  histB = regionBinValues(regions, values[:,1].astype(float64), binSize)
  
  starts = binSize * arange(0, len(histA))
  
  idx = histA.nonzero()
  histA = histA[idx]
  histB = histB[idx]
  starts = starts[idx]

  ends = (binSize-1) + starts
  
  regions = array([starts, ends], int32).T
  values = array([histA, histB], float32).T
  
  
  return regions, values
