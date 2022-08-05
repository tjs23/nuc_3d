from libc.math cimport log, abs
import sys, cython, gc
from numpy cimport ndarray, uint8_t
from numpy import zeros, ones, uint8, int32, vstack

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def readPairedSam(fileName, fileMode):

  cdef int a, b, c, d, e, i, j, n
  cdef int nBytes
  cdef samfile_t *samFile
  cdef bam_header_t *samHead
  cdef int buf = 5000
  
  samFile = samopen(fileName, fileMode, NULL)
  
  if samFile == NULL:
    msg = 'Could not open BAM/SAM file "%s"' % fileName
    raise IOError(msg)
  
  alignBuff = <bam1_t*>calloc(1, sizeof(bam1_t))
  samHead = samFile.header
  n = samHead.n_targets
  
  cdef ndarray[int, ndim=2] lims = zeros((n,2), int32)
  cdef ndarray[int, ndim=2] last = zeros((n,n), int32)
  cdef ndarray[int, ndim=2] sizes = ones((n,n), int32)
  cdef ndarray[int, ndim=2] contacts
  chrNames = []
  
  for c in range(n):
    chrNames.append(str(samHead.target_name[c]))
    lims[c,0] = samHead.target_len[c]
    
  matrices = {}
  for a in range(n):
    chrA = chrNames[a]
    matrices[chrA] = {}
  
    for b in range(a, n):
      chrB = chrNames[b]
      matrices[chrA][chrB] = zeros((1,3), int32)
      
  nBytes = samread(samFile, alignBuff)
  
  while nBytes > 0:
    c = alignBuff.core.tid
    d = alignBuff.core.mtid
    
    if c == d:
      a = c
      b = d
      i = <int>alignBuff.core.pos
      j = <int>alignBuff.core.mpos
      
      if i > j:
        e = i
        i = j
        j = e
    
    elif c > d:
      a = d
      b = c
      i = <int>alignBuff.core.mpos
      j = <int>alignBuff.core.pos
    
    else:
      a = c
      b = d
      i = <int>alignBuff.core.pos
      j = <int>alignBuff.core.mpos
      
    chrA = chrNames[a]
    chrB = chrNames[b]
    
    contacts = matrices[chrA][chrB]
    e = last[a,b]
    
    if e > 0:
      c = contacts[e-1,0]
      d = contacts[e-1,1]
    
      if (c==i) and (d==j):
        nBytes = samread(samFile, alignBuff)
        continue
    
    if e >= sizes[a,b]: 
      contacts = vstack([contacts, zeros((buf,3), int32)])
      matrices[chrA][chrB] = contacts
      sizes[a,b] += buf
    
    if i > lims[a,1]:
      lims[a,1] = i
    
    if i < lims[a,0]:
      lims[a,0] = i

    if j > lims[b,1]:
      lims[b,1] = j
    
    if j < lims[a,0]:
      lims[b,0] = j
    
    contacts[e,0] = i
    contacts[e,1] = j
    contacts[e,2] = 1
    
    last[a,b] += 1
    
    nBytes = samread(samFile, alignBuff)
    
  posDict = {}
  for a in range(n):
    chrA = chrNames[a]
  
    for b in range(a, n):
      chrB = chrNames[b]
      
      if last[a,b] > 0:
        matrices[chrA][chrB] = matrices[chrA][chrB][:last[a,b]]
        posDict[chrA] = lims[a]
        posDict[chrB] = lims[b]
      else:
        del matrices[chrA][chrB]
      
    if not matrices[chrA]:
      del matrices[chrA]
    
  samclose(samFile)
  
  return matrices, posDict
  
  
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def readBinPairedSam(fileName, fileMode, int binSize=50000):

  cdef int a, b, c, d, i, j, n
  cdef int nBytes, nNonZero
  cdef samfile_t *samFile
  cdef bam_header_t *samHead
  cdef char *chrName
  cdef int chrLen, maxLen=0
  cdef uint64_t fStart
  
  samFile = samopen(fileName, fileMode, NULL)
  
  if samFile == NULL:
    msg = 'Could not open BAM/SAM file "%s"' % fileName
    raise IOError(msg)
  
  fStart = bam_tell(samFile.x.bam)
  alignBuff = <bam1_t*>calloc(1, sizeof(bam1_t))
  samHead = samFile.header
  n = samHead.n_targets
  
  cdef ndarray[int, ndim=1] chrLens = zeros((n,), int32)
  chrNames = []
  
  for c in range(n):
    chrName = samHead.target_name[c]
    chrNames.append(str(chrName))
    chrLen = 1 + <int>samHead.target_len[c]/binSize
    chrLens[c] = chrLen
    if chrLen > maxLen:
      maxLen = chrLen
    
  cdef ndarray[uint8_t, ndim=2] matrix = zeros((maxLen,maxLen), uint8)
    
  matrices = {}
    
  for c in range(n):
    chrA = chrNames[c]
    matrices[chrA] = {}
  
    for d in range(c, n):
      chrB = chrNames[d]
      
      bam_seek(samFile.x.bam, fStart, 0)
 
      nNonZero = 0
      nBytes = samread(samFile, alignBuff)
      while nBytes > 0:
        a = alignBuff.core.tid
        b = alignBuff.core.mtid
 
        if (a == c) and (b == d):
          i = <int>alignBuff.core.pos/binSize
          j = <int>alignBuff.core.mpos/binSize
          if matrix[i,j] == 0:
            nNonZero += 1
            
          matrix[i,j] += 1
        
        elif (a == d) and (b == c):
          i = <int>alignBuff.core.mpos/binSize
          j = <int>alignBuff.core.pos/binSize
          if matrix[i,j] == 0:
            nNonZero += 1
            
          matrix[i,j] += 1
        
        nBytes = samread(samFile, alignBuff)
      
      if nNonZero > 0:
        contacts = zeros((nNonZero, 3), uint8)
        print chrA, chrB, nNonZero
 
        a = 0
        for i in range(chrLens[c]):
          for j in range(chrLens[d]):
            if matrix[i,j] > 0:
              contacts[a,0] = i
              contacts[a,1] = j
              contacts[a,2] = matrix[i,j]
 
              a += 1
              matrix[i,j] = 0
 
        matrices[chrA][chrB] = contacts

  samclose(samFile)
  
  return matrices
  
  
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def pairedSamSeqSepHist(fileName, fileMode):

  cdef int a, b, i, j, k, n
  cdef int nBytes
  cdef samfile_t *samFile
  cdef bam_header_t *samHead
  cdef int buf = 100000
  cdef double sep
  
  samFile = samopen(fileName, fileMode, NULL)
  
  if samFile == NULL:
    msg = 'Could not open BAM/SAM file "%s"' % fileName
    raise IOError(msg)
  
  alignBuff = <bam1_t*>calloc(1, sizeof(bam1_t))
  samHead = samFile.header
  n = samHead.n_targets
  
  cdef ndarray[int, ndim=2] hist = zeros((n,100), int32)
  chrNames = []
  
  for c in range(n):
    chrNames.append(str(samHead.target_name[c]))   
      
  nBytes = samread(samFile, alignBuff)
  
  while nBytes > 0:
    a = alignBuff.core.tid
    b = alignBuff.core.mtid
    
    if a != b:
      nBytes = samread(samFile, alignBuff)
      continue
    
    i = <int>alignBuff.core.pos
    j = <int>alignBuff.core.mpos
    
    sep = <double>abs(i-j)
    
    if sep == 0.0:
      nBytes = samread(samFile, alignBuff)
      continue
    
    k = <int>(10.0*log(sep)/2.302585092994046)
   
    #k = <int>(sep/1e4)
    if k < 100:
      hist[a,k] += 1   
    
    nBytes = samread(samFile, alignBuff)
    
  chrHist = {}
  for a in range(n):
    chrHist[chrNames[a]] = hist[a]
    
  samclose(samFile)
  
  return chrHist
  
  
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def pairedSamReadBasic(fileName, fileMode):

  cdef int a, b, i, j, k, n
  cdef int nBytes
  cdef samfile_t *samFile
  cdef bam_header_t *samHead
  cdef double sep
  
  samFile = samopen(fileName, fileMode, NULL)
  
  if samFile == NULL:
    msg = 'Could not open BAM/SAM file "%s"' % fileName
    raise IOError(msg)
  
  alignBuff = <bam1_t*>calloc(1, sizeof(bam1_t))
  samHead = samFile.header
  n = samHead.n_targets
  
  chrNames = []
  data = []
  
  for c in range(n):
    chrNames.append(str(samHead.target_name[c]))   
      
  nBytes = samread(samFile, alignBuff)
  
  while nBytes > 0:
    a = alignBuff.core.tid
    b = alignBuff.core.mtid
    
    if a != b:
      nBytes = samread(samFile, alignBuff)
      continue
    
    i = <int>alignBuff.core.pos
    j = <int>alignBuff.core.mpos
    
    sep = <double>abs(i-j)
    
    if sep == 0.0:
      nBytes = samread(samFile, alignBuff)
      continue
    
    data.append((chrNames[a], chrNames[b], i, j)) 
    
    nBytes = samread(samFile, alignBuff)
     
  samclose(samFile)
  
  return data

