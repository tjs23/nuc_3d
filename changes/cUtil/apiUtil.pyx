from libc.math cimport abs, sqrt, ceil
from numpy cimport ndarray
import numpy as np
from numpy import ones, zeros, int32, float32, sort, empty, array, arange, concatenate
import cython

def interpolateChromoModelCoords(posDict, prevPosDict, ndarray[double, ndim=2] coords):
  """
  Interpolate x,y,z particle positions for an array of seq positions to a new  
  seq positions e.g. for a change in binned resolution.
  """
  
  cdef int i, j, i0, j0, p1, p2, n, m, d, dMin
  cdef double f, g
  cdef ndarray[int, ndim=1] positions
  cdef ndarray[int, ndim=1] prevPositions
  
  cdef int a, b
  
  # Get index offsets for each chromo
  
  i = 0
  offsets = {}
  for chrA in sorted(posDict):
    offsets[chrA] = i
    i += len(posDict[chrA])
  
  a = i
  cdef ndarray[double, ndim=2] newCoords = empty((i, 3), float)
  
  i = 0
  prevOffsets = {}
  for chrA in sorted(prevPosDict):
    prevOffsets[chrA] = i
    i += len(prevPosDict[chrA])  
  
  b = i
  
  for chrA in posDict:
    i0 = offsets[chrA]
    j0 = prevOffsets[chrA]
    positions = posDict[chrA]
    prevPositions = prevPosDict[chrA]
    n = len(positions)
    m = len(prevPositions)
    
    for i in range(n):

      #find closest old positions for coordinate interpolation
      p1 = 0
      dMin = positions[i]-prevPositions[0]
      
      for j in range(1,m):
        d = positions[i]-prevPositions[j]
        
        if abs(d) < abs(dMin):
          p1 = j
          dMin = d #closest pos
          
        elif abs(d) > abs(dMin): # Seq positions were in order
          break  
    
      if dMin == 0: #new position coincides with an old position
        p2 = p1
      elif dMin > 0: #new pos is above p1
        if (p1+1 > m-1):
          print 'warning: InterpolatedCoords: query pos ',positions[i],' is not within the position limits [',prevPositions[0],',',prevPositions[m-1],'], mapping onto the last bead.'
        p2 = min(p1+1, m-1)
      else: #new pos is below p1
        if (p1-1 < 0):
          print 'warning: InterpolatedCoords: query pos ',positions[i],' is not within the position limits [',prevPositions[0],',',prevPositions[m-1],'], mapping onto the first bead.'
        p2 = p1
        p1 = max(0, p1-1)
        dMin = positions[i] - prevPositions[p1]
      #print 'new resolution pos ',positions[i],' is in the interval [a,b] = [',prevPositions[p1],',',prevPositions[p2],'] with a distance from a ',dMin

      #p1 <= newpos <= p2, dMin = pos[newpos] - pos[p1]
      #calculate coordinates
      if p1 == p2:
        newCoords[i0+i, 0] = coords[j0+p1, 0]
        newCoords[i0+i, 1] = coords[j0+p1, 1]
        newCoords[i0+i, 2] = coords[j0+p1, 2]
 
      else: #interpolate
        f = <float>dMin/<float>(prevPositions[p2]-prevPositions[p1])
        g = 1.0 - f
        
        newCoords[i0+i, 0] = g * coords[j0+p1, 0] + f * coords[j0+p2, 0]
        newCoords[i0+i, 1] = g * coords[j0+p1, 1] + f * coords[j0+p2, 1]
        newCoords[i0+i, 2] = g * coords[j0+p1, 2] + f * coords[j0+p2, 2]

      #print "Interpolated coordinate ",newCoords[i0+i,0], " is between ", coords[j0+p1,0], " and ", coords[j0+p2,0]
      #print "Interpolated coordinate ",newCoords[i0+i,1], " is between ", coords[j0+p1,1], " and ", coords[j0+p2,1]
      #print "Interpolated coordinate ",newCoords[i0+i,2], " is between ", coords[j0+p1,2], " and ", coords[j0+p2,2]
  
  return newCoords
  
  

def concatenateRestraints(restraintDict, posDict, seqScale, backboneLower=0.1, backboneUpper=1.1):
  """
  Joins restraints stored in a dict according to chromo pairs into long concatenated arrays.
  Indices of restraints relate to concatenated chromo seq pos.
  Add-in all the backbone restraints for sequential particles.
  """
  
  cdef int i, n, nRest, m
  cdef int startA, startB
  cdef int bbl = int32(backboneLower)
  cdef int bbu = int32(backboneUpper)
  cdef double v, sScale = seqScale
  
  # Get total max number restraints and final restraint index offset for all chromos
    
  offsets = {}
  nRest = 0
  i = 0
  
  for chrA in sorted(posDict):
    n = len(posDict[chrA])
    offsets[chrA] = i
    nRest += n-1 # One fewer because at end of chain
    i += n
        
  for chrA in restraintDict:
    for chrB in restraintDict[chrA]:
      nRest += restraintDict[chrA][chrB].shape[1]
 
  cdef ndarray[int, ndim=1] positions
  cdef ndarray[double, ndim=2] restraints
  cdef ndarray[int, ndim=2] indicesArray   = empty((nRest, 2), int32)
  cdef ndarray[double, ndim=2] boundsArray = empty((nRest, 2), float)
  cdef ndarray[int, ndim=1] ambigArray     = empty(nRest, int32)
  
  # Add backbone path restraints   
  
  m = 0
  for chrA in posDict:
    
    positions = posDict[chrA]
    n = len(positions)
    startA = offsets[chrA]
    
    for i in range(n-1):
      
      v = (positions[i+1]-positions[i])/sScale
      
      indicesArray[m,0] = i + startA
      indicesArray[m,1] = i + startA + 1
      
      boundsArray[m,0] = v * bbl # lower
      boundsArray[m,1] = v * bbu # upper
      
      ambigArray[m] = 1
      
      m += 1
      
  # Add regular restraints for chromo pairs 
  
  for chrA in restraintDict:
    startA = offsets[chrA]
   
    for chrB in restraintDict[chrA]:
      startB = offsets[chrB]
      
      restraints = restraintDict[chrA][chrB]
      n = restraints.shape[1]
      
      for i in range(n):
        indicesArray[m,0] = <int>restraints[0,i] + startA
        indicesArray[m,1] = <int>restraints[1,i] + startB
 
        boundsArray[m,0] = restraints[4,i] # lower
        boundsArray[m,1] = restraints[5,i] # upper
 
        ambigArray[m] = 1
        
        m += 1

 
  return indicesArray, boundsArray, ambigArray


def fillInDomLst(ndarray[int, ndim=2] domLst,
                 ndarray[int, ndim=1] binLimits,
                 ndarray[int, ndim=1] limits):
  """
  NEED TO DO PER CHRM
  Take an empty or partially filled domLst, as well as the binLimits:
  gaps smaller than first val are ignored, second val is max spacer size
  before sub division. Returns a filled domLst.
  """

  cdef ndarray[int, ndim=2] filledDomLst
  cdef ndarray[int, ndim=2] gap

  if array(domLst).shape == (1,0):
    filledDomLst = fillGap(limits[1], limits[0], binLimits)
    
  else:
    gap = fillGap(domLst[0,0], limits[0], binLimits)
    initDom = array([domLst[0]])
    filledDomLst = (utilConcat([gap, initDom]))
    
    for i in range(1, len(domLst)):
      if i == (len(domLst) - 1):
        endGap = array(fillGap(limits[1], domLst[i, 1] + 1, binLimits))
        filledDomLst = utilConcat([filledDomLst,
                                   fillGap((domLst[i,0]), (domLst[i-1,1]), binLimits) , array([domLst[i]]),
                                   endGap])
      else:
        filledDomLst = utilConcat([filledDomLst,
                                   fillGap((domLst[i,0])-1, (domLst[i-1,1]) + 1, binLimits),
                                   array([domLst[i]]) ])

  return filledDomLst


def utilConcat(concList):
  
  newLst = [x for x in concList if x.shape != (0,2)]
  
  if len(newLst) == 1:
    return newLst[0]
  else:
    return concatenate(newLst)
    

def fillGap(int upper, int lower, ndarray[int, ndim=1] binLimits):
  """
  Tweak this as necessary. 
  bins defined as [start, end]
  """
  cdef ndarray [int, ndim=1] current
  cdef ndarray [int, ndim=2] outArray
  cdef int k
  
  current = zeros(2, int32)
  outArray = zeros(((((upper-lower)/binLimits[0]) + 5), 2), int32)

  current[0] = lower
  current[1] = lower+binLimits[1] - 1
  # current[0] = lower+binLimits[1] - 1
  # current[1] = lower+(2 * binLimits[1]) - 1
  k = 0

  while current[1] < upper:
    outArray[k] = current
    k += 1
    newLowLim = current[0] + binLimits[1]
    newUpLim  = current[1] + binLimits[1]
    
    if newUpLim < upper:
      current[0] = newLowLim
      current[1] = newUpLim
      
    else:
      """
      if new lowlim >= upper, return, else return newLowLim / upper. 
      """
      # current[1] = upper
      # if upper - current[1] <= binLimits[0]:
      #   return outArray[:k]
        
      # else:
      # newLowLim = current[1] + 1
        # current = [newLowLim, newUpLim]
      current[0] = newLowLim
      current[1] = upper

  return outArray[:k]

"""
def fillGapTest():
  # cdef ndarray [int, ndim=2] filledDomLst
  # filledDomLst = zeros((100, 2), int32)
  cdef ndarray [int, ndim=1] binLimits
  cdef int k
  binLimits = zeros(2, int32)
  binLimits[0] = 5
  binLimits[1] = 10
  # k = 0
  # filledDomLst = fillGap(100, 10, binLimits)
  filledDomLst = fillGap(100, 99, binLimits)
  print(filledDomLst)

def fillDomLstTest():
  cdef ndarray [int, ndim=2] domLst
  cdef ndarray [int, ndim=1] binLimits
  cdef ndarray [int, ndim=1] limits
  domLst = array([[10,30], [45,80], [81,85]], dtype=int32)
  # domLst = array([[]], dtype=int32)
  binLimits = array([3,10], dtype=int32)
  limits = array([10,100], dtype=int32)
  print(fillInDomLst(domLst, binLimits,limits))
"""

def bpToBin(int bp, ndarray[int, ndim=2] binLst, int n):
  for i in range(n):
    if binLst[i,0] <= bp <= binLst[i,1]:
      return i



def calcRestraints(chromosomes, contactDict, pyIsSingleCell,
                   pyBboneSep=[500000, 500000], domLstDict={},
                   float scale=1.0, float exponent=-0.33,
                   float lower=0.8, float upper=1.2,
                   pyMinCount=2, float maxPopDist=5.0,
                   pyModel=-1):
    
  """
  Liam's notes:
  bboneSep is now a list of two values, that represent the range of acceptable bin sizes.
  (e.g. ignore a gap smaller than the min, don't allow a bin larger than the max)
  It's acceptable values are None and [a, b] where b>=a.
  domLst is a list of domains. It can take the values None, the empty list, or a list of length 2 lists (start / end of the domain).
  Expect domLst to be well formed (e.g. [[a,b], [c,d]] where a<b, and b<c etc)
  The combination of pyBboneSep and domLst controls what the function does:
  None / None -> Native resolution
  [a,b] / [] -> Binned at the average of a and b.
  [a,b] / [..] -> Using a and b, attempt to fill in the gaps in the domLst (if necessary). If this fails, return some kind of error.
  
  
  Input dict of contact arrays keyed contactDict[chrA][chrB]
  
  Contact array shape (3:[bpA, bpB, numObs], nContacts)
  
  Calculates pairwise distance restraints for a pair of chromosomes
  given contact information and also extracts a list of particle seq
  positions which will be used to construct the coodinate models.
  
  Goes through all contacts for all chromosome pairs
  - Need to know the extent of the seq positions to know the particle ranges
  - Need only the unique seq positions
  - Particle seq pos made according to range and backbone options
  - For each contact the particle index is known
  - Restraints may be binned onto the particle
  
  Return distance restraint arrays, chromosomal positions, backbone spacer positions. 
  
  """  
  
  cdef int isSingleCell = 1 if pyIsSingleCell else 0
  cdef int binned = 1 if domLstDict else 0
  cdef int minCount = pyMinCount
  cdef int model = pyModel
     
  cdef int i, j, k, a, b, c, n, na, nb
  cdef int mFilter
  cdef double v
  cdef ndarray[uint, ndim=2] contacts      # Contact matrix (3:(posA, posB, nObs), nContacts)
  cdef ndarray[double, ndim=2] restraints  # Distance restraints (6:(), nRestraints)
  cdef ndarray[int, ndim=2] binMatrix      # Temp array for binned contacts
  cdef ndarray[int, ndim=2] indices        # Restraint indices for seq pos (nPos, 2:(idxA, idxB))
  cdef ndarray[int, ndim=1] counts 
  cdef ndarray[int, ndim=1] positionsA
  cdef ndarray[int, ndim=1] positionsB
  cdef ndarray[int, ndim=1] backboneA
  cdef ndarray[int, ndim=1] backboneB
  cdef ndarray[int, ndim=1] offsets
  cdef ndarray[int, ndim=2] limits # shape: (chromoId, 2:[start, end])
  cdef ndarray[int, ndim=1] bboneSep = array(pyBboneSep, int32)

  binLstDict = {}
  nContDict  = {}
  posDict    = {}
  bboneDict  = {}
  indexDict  = {}
  restraintDict = {}
  chromos = set(chromosomes)
  
  # Num contacts per chromo for mem allocations
  for chrA in contactDict:
    if chrA not in chromos:
      continue
    
    for chrB in contactDict[chrA]:
      if chrB not in chromos:
        continue
        
      contacts = contactDict[chrA][chrB]
      n = len(contacts[0])
      
      if chrA in nContDict:
        nContDict[chrA] += n
      else:
        nContDict[chrA] = n
      
      if chrB in nContDict:
        nContDict[chrB] += n
      else:
        nContDict[chrB] = n
  
  chrIdx  = {chrA:i for i, chrA in enumerate(chromosomes)}
  c = len(chromosomes)
  
  limits = zeros((c,2), int32)
  counts = zeros(c, int32)
  
  for chrA in chromos:
    posDict[chrA] = empty(nContDict[chrA], int32)
    if domLstDict is not None:
      if chrA in domLstDict:
        domLstDict[chrA] = array(domLstDict[chrA], int32)
      else:
        domLstDict[chrA] = array([[]], int32)
  
  # Get chromosome ranges and used seq positions
  for chrA in contactDict:
    if chrA not in chromos:
      continue

    a = chrIdx[chrA]
    indexDict[chrA] = {}
    
    for chrB in contactDict[chrA]:
      if chrB not in chromos:
        continue
        
      b = chrIdx[chrB]
      
      positionsA = posDict[chrA]
      positionsB = posDict[chrB]
      contacts = contactDict[chrA][chrB]
      n = len(contacts[0])
      
      indices = empty((n,2), int32)
      indexDict[chrA][chrB] = indices
      
      for i in range(n):
        
        if limits[a,0] == 0: # zero is not a valid seq pos anyhow
          limits[a,0] = contacts[0,i]
        elif contacts[0,i] < limits[a,0]:
          limits[a,0] = contacts[0,i]
        
        if limits[b,0] == 0:
          limits[b,0] = contacts[1,i]
        elif contacts[1,i] < limits[b,0]:
          limits[b,0] = contacts[1,i]
        
        if contacts[0,i] > limits[a,1]:
          limits[a,1] = contacts[0,i]
      
        if contacts[1,i] > limits[b,1]:
          limits[b,1] = contacts[1,i]
      
        positionsA[counts[a]] = contacts[0,i]
        positionsB[counts[b]] = contacts[1,i]
        
        indices[i,0] = counts[a]
        indices[i,1] = counts[b]
        
        counts[a] += 1
        counts[b] += 1
  
  if (domLstDict is None) and (bboneSep is not None):
    avSep = (bboneSep[0] + bboneSep[1]) / 2
    for a in range(c):
      limits[a,0] = avSep * (limits[a,0]/avSep)
      limits[a,1] = avSep * <int>(ceil(<float>limits[a,1]/avSep))

  elif (domLstDict is not None) and (bboneSep is not None):
    for a in range(c):
      limits[a,0] = bboneSep[0] * (limits[a,0]/bboneSep[0])
      limits[a,1] = bboneSep[1] * <int>(ceil(<float>limits[a,1]/bboneSep[1]))     

  offsets = zeros(c, int32)
  for a in range(1, c):
    offsets[a] = offsets[a-1] + limits[a-1,1] - limits[a-1,0]
  
  # remove duplicate positions, set backbone arrays
  for chrA in chromosomes:
    a = chrIdx[chrA]
  
    if domLstDict is not None:
      binLstDict[chrA] = fillInDomLst(domLstDict[chrA], bboneSep, limits[a])

    if (domLstDict is not None) and (bboneSep is not None): # Set regularly spaced seq positions only
      spacer = ((limits[a,1] + bboneSep[1]) - limits[a,0]) / (len(binLstDict[chrA]) - 1)
      posDict[chrA] = arange(limits[a,0], limits[a,1] + bboneSep[1], spacer, int32)
        
      bboneDict[chrA] = ones(len(posDict[chrA]), int32)
    
    else:
      positionsA = posDict[chrA]
      
      if (bboneSep is not None) and (domLstDict is None): # inject backbone spacers
        sep = (bboneSep[0] + bboneSep[1]) / 2
        positionsA = positionsA + arange(limits[a,0], limits[a,1]+sep, sep, int32)
      
      sort(positionsA)
      n = len(positionsA)
 
      if n ==0:
        continue
 
      positionsB = empty(n, int32)
      positionsB[0] = positionsA[0]
      
      if (bboneSep is not None) and (domLstDict is None):
        sep = (bboneSep[0] + bboneSep[1]) / 2
        backboneA = empty(n, int32)
        
        if positionsA[0] % sep == 0:
          backboneA[0] = 1
        else:
          backboneA[0] = 0
        
        j = 1
        for i in range(1,n):
          if positionsA[i] != positionsA[i-1]:
            positionsB[j] = positionsA[i]
            
            if positionsA[i] % sep == 0:
              backboneA[j] = 1
            else:
              backboneA[j] = 0
            
            j += 1
        
        posDict[chrA] = positionsB[:j]
        bboneDict[chrA] = backboneA[:j]
      
      else:
        j = 1
        for i in range(1,n):
          if positionsA[i] != positionsA[i-1]:
            positionsB[j] = positionsA[i]
            j += 1
 
        posDict[chrA] = positionsB[:j]
        bboneDict[chrA] = zeros(j, int32)
  
  # Get indices, do any binning of observations
  
  for chrA in contactDict:
    if chrA not in chromos:
      continue
      
    a = chrIdx[chrA]
    backboneA = bboneDict[chrA]
    restraintDict[chrA] = {}
    
    for chrB in contactDict[chrA]:
      if chrB not in chromos:
        continue
        
      b = chrIdx[chrB]
      backboneB = bboneDict[chrB]
      
      contacts = contactDict[chrA][chrB]
      n = len(contacts[0])  
      restraints = empty((6, n), float)
  
      if model >= 0 and len(contacts) == 4:
        mFilter = 1
      else:
        mFilter = 0
 
      if (domLstDict is not None) and (bboneSep is not None): 
        """
        1. Create the binDict from the domLstDict
        2. Create the bbonedict (e.g. everything is a bbone elem. (TODO: ASK TIM))
        3. Add the restraints based on the domLst.
        """
        #build contact matrix at this resolution: binMatrix
        na = (binLstDict[chrA].shape)[0]
        nb = (binLstDict[chrB].shape)[0]
        # na = 1 + (limits[a,1]-limits[a,0])/bboneSep 
        # nb = 1 + (limits[b,1]-limits[b,0])/bboneSep 
        binMatrix = zeros((na, nb), int32)
  
  return restraintDict, posDict, bboneDict, binLstDict
   
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
#For a (list of) query position(s), find the flanking positions in a sequence of positions,
#and use linear interpolation of the coordinates to get the coordinates of the query position(s).
def getInterpolatedCoords(ndarray[int, ndim=1] seqPos,
                          ndarray[double, ndim=2] coords,
                          ndarray[int, ndim=1] queryPos):
  
  cdef int i, j, k, a, b
  cdef int nPos = len(seqPos)
  cdef int nQuery = len(queryPos)
  cdef int d, dMin, dMax=seqPos.max()
  cdef double f
  cdef ndarray[double, ndim=2] p = zeros((nQuery,3))  
  
  if len(seqPos) != len(coords):
    print 'error: getInterpolatedCoords: sequence position and coordinate list lengths mismatch.'

  for k in range(nQuery):
    dMin = dMax
  
    for i in range(nPos):
      d = seqPos[i]-queryPos[k]
 
      if abs(d) < abs(dMin):
        a = i
        dMin = d #closest pos
      else:
        break  

    #find the other end of the interval such that a <= query <= b
    if dMin == 0: # query is a
      b = a
    elif dMin > 0: # query is below a
      if (a-1 < 0):
        print 'warning: getInterpolatedCoords: query pos ',queryPos[k],' is not within the position limits [',seqPos[0],',',seqPos[nPos-1],'], mapping onto the first bead.'
      b = a
      a = max(0, a-1) # use 1st if query is below 1st
      dMin = seqPos[b]-seqPos[a] - dMin
    else: # query is above a
      if (a+1 > nPos-1):
        print 'warning: getInterpolatedCoords: query pos ',queryPos[k],' is not within the position limits [',seqPos[0],',',seqPos[nPos-1],'], mapping onto the last bead.'
      b = min(a+1, nPos-1) #use last if query is beyond last
      dMin = - dMin
    #print 'query pos ',queryPos[k],' is in the interval [a,b] = [',seqPos[a],',',seqPos[b],'] with a distance from a ',dMin
 
    #a <= newpos <= b, dMin = pos[b] - pos[query]
    #calculate coordinates
    if a == b:
      for j in range(3):
        p[k,j] =  coords[a,j]

    else: #interpolate
      f = <float>dMin/<float>(seqPos[b]-seqPos[a])
 
      for j in range(3):
        p[k,j] = (1.0-f) * coords[a,j] + f * coords[b,j] 
  
  return p
  

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def getClosestPoints(ndarray[double, ndim=1] refCoord,
                     ndarray[double, ndim=2] coords,
                     int numPoints=10):
  
  cdef int i, j, k, n = len(coords)
  cdef float dx, dy, dz, d2, dlim = 0.0
  cdef float dlim2 = 0.0
  cdef float x, y, z
  cdef ndarray[int, ndim=1] closest = zeros(numPoints, int32)
  cdef ndarray[double, ndim=1] dists2 = zeros(numPoints, float)
  
  if n < numPoints:
    raise Exception('Number of points in getClosestPoints cannot be larger than number of coords')
  
  x = refCoord[0]
  y = refCoord[1]
  z = refCoord[2]
  
  for i in range(numPoints):
    dx = coords[i,0]-x
    dy = coords[i,1]-y
    dz = coords[i,2]-z
    d2 = dx*dx + dy*dy + dz*dz
    dists2[i] = d2
    closest[i] = i
    
    if d2 > dlim2:
      dlim2 = d2
      dlim = sqrt(d2)
      k = i
    
  for i in range(numPoints,n):
    dx = abs(coords[i,0]-x)
    if dx > dlim:
      continue
    
    dy = abs(coords[i,1]-y)
    if dy > dlim:
      continue
    
    dz = abs(coords[i,2]-z)
    if dz > dlim:
      continue
      
    d2 = dx*dx + dy*dy + dz*dz
    if d2 > dlim2:
      continue
    
    if d2 == 0.0:
      continue
    
    dists2[k] = d2
    closest[k] = i
    dlim2 = d2
    
    for j in range(numPoints):
      if dists2[j] > dlim2:
        dlim2 = dists2[j]
        dlim = sqrt(dlim2)
        k = j
      
  return closest


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def getIsolatedPairs(ndarray[int, ndim=2] positions,
                     int threshold=2000000, int posErr=100):
                  
  cdef int i, j, n = len(positions)
  cdef int pA, pB, pC, pD
  cdef ndarray[int, ndim=1] unsupported = ones(n, int32)
  
  for i in range(n):
    pA = positions[i,0]
    pB = positions[i,1]

    for j in range(n):
      if j == i:
        continue
      
      pC = positions[j,0]
      pD = positions[j,1]
      
      if (posErr < abs(pC-pA) < threshold) and (posErr < abs(pD-pB) < threshold):
        unsupported[i] = 0
        break
 
      elif (posErr < abs(pD-pA) < threshold) and (posErr < abs(pC-pB) < threshold):
        unsupported[i] = 0
        break
    
  indices = unsupported.nonzero()[0]
  
  return indices
  
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def binContacts(ndarray[int, ndim=2] contacts,
                ndarray[int, ndim=2] binMatrix,
                int offsetA, int offsetB,
                int binSize=2000000,
                int symm=0, int transpose=0):
  
  cdef int i, a, b
  cdef int n, m, nCont = len(contacts[0])
  
  n = len(binMatrix)
  m = len(binMatrix[0])
  
  if transpose:
    for i in range(nCont):
      b = (contacts[0,i]-offsetA)/binSize
      a = (contacts[1,i]-offsetB)/binSize
 
      if (0 <= a < n) and (0 <= b < m):
        binMatrix[a,b] += contacts[2,i]
 
        if symm and (a != b):
          binMatrix[b,a] += contacts[2,i]
  
  else:
    for i in range(nCont):
      a = (contacts[0,i]-offsetA)/binSize
      b = (contacts[1,i]-offsetB)/binSize
 
      if (0 <= a < n) and (0 <= b < m):
        binMatrix[a,b] += contacts[2,i]
 
        if symm and (a != b):
          binMatrix[b,a] += contacts[2,i]
  
  return binMatrix
  

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def intMatrixToSparse(ndarray[int, ndim=2] matrix,
                      int scaleA, int scaleB,
                      int offsetA=0, int offsetB=0):
 
  cdef int i, j, n, m, a, c=0
  
  n = len(matrix)
  m = len(matrix[0])
  
  for i in range(n):
    for j in range(m):
      if matrix[i,j] != 0:
        c += 1
        
  cdef ndarray[int, ndim=2] sparse = zeros((3, c), int32)
  
  c = 0
  for i in range(n):
    a = offsetA + (i * scaleA)
    
    for j in range(m):
      
      if matrix[i,j] != 0:
        sparse[0,c] = a
        sparse[1,c] = offsetB + (j * scaleB)
        sparse[2,c] = matrix[i,j]
        c += 1 
  
  return sparse
  
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def calcCoordDensity(ndarray[double, ndim=2] coords,
                     double rMax=10.0):
  
  cdef int i, j, k, l, a, b, c, n=len(coords)
  cdef int mx, my, mz, m=0
  cdef ndarray[double, ndim=1] densities = zeros(n, float)
  cdef double rBin, d, minX, minY, minZ, maxX, maxY, maxZ
  
  minX = coords[0,0]
  minY = coords[0,1]
  minZ = coords[0,2]
  maxX = coords[0,0]
  maxY = coords[0,1]
  maxZ = coords[0,2]
  
  for i in range(1, n):
    
    if coords[i,0] < minX:
      minX = coords[i,0]

    elif coords[i,0] > maxX:
      maxX = coords[i,0]
  
    if coords[i,1] < minY:
      minY = coords[i,1]

    elif coords[i,1] > maxY:
      maxY = coords[i,1]
      
    if coords[i,2] < minZ:
      minZ = coords[i,2]

    elif coords[i,2] > maxZ:
      maxZ = coords[i,2]
  
  rMax = min(rMax, max((maxX-minX), max((maxY-minY), (maxZ-minZ)))/10.0)
  
  rBin = rMax / 5.0
  mx = int(ceil((maxX-minX)/rBin))
  my = int(ceil((maxY-minY)/rBin))
  mz = int(ceil((maxZ-minZ)/rBin))
  d = 1.0/<float>(515*n)
  
  cdef ndarray[double, ndim=3] binDensity = zeros((mx, my, mz), float)
  
  for i in range(n):
    a = <int>((coords[i,0]-minX)/rBin)
    b = <int>((coords[i,1]-minY)/rBin)
    c = <int>((coords[i,2]-minZ)/rBin)
    
    for j in range(-5, 6):
      if a+j < 0:
        continue
      if a+j >= mx:
        continue
    
      for k in range(-5, 6):
        if b+k < 0:
          continue
        if b+k >= my:
          continue
          
        if (j*j) + (k*k) > 25:
          continue
          
        for l in range(-5, 6):
          if c+l < 0:
            continue
          if c+l >= mz:
            continue  
          if (j*j) + (k*k) + (l*l)> 25:
            continue
    
          binDensity[a+j,b+k,c+l] += d

  for i in range(n):
    a = <int>((coords[i,0]-minX)/rBin)
    b = <int>((coords[i,1]-minY)/rBin)
    c = <int>((coords[i,2]-minZ)/rBin)
    densities[i] = binDensity[a,b,c]

  
  return densities

    
