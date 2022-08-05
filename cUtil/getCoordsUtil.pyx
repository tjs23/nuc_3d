from libc.math cimport abs
import numpy as np
cimport numpy as np
import cython

DTYPE = np.int32
ctypedef np.int32_t DTYPE_t

FTYPE = np.float64
ctypedef np.float64_t FTYPE_t

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
#For a (list of) query position(s), find the flanking positions in a sequence of positions,
#and use linear interpolation of the coordinates to get the coordinates of the query position(s).
def getInterpolatedCoordsSort(np.ndarray[dtype=DTYPE_t, ndim=1] seqPos,
                              np.ndarray[dtype=FTYPE_t, ndim=2] coords,
                              np.ndarray[dtype=DTYPE_t, ndim=1] queryPos):

  cdef unsigned int i, j, k, a, b
  cdef int nPos = seqPos.shape[0]
  cdef int nQuery = queryPos.shape[0]
  cdef int d, dMin, dA, dB
  cdef int dMax = seqPos.max()
  cdef float f
  cdef np.ndarray[dtype=FTYPE_t, ndim=2] p = np.zeros((nQuery,3))
  cdef unsigned int imin, imax, imid

  if nPos != coords.shape[0]:
    print 'error: getInterpolatedCoords: sequence position and coordinate list lengths mismatch.'


  for k in range(nQuery):
    d = dMax
    imin = 0
    imax = nPos
    imid = (imax + imin) / 2
    while (imax > imin + 1):
      d = seqPos[imid]-queryPos[k]
      if (d==0):
        imin = imid
        imax = imid
        break
      elif (d < 0):
        imin = imid
      elif (d > 0):
        imax = imid
      imid = (imax + imin) / 2

    # necessary because of how we define the while loop. could probably make more elegant
    if ((imax - imin) == 2):
      imax = imin + 1

    a = imin
    b = imax

    dA = seqPos[a] - queryPos[k]
    dB = seqPos[b] - queryPos[k]

    dMin = min(abs(dA), abs(dB))

    # #find the other end of the interval such that a <= query <= b
    # if d == 0: # query is a
    #   pass
    # elif d > 0: # query is below a
    #   if (a-1 < 0):
    #     print 'warning: getInterpolatedCoords: query pos ',queryPos[k],' is not within the position limits [',seqPos[0],',',seqPos[nPos-1],'], mapping onto the first bead.'
    #     a = 0
    #   # b = a
    #   # a = max(0, a-1) # use 1st if query is below 1st
    # else: # query is above a
    #   if (b+1 > nPos-1):
    #     print 'warning: getInterpolatedCoords: query pos ',queryPos[k],' is not within the position limits [',seqPos[0],',',seqPos[nPos-1],'], mapping onto the last bead.'
    #   # b = min(b+1, nPos-1) #use last if query is beyond last
    # #print 'query pos ',queryPos[k],' is in the interval [a,b] = [',seqPos[a],',',seqPos[b],'] with a distance from a ',dMin

    # dOldMin = dMax
    # for i in range(nPos):
    #   dOld = seqPos[i]-queryPos[k]
      
    #   if abs(dOld) < abs(dOldMin):
    #     aOld = i
    #     dOldMin = dOld #closest pos
    #   else:
    #     break  
    # print(dOldMin)

    # #find the other end of the interval such that a <= query <= b
    # if dOldMin == 0: # query is a
    #   bOld = aOld
    # elif dOldMin > 0: # query is below a
    #   if (aOld-1 < 0):
    #     print 'warning: getInterpolatedCoords: query pos ',queryPos[k],' is not within the position limits [',seqPos[0],',',seqPos[nPos-1],'], mapping onto the first bead.'
    #   bOld = aOld
    #   aOld = max(0, aOld-1) # use 1st if query is below 1st
    #   dOldMin = seqPos[bOld]-seqPos[aOld] - dMin
    # else: # query is above a
    #   if (aOld+1 > nPos-1):
    #     print 'warning: getInterpolatedCoords: query pos ',queryPos[k],' is not within the position limits [',seqPos[0],',',seqPos[nPos-1],'], mapping onto the last bead.'
    #   bOld = min(aOld+1, nPos-1) #use last if query is beyond last
    #   dOldMin = - dOldMin

    # if (dMin != dOldMin):
    #   print("a: {}, b: {}, dMin: {}".format(a, b, dMin))
    #   print("aOld: {}, bOld: {}, dOldMin: {}".format(aOld, bOld, dOldMin))


    #calculate coordinates
    if a == b:
      for j in range(3):
        p[k,j] =  coords[a,j]

    else: #interpolate
      f = <float>dMin/<float>(seqPos[b]-seqPos[a])
 
      for j in range(3):
        p[k,j] = (1.0-f) * coords[a,j] + f * coords[b,j] 
  
  return p
