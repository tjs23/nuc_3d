from __future__ import division
import numpy as np
cimport numpy as np
import cython

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def getPercOL(np.ndarray[np.int32_t, ndim=2] coordsA,
              np.ndarray[np.int32_t, ndim=2] coordsB):
    """
    (gets the number of bps in A that overlap with B, divides by total bps in A)
    Both coords arrays should be for a specific chromosome from the data track.
    Go through each coord in coordsA, find it's overlap with coordsB.
    """

    cdef np.int32_t aS
    cdef np.int32_t aE

    cdef ssize_t aCoordsLen = coordsA.shape[0]
    cdef ssize_t i 
    cdef np.int32_t totalBP = 0

    cdef np.int32_t olBP = 0

    for i in range(aCoordsLen):
        aS = coordsA[i,0]
        aE = coordsA[i,1]
        totalBP += aE - aS + 1
        olBP += getBPOL(aS, aE, coordsB)

    return(olBP / totalBP)

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
cdef getBPOL(np.int32_t aS,
             np.int32_t aE,
             np.ndarray[np.int32_t, ndim=2] coordsB):
    """
    For a specific bp range, get the number of bps that overlap with the provided array.
    """
        
    cdef np.int32_t bS 
    cdef np.int32_t bE 

    cdef ssize_t bCoordsLen = coordsB.shape[0]
    cdef ssize_t i 

    for i in range(bCoordsLen):

        bS = coordsB[i,0]
        bE = coordsB[i,1]

        if (aS == bE):
            return(1)

        # Does A intersect B? 
        if (aE >= bS) and (bE >= aS):
            # A is 'in' B.
            # A:   ---
            # B: -------
            if (aS >= bS) and (aE <= bE):
                return(aE - aS + 1)
            # A overlaps B, and sticks out on the rhs.
            # A:    -----
            # B: ----- 
            elif (aS >= bS) and (aE >= bE):
                if aE > bE:
                    # This might be wrong..
                    return((bE - aS + 1) + getBPOL(bE+1, aE, coordsB))
                else:
                    return(bE - aS + 1)
            # A overlaps B, and sticks out on the lhs.
            # A: ----
            # B:   -----
            elif (aS <= bS) and (aE <= bE):
                return (aE - bS + 1) 
            # B falls within A.
            # A: -------
            # B:   ---
            elif (aS <= bS) and (aE >= bE):
                if aE > bE:
                    return((bE - bS + 1) + getBPOL(bE+1, aE, coordsB))
                else:
                    return(bE - bS + 1)

    # No overlap found.
    return(0)

