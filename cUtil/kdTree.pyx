from libc.math cimport abs, sqrt, ceil
from numpy cimport ndarray
from numpy import empty, array, int32, zeros
import cython, sys

cdef double BIG = sys.float_info.max

cdef void arraySwapInt1d(ndarray[int, ndim=1] vals, int a, int b):
  
  cdef int c
  c = vals[a]
  vals[a] = vals[b]
  vals[b] = c
  

cdef double dist2(int ndim,
                  ndarray[double, ndim=1] v1,
                  ndarray[double, ndim=1] v2):
     
  cdef double d, d2 = 0.0
  
  for i in range(ndim):
    d = v1[i] - v2[i]
    d2 += d*d
        
  return d2
  
  
cdef double dist2box(int ndim, 
                     ndarray[double, ndim=2] region,
                     ndarray[double, ndim=1] v):
    
  cdef double d, d2 = 0.0
  
  for i in range(ndim):
    if v[i] < region[0,i]:
      d = v[i] - region[0,i]
      d2 += d*d
      
    if v[i] > region[1,i]:
      d = v[i] - region[1,i]
      d2 += d*d
  
  return d2
    

cdef selecti(int k, ndarray[int, ndim=1] indx,
            int n, ndarray[double, ndim=1] arr):

  cdef int i, ia, ir, j, l, mid
  cdef double a
  
  l = 0
  ir = n-1
  
  while 1:
    if ir <= l+1:
      if (ir == l+1) and (arr[indx[ir]] < arr[indx[l]]):
        arraySwapInt1d(indx, l, ir)
        
      return indx[k]
        
    else:
      mid = (l+ir) >> 1 # /2
      arraySwapInt1d(indx, mid, l+1)
      
      if arr[indx[l]] > arr[indx[ir]]:
        arraySwapInt1d(indx, l, ir)

      if arr[indx[l+1]] > arr[indx[ir]]:
        arraySwapInt1d(indx, l+1, ir)

      if arr[indx[l]] > arr[indx[l+1]]:
        arraySwapInt1d(indx, l, l+1)
        
      i = l+1
      j = ir
      ia = indx[l+1]
      a = arr[ia]
      
      while 1:
        i += 1
        while arr[indx[i]] < a:
          i += 1          
          
        j -= 1
        while arr[indx[j]] > a:
          j -= 1
          
        if j < i:
          break
          
        arraySwapInt1d(indx, i, j) 
      
      indx[l+1] = indx[j]
      indx[j] = ia
      
      if j >= k:
        ir = j-1
        
      if j <= k:
        l = i  


class KdTree(object):

  def __init__(self, points):
    
    self.points = array(points, float)
    self.ndim = len(points[0])
    
    ptindx, rptindx, regions, parent, child1, child2, ptlo, pthi = _construct(self.points)
    
    self.ptindx = ptindx
    self.rptindx = rptindx
    self.regions = regions
    self.parent  = parent
    self.child1 = child1
    self.child2 = child2
    self.ptlo = ptlo
    self.pthi = pthi
  
 
  def getNearest(self, ndarray[double, ndim=2] coords):

    cdef i, j, n = len(coords) 
    cdef ndarray[double, ndim=2] nearest = empty((n, 3), float)
    cdef ndarray[int, ndim=1] indices = empty(n, int32)
    
    cdef ndarray[double, ndim=2] points = self.points
    cdef ndarray[int, ndim=1] ptindx = self.ptindx
    cdef ndarray[int, ndim=1] rptindx = self.rptindx
    cdef ndarray[double, ndim=3] regions = self.regions
    cdef ndarray[int, ndim=1] parent = self.parent
    cdef ndarray[int, ndim=1] child1 = self.child1
    cdef ndarray[int, ndim=1] child2 = self.child2
    cdef ndarray[int, ndim=1] ptlo = self.ptlo
    cdef ndarray[int, ndim=1] pthi = self.pthi
    
    indices = _nearest(coords, points, ptindx, regions,
                       child1, child2, ptlo, pthi)
    
    for i in range(n):
      nearest[i] = points[indices[i]]
      
    return nearest


cdef _construct(ndarray[double, ndim=2] points):      
  
  cdef int npts = len(points)
  cdef ndim = len(points[0])
  cdef ndarray[int, ndim=1] ptindx = empty(npts, int32)
  cdef ndarray[int, ndim=1] rptindx = empty(npts, int32)
  
  cdef int ntmp, m, k, kk, j, nowtask, jbox, tmom, tdim, ptlo, pthi
  cdef ndarray[int, ndim=1] hp     
  cdef ndarray[double, ndim=1] cp     

  cdef ndarray[int, ndim=1] taskmom = empty(50, int32)    
  cdef ndarray[int, ndim=1] taskdim = empty(50, int32)    
  
  for k in range(0, npts):
    ptindx[k] = k
    
  m = 1
  
  ntmp = npts
  while ntmp:
    m <<= 1    # *= 2
    ntmp >>= 1 # /= 2
  
  cdef int nboxes = 2*npts - (m >> 1)
  if m < nboxes:
    nboxes = m
  
  nboxes -= 1
  
  cdef ndarray[double, ndim=3] regions = zeros((nboxes, 2, 3), float)
  cdef ndarray[int, ndim=1] parent  = zeros(nboxes, int32)
  cdef ndarray[int, ndim=1] child1 = zeros(nboxes, int32)
  cdef ndarray[int, ndim=1] child2 = zeros(nboxes, int32)
  cdef ndarray[int, ndim=1] boxptlo = zeros(nboxes, int32)
  cdef ndarray[int, ndim=1] boxpthi = zeros(nboxes, int32)
  
  cdef ndarray[double, ndim=1] coords = empty(ndim*npts, float)
  
  j = 0
  kk = 0
  while j < ndim:
    for k in range(npts):    
      coords[kk+k] = points[k,j]
    
    j += 1
    kk += npts
  
  cdef ndarray[double, ndim=2] region = empty((2, ndim), float)
  
  for k in range(ndim):
    region[0,k] = -BIG
    region[1,k] = BIG
    regions[0,0,k] = -BIG
    regions[0,1,k] = BIG   
    
  parent[0] = 0
  child1[0] = 0
  child2[0] = 0
  boxptlo[0] = 0
  boxpthi[0] = npts-1
  
  jbox = 0
  taskmom[1] = 0
  taskdim[1] = 0
  nowtask = 1
  
  while nowtask:
    tmom = taskmom[nowtask]
    tdim = taskdim[nowtask]
    nowtask -= 1
    ptlo = boxptlo[tmom]
    pthi = boxpthi[tmom]
    hp = ptindx[ptlo:]
    cp = coords[tdim*npts:]
    np = pthi-ptlo + 1
    kk = (np-1)/2
    
    selecti(kk, hp, np, cp)
    
    for k in range(ndim):
      region[0,k] = regions[tmom,0,k]
      region[1,k] = regions[tmom,1,k]
      
    region[0,tdim] = coords[tdim*npts + hp[kk]]
    region[1,tdim] = coords[tdim*npts + hp[kk]]
  
    jbox += 1
    for k in range(ndim):
      regions[jbox,0,k] = regions[tmom,0,k]
      regions[jbox,1,k] = region[1,k]
    
    parent[jbox] = tmom
    child1[jbox] = 0
    child2[jbox] = 0
    boxptlo[jbox] = ptlo
    boxpthi[jbox] = ptlo+kk
    
    jbox += 1
    for k in range(ndim):
      regions[jbox,0,k] = region[0,k]
      regions[jbox,1,k] = regions[tmom,1,k]
    
    parent[jbox] = tmom
    child1[jbox] = 0
    child2[jbox] = 0
    boxptlo[jbox] = ptlo+kk+1
    boxpthi[jbox] = pthi
    
    child1[tmom] = jbox-1
    child2[tmom] = jbox
    
    if kk > 1:
      nowtask += 1
      taskmom[nowtask] = jbox-1
      taskdim[nowtask] = (tdim+1) % ndim
      
    if (np-kk) > 3:
      nowtask += 1
      taskmom[nowtask] = jbox
      taskdim[nowtask] = (tdim+1) % ndim
      
  
  for j in range(npts):
    rptindx[ptindx[j]] = j
    
  del coords 
  
  return ptindx, rptindx, regions, parent, child1, child2, boxptlo, boxpthi

    
cdef _nearest(ndarray[double, ndim=2] query, ndarray[double, ndim=2] points,
              ndarray[int, ndim=1] ptindx, ndarray[double, ndim=3] regions,
              ndarray[int, ndim=1] child1, ndarray[int, ndim=1] child2,
              ndarray[int, ndim=1] ptlo, ndarray[int, ndim=1] pthi):
   
  cdef int i, j, k, nrst, ntask, ndim = points.shape[1]
  cdef int d1, jdim, nq = len(query)
  cdef ndarray[int, ndim=1] indices = empty(nq, int32)
  cdef ndarray[int, ndim=1] task = empty(50, int32)
  cdef ndarray[double, ndim=1] pt = empty(3, float)
  cdef double d, dnrst  
  
  for j in range(nq):
    pt = query[j]
    dnrst = BIG
    k = jdim = 0
 
    while child1[k]:
      d1 = child1[k]
 
      if pt[jdim] <= regions[d1, 1, jdim]:
        k = d1
      else:
        k = child2[k]
 
      jdim += 1
      jdim = jdim % ndim
 
    for i in range(ptlo[k], pthi[k]+1):
 
      d = dist2(ndim, points[ptindx[i]], pt)
 
      if d < dnrst:
        nrst = ptindx[i]
        dnrst = d
 
    task[1] = 0
    ntask = 1
 
    while ntask:
      k = task[ntask]
      ntask -= 1
 
      if dist2box(ndim, regions[k], pt) < dnrst:
        if child1[k]:
          ntask += 1
          task[ntask] = child1[k]
          ntask += 1
          task[ntask] = child2[k]
 
        else:
          for i in range(ptlo[k], pthi[k]+1):
            d = dist2(ndim, points[ptindx[i]], pt)
 
            if d < dnrst:
              nrst = ptindx[i]
              dnrst = d

    indices[j] = nrst
    
    
  return indices  
    
  
