from libc.math cimport abs, sqrt, ceil, floor, log, log2, acos, cos
from numpy cimport ndarray
import numpy as np
from numpy import ones, zeros, int32, float32, uint8, fromstring
from numpy import sort, empty, array, arange, concatenate, searchsorted
import cython, sys
from scipy.spatial import cKDTree
from time import time

from .kdTree import KdTree

cdef double BIG = sys.float_info.max
cdef int BIGINT = 4294967295
    
cdef void cross3d(double v1[3], double v2[3], double v3[3]):

  v1[0] = v2[1]*v3[2] - v2[2]*v3[1]
  v1[1] = v2[2]*v3[0] - v2[0]*v3[2]
  v1[2] = v2[0]*v3[1] - v2[1]*v3[0]


cdef double dot3d(double v1[3], double v2[3]):
  
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]


cdef void diff3d(double v1[3], double v2[3], double v3[3]):
  
  v1[0] = v2[0] - v3[0]
  v1[1] = v2[1] - v3[1]
  v1[2] = v2[2] - v3[2]
    
    
cdef calcTorsionAngle(double v1[3], double v2[3], double v3[3], double v4[3]):
  
  cdef double b12[3]
  cdef double b32[3]
  cdef double b43[3]
  cdef double p13[3]
  cdef double p24[3]
  cdef double angle, sz
  
  diff3d(b12, v1, v2)
  diff3d(b32, v3, v2)
  diff3d(b43, v4, v3)

  cross3d(p13, b12, b32)
  cross3d(p24, b43, b32)
  
  sz =  sqrt( dot3d(p13, p13) * dot3d(p24, p24))
  
  if sz > 0.0:
    angle = dot3d(p13, p24) / sz
    angle = min(1.0, max(-1.0, angle))
    angle = acos(angle)
  
  else:
    angle = 0.0
  
  cross3d(b43, p24, b32)
  
  if dot3d(p13, b43) < 0:
    angle = -angle

  return angle 
  
  
  
def getDihedralAngles(ndarray[double, ndim=3] coords):

  cdef int i, j, k
  cdef int m = len(coords)
  cdef int n = len(coords[0])
  cdef double v1[3]
  cdef double v2[3]
  cdef double v3[3]
  cdef double v4[3]
  cdef ndarray[double, ndim=2] angles = empty((m, n-4), float)
  
  for i in range(m):
    for j in range(n-4):
      for k in range(3):
        v1[k] = coords[i,j,  k]
        v2[k] = coords[i,j+1,k]
        v3[k] = coords[i,j+2,k]
        v4[k] = coords[i,j+3,k]
      
      angles[i,j] = calcTorsionAngle(v1, v2, v3, v4)
  
  return angles
  
  
def getRadialClipRegions(ndarray[double, ndim=3] coords,
                         ndarray[double, ndim=1] origin,
                         double radius):
  
  cdef int i, j, k = 0
  cdef int nm = coords.shape[0]
  cdef int np = coords.shape[1]
  cdef int nr = 0
  cdef double dx, dy, dz, d2
  cdef double r2 = radius ** 2
  
  cdef ndarray[int, ndim=2] regions = zeros((np, 2), int32)
  
  for i in range(np):
  
    for j in range(nm):
      dx = coords[j,i,0] - origin[0]
      dy = coords[j,i,1] - origin[1]
      dz = coords[j,i,2] - origin[2]
      d2 = dx*dx + dy*dy + dz*dz
      
      if d2 > r2:
        break
      
    else:
      
      if nr > 0:
        if regions[k,1] == i-1:
          regions[k,1] = i
          
        else:
          if regions[k,1] > regions[k,0]:
            k += 1
            nr += 1
 
          regions[k,0] = i
          regions[k,1] = i
      
      else:
        regions[k,0] = i
        regions[k,1] = i
        nr += 1
  
  if regions[k,1] == regions[k,0]:
    nr -= 1
        
  return regions[:nr]

def getMinDists(ndarray[double, ndim=2] coords_a,
                ndarray[double, ndim=2] coords_b):

  cdef int i
  cdef int na = len(coords_a)
  cdef int nb = len(coords_b)

  cdef ndarray[double, ndim=1] min_dists = zeros(nb, float)
  
  kt = cKDTree(coords_b, 10)
  
  dists, idx = kt.query(coords_a, k=1)
  
  
  
  
def countDistances(ndarray[double, ndim=2] coords_a,
                   ndarray[double, ndim=2] coords_b,
                   ndarray[long, ndim=1] seq_a,
                   ndarray[long, ndim=1] seq_b,
                   ndarray[double, ndim=1] max_dists,
                   int min_seq_sep):

  cdef int i, j, k, nb
  cdef double dx, dy, dz, dd
  cdef int na = len(coords_a)
  cdef int nc = len(max_dists)
  
  cdef ndarray[double, ndim=1] max_dists_2 = max_dists ** 2
  cdef ndarray[int, ndim=2] counts = zeros((nc, na), int32)
  cdef ndarray[int, ndim=1] idx

  kt = cKDTree(coords_b, 10)
  idx_list = kt.query_ball_point(coords_a, max_dists.max())
      
  for i in range(na):
    idx = array(idx_list[i], int32)
    nb = len(idx)
    
    for j in range(nb):
      
      if abs(seq_a[i] - seq_b[idx[j]]) < min_seq_sep:
        continue
      
      dx = coords_a[i, 0] - coords_b[idx[j], 0]
      dy = coords_a[i, 1] - coords_b[idx[j], 1]
      dz = coords_a[i, 2] - coords_b[idx[j], 2]
      dd = dx*dx + dy*dy + dz*dz
      
      for k in range(nc):
        if dd <= max_dists_2[k]:
          counts[k,i] +=1
      
  return counts
  

def binPointValueSeqSeps(ndarray[int, ndim=1] ref_points,
                         ndarray[int, ndim=1] exp_points,
                         ndarray[double, ndim=1] exp_values,
                         py_bin_size, py_n_bins):

  cdef int i, j, ref, delta
  cdef int bin_size = py_bin_size
  cdef int n_bins = py_n_bins
  cdef int ne = len(exp_points)
  cdef int nr = len(ref_points)
  cdef int seq_width = n_bins * bin_size
  
  cdef ndarray[int, ndim=1] closest = searchsorted(exp_points, ref_points).astype(int32)
  cdef ndarray[double, ndim=1] hist = zeros(n_bins, float)


  for i in range(nr):
    ref = ref_points[i]
  
    j = closest[i]
    
    if j < ne-1:
    
      delta = abs(exp_points[j] - ref)
 
      while delta < seq_width:
        hist[delta/bin_size] += exp_values[j]
 
        j += 1
        
        if not j < ne:
          break
        
        delta = exp_points[j] - ref
        
    j = closest[i]-1
    
    if j > 0:
      delta = abs(ref - exp_points[j])
      
      while delta < seq_width:
        hist[delta/bin_size] += exp_values[j]
 
        j -= 1
        
        if j < 0:
          break
        
        delta = ref - exp_points[j]
  
  return hist
  

def getNeighbourSeqSeps(ndarray[int, ndim=1] seq_pos_a,
                        ndarray[int, ndim=1] seq_pos_b,
                        max_neighbour):

  cdef int i, j, j_min, k, d1, d2, d_min, step = 100
  cdef int na = len(seq_pos_a)
  cdef int nb = len(seq_pos_b)
  cdef int m = max_neighbour

  cdef ndarray[int, ndim=1] order = array(seq_pos_b.argsort(), int32)   
  cdef ndarray[int, ndim=2] seq_seps = zeros((m, na), int32)
  
  for i in range(na):
    j_min = 0

    k = 0
    
    while seq_pos_a[i] < seq_pos_b[order[k]]:
      k += step
      
      if k >= nb:
        break
    
    k = max(0, k-step)
    d_min = abs(seq_pos_a[i] - seq_pos_b[order[k]])
    
    for j in range(k, nb):
      
      d1 = abs(seq_pos_a[i] - seq_pos_b[order[j]])
      
      if d1 == 0:
        continue
      
      if d1 < d_min:
        d_min = d1
        j_min = j
      
      if d1 > d_min: # gone past, getting worse
        break
    
    seq_seps[0, i] = d_min # best
    
    # next index offsets to check
    p = 1 # downstream
    q = 1 # upstream
    
    for k in range(1, m): # other ranks
      
      if j_min+q < nb:
        d2 = abs(seq_pos_a[i] - seq_pos_b[order[j_min+q]])
        
      else:
        d2 = BIGINT
        
      if j_min-p > 0:
        d1 = abs(seq_pos_a[i] - seq_pos_b[order[j_min-p]])
      else:
        d1 = BIGINT
      
      if d1 == 0:
        seq_seps[k, i] = d2
        q += 1
      
      elif d2 == 0:
        seq_seps[k, i] = d1
        p += 1
      
      elif d1 < d2:
        seq_seps[k, i] = d1
        p += 1
      
      else:
        seq_seps[k, i] = d2
        q += 1
        
  return seq_seps  
  

def getInvDistSums(ndarray[double, ndim=2] coords_a,
                   ndarray[double, ndim=2] coords_b,
                   py_seq_pos_a=None,
                   py_seq_pos_b=None,
                   py_min_seq_sep=None, 
                   py_values_b=None,
                   power_adj=1):

  cdef int i, j
  cdef double dx, dy, dz, d
  cdef int p = power_adj
  cdef int na = len(coords_a)
  cdef int nb = len(coords_b)
  cdef long min_seq_sep
  
  cdef ndarray[double, ndim=1] inv_dist_sums = zeros(na, float)
  cdef ndarray[double, ndim=1] values_b
  cdef ndarray[long, ndim=1] seq_pos_a
  cdef ndarray[long, ndim=1] seq_pos_b
  
  if py_values_b is None:
    values_b = ones(nb, float)
  else:
    values_b = py_values_b
  
  if py_min_seq_sep is None:
  
    for i in range(na):
      for j in range(nb):
 
        dx = coords_a[i, 0] - coords_b[j, 0]
        dy = coords_a[i, 1] - coords_b[j, 1]
        dz = coords_a[i, 2] - coords_b[j, 2]
 
        #d = sqrt(dx*dx + dy*dy + dz*dz)
        d = dx*dx + dy*dy + dz*dz
        
        if p > 1:
          d = d ** p
        
        if d > 0:
          inv_dist_sums[i] += values_b[j]/d
  
  else:
    min_seq_sep = py_min_seq_sep
    seq_pos_a = py_seq_pos_a
    seq_pos_b = py_seq_pos_b
  
    for i in range(na):
      for j in range(nb):
        if abs(seq_pos_a[i]-seq_pos_b[j]) < min_seq_sep:
          continue
        
        dx = coords_a[i, 0] - coords_b[j, 0]
        dy = coords_a[i, 1] - coords_b[j, 1]
        dz = coords_a[i, 2] - coords_b[j, 2]
 
        d = dx*dx + dy*dy + dz*dz
        
        if p > 1:
          d = d ** p
 
        if d > 0:
          inv_dist_sums[i] += values_b[j]/d  
  
  return inv_dist_sums


def getDistanceMatrix(ndarray[double, ndim=2] coords_a,
                      ndarray[double, ndim=2] coords_b,
                      max_dist=BIG):

  cdef int i, j
  cdef double dx, dy, dz
  cdef int na = len(coords_a)
  cdef int nb = len(coords_b)
  
  cdef ndarray[double, ndim=2] matrix = empty((na, nb), float)
  
  for i in range(na):
  
    for j in range(nb):
      
      dx = coords_a[i, 0] - coords_b[j, 0]
      dy = coords_a[i, 1] - coords_b[j, 1]
      dz = coords_a[i, 2] - coords_b[j, 2]
      
      matrix[i,j] = min(sqrt(dx*dx + dy*dy + dz*dz), max_dist)
  
  return matrix
  

def getDistanceQuantiles(py_seq_seps, py_dists, py_is_cis,
                         ndarray[double, ndim=2] cis_dist_distrib,
                         ndarray[double, ndim=1] trans_dist_distrib,
                         particle_sep=100000, dist_width=0.1,):
                         
  cdef ndarray[int, ndim=1] seq_seps = array(py_seq_seps, int32)
  cdef ndarray[double, ndim=1] dists = array(py_dists, float)
  cdef ndarray[int, ndim=1] is_cis = array(py_is_cis, int32)
  
  cdef int i, j, k, b
  cdef int nd = dists.shape[0]
  cdef int nb = cis_dist_distrib.shape[0]
  cdef int nc = cis_dist_distrib.shape[1]
  cdef int nt = trans_dist_distrib.shape[0]
  
  cdef double dw = dist_width
  cdef double f, q, t
  cdef double ps = float(particle_sep)
  cdef ndarray[double, ndim=1] quantiles = empty(nd, float)
  
  for i in range(nd):
    
    f = dists[i] / dw
    b = <int>(f) # bin of dist within the distribution
    
    if seq_seps[i] < 1:
      q = 0.0
      t = 1.0
    
    elif is_cis[i]:

      k = <int>ceil(seq_seps[i]/ps) # bead separation
      
      if k < 1:
        k = 1
      
      k = <int>floor(log2(k) * 10.0) # log2 bead bin
      
      if k >= nb:
        k = nb - 1
      
      q = (f-b) * cis_dist_distrib[k, b]
      t = 0.0
      
      for j in range(nc):
        if j < b:
          q += cis_dist_distrib[k, j]
        
        t += cis_dist_distrib[k, j] 
            
    else:
      
      q = (f-b) * trans_dist_distrib[b]
      t = 0.0
      
      for j in range(nt):
        if j < b:
          q += trans_dist_distrib[j]
        
        t += trans_dist_distrib[j]

      
    quantiles[i] = 100.0 * q / t  
  
  
  return quantiles


def pairsDualRegionIntersection(ndarray[int, ndim=2] point_pairs,
                                ndarray[int, ndim=2] regions_a,
                                ndarray[int, ndim=2] regions_b,
                                exclude=False):
  
  cdef int i, j, k, a, b
  cdef int sel_overlap = 1 - int(exclude)
  cdef int ni = 0
  cdef int np = len(point_pairs)
  cdef int nra = len(regions_a)
  cdef int nrb = len(regions_b)
  
  cdef ndarray[int, ndim=1] indices = empty(np, int32)
  cdef ndarray[int, ndim=1] order_a = array(regions_a[:,0].argsort(), int32)  
  cdef ndarray[int, ndim=1] order_b = array(regions_b[:,0].argsort(), int32)  
  
  for i in range(np):

    if point_pairs[i,0] < regions_a[order_a[0],0]:
      if not sel_overlap:
        indices[ni] = i
        ni += 1
      
      continue
    
    if point_pairs[i,1] < regions_b[order_b[0],0]:
      if not sel_overlap:
        indices[ni] = i
        ni += 1
      
    a = 0
    b = 0
    
    for k in range(nra):
      j = order_a[k]
       
      if (regions_a[j,0] <= point_pairs[i,0]) and (point_pairs[i,0] <= regions_a[j,1]):
        a = 1
        break
 
      if point_pairs[i, 0] < regions_a[j, 0]:
        break

    for k in range(nrb):
      j = order_b[k]
      
      if (regions_b[j,0] <= point_pairs[i,1]) and (point_pairs[i,1] <= regions_b[j,1]):
        b = 1
        break
 
      if point_pairs[i, 1] < regions_b[j, 0]:
        break
    
    if sel_overlap == a & b:
      indices[ni] = i
      ni += 1
  
  return indices[:ni]


def pointRegionsIntersection(ndarray[int, ndim=1] pos,
                             ndarray[int, ndim=2] regions,
                             exclude=False):
  
  cdef int i, j, k, a, b
  cdef int sel_overlap = 1 - int(exclude)
  cdef int ni = 0
  cdef int np = len(pos)
  cdef int nr = len(regions)
  cdef ndarray[int, ndim=1] indices = empty(np, int32)
  cdef ndarray[int, ndim=1] order = array(regions[:,0].argsort(), int32)  
  
  for i in range(np):
    
    if pos[i] < regions[order[0],0]:
      if not sel_overlap:
        indices[ni] = i
        ni += 1
      
      continue
      
    a = 0
    for k in range(nr):
      j = order[k]
      
      if (regions[j,0] <= pos[i]) and (pos[i] <= regions[j,1]):
        a = 1
        break
 
      if pos[i] < regions[j, 0]:
        break
        
    if sel_overlap == a:
      indices[ni] = i
      ni += 1
  
  return indices[:ni]


def pairRegionsIntersection(ndarray[int, ndim=2] pairs,
                            ndarray[int, ndim=2] regions,
                            exclude=False, allow_partial=False,
                            region_indices=False):
  
  cdef int i, j, k, a, b
  cdef int exc = int(exclude)
  cdef int partial = int(allow_partial)
  cdef int ni = 0
  cdef int np = len(pairs)
  cdef int nr = len(regions)
  cdef ndarray[int, ndim=1] indices = empty(np, int32)
  cdef ndarray[int, ndim=1] indices_reg = empty(np, int32)
  cdef ndarray[int, ndim=1] order = array(regions[:,0].argsort(), int32)  
  
  for i in range(np):
    
    if pairs[i,1] < regions[order[0],0]:
      if exc:
        indices[ni] = i
        ni += 1
      
      continue

    if pairs[i,0] < regions[order[0],0]:
      if exc and partial:
        indices[ni] = i
        ni += 1
      
      continue
      
    a = 0
    b = 0
    
    for k in range(nr):
      j = order[k]
      #print i, j, k
      
      if (regions[j,0] <= pairs[i,0]) and (pairs[i,0] <= regions[j,1]):
        a = 1
    
      if (regions[j,0] <= pairs[i,1]) and (pairs[i,1] <= regions[j,1]):
        b = 1
 
      if (pairs[i, 0] < regions[j, 0]) and (pairs[i, 1] < regions[j, 0]):
        break
      
      if partial & (a | b):
        break
      elif a & b:
        break
    
    if partial:
      if exc and not (a & b):
        indices[ni] = i
        indices_reg[ni] = j
        ni += 1
      
      elif a | b:
        #print a, b, i, ni, regions[j,0], pairs[i,0], pairs[i,1], regions[j,1]
        indices[ni] = i
        indices_reg[ni] = j
        ni += 1
      
    else:
      if exc and not (a | b):
        indices[ni] = i
        indices_reg[ni] = j
        ni += 1
      
      elif a & b:
        indices[ni] = i
        indices_reg[ni] = j
        ni += 1

  
  if region_indices:
    return indices[:ni], indices_reg[:ni]
    
    
  else:
    return indices[:ni]

  
def calcSeqEntropy(seq, window=100):
  
  cdef int i, j, w = window
  cdef int n = len(seq)
  cdef int m = 5
  cdef int n2 = n-w
  
  cdef double c, h, b = 0.25
  cdef double fw = float(window)
  
  cdef ndarray[int, ndim=1] counts = zeros(m, int32)
  cdef ndarray[int, ndim=1] idx = zeros(26, int32)
  cdef ndarray[double, ndim=1] entropy = ones(n, float)
  cdef ndarray[char, ndim=1] iseq = fromstring(seq, dtype=uint8) - 65
  
  for k, a in enumerate('ACGNT'):
    idx[ord(a) - ord('A')] = k
  
  for i in range(w):
    counts[idx[iseq[i]]] += 1
  
  h = 0.0
  for j in range(m):
    if counts[j] > 0:
      c = counts[j]/fw
      h += c * log2(c/b)
  
  for i in range(w):
    entropy[i] = h
    
  for i in range(1, n2):
    
    if i % 1000000 == 0:
      print i
    
    counts[idx[iseq[i-1]]] -= 1
    counts[idx[iseq[i+w]]] += 1
    
    h = 0.0
    for j in range(m):
      if counts[j] > 0:
        c = counts[j]/fw
        h += c * log2(c/b)
    
    entropy[i] = h

  for i in range(n2,n):
    entropy[i] = h
  
  return entropy  


def calcEnsembleRadGyration(ndarray[double, ndim=3] coords, report_mean=True):
  
  cdef int i, j
  cdef int n = len(coords[0])
  cdef int m = len(coords)
  cdef int do_mean = 1 if report_mean else 0
  cdef double dx, dy, dz, xm, ym, zm, rg
  cdef double nf = float(n)
  cdef ndarray[double, ndim=1] rogs = zeros(m, float)
  
  
  for j in range(m):
    rg = 0.0
    xm = 0.0
    ym = 0.0
    zm = 0.0
 
    for i in range(n):
      xm += coords[j, i, 0]
      ym += coords[j, i, 1]
      zm += coords[j, i, 2]

    xm /= nf
    ym /= nf
    zm /= nf

    for i in range(n):
      dx = coords[j, i, 0] - xm
      dy = coords[j, i, 1] - ym
      dz = coords[j, i, 2] - zm
 
      rg += dx*dx + dy*dy + dz*dz
    
    rogs[j] = rg/nf
      
  if report_mean:
    return rogs.mean()
  else:
    return rogs 
  

def calcCoordsRadGyration(ndarray[double, ndim=2] coords):
  
  cdef int i
  cdef int n = len(coords)
  cdef double dx, dy, dz, xm, ym, zm, rg
  cdef double nf = float(n)

  rg = 0.0
  xm = 0.0
  ym = 0.0
  zm = 0.0
  
  for i in range(n):
    xm += coords[i, 0]
    ym += coords[i, 1]
    zm += coords[i, 2]

  xm /= nf
  ym /= nf
  zm /= nf

  for i in range(n):
    dx = coords[i, 0] - xm
    dy = coords[i, 1] - ym
    dz = coords[i, 2] - zm
  
    rg += dx*dx + dy*dy + dz*dz
    
  rg /= nf

  return rg


def calcChainRadGyration(ndarray[double, ndim=2] coords, numPoints=11):

  # sliding along the coordinate chain
  # for a given number of particles  
  # get mean position as centre of mass
  # RMS distance to centre of mass of those points
  
  cdef int i, j
  cdef int np = numPoints
  cdef int n = len(coords)
  cdef int mp = np/2
  cdef double dx, dy, dz, xm, ym, zm, rg
  cdef double npf = float(numPoints)
  cdef ndarray[double, ndim=1] rgArray = zeros(n, float)
  
  for i in range(n-np):
    rg = 0.0
    xm = 0.0
    ym = 0.0
    zm = 0.0
    
    for j in range(i, i+np):
      xm += coords[j, 0]
      ym += coords[j, 1]
      zm += coords[j, 2]

    xm /= npf
    ym /= npf
    zm /= npf

    for j in range(i, i+np):
      dx = coords[j, 0] - xm
      dy = coords[j, 1] - ym
      dz = coords[j, 2] - zm
   
      rg += dx*dx + dy*dy + dz*dz
      
    rg /= npf
    rgArray[i+mp] = rg
  
  # Extend values to edges, otherwise value represent middle of window
  for i in range(mp):
    rgArray[i] = rgArray[mp]
  
  for i in range(n-mp-1, n):
    rgArray[i] = rgArray[n-mp-2]

  return rgArray
      

def calcCoordMesh(ndarray[double, ndim=2] coords,
                  ndarray[double, ndim=2] inColors,
                  radius=5.0, nVoxels=100):
  
  from cUtil.contour import contourer3d
  
  cdef int i, j, k, x, y, z, n, nv = nVoxels
  cdef double size, r, g, b, a
  cdef double xf, yf, zf, s, rad = radius
  cdef ndarray[double, ndim=3] voxels
  cdef ndarray[double, ndim=4] vcolors
  cdef ndarray[double, ndim=4] vnorms
  cdef ndarray[double, ndim=1] offset, vertices, colors, normals
  
  #t0 = time()  
  voxels, vcolors, vnorms, offset, size = calcVoxels(coords, inColors, rad,  nv)
  
  #t1 = time()  
  vertices = contourer3d(voxels, 0.51)
  
  #t2 = time()  
  n = len(vertices)
  normals = empty(n, float)
  colors  = empty((n*4)/3, float)
   
  i = 0
  j = 0
  while i < n-2:
    xf = vertices[i]  
    yf = vertices[i+1]
    zf = vertices[i+2]
    
    vertices[i]   = (xf * size) + offset[0]
    vertices[i+1] = (yf * size) + offset[1]
    vertices[i+2] = (zf * size) + offset[2]

    x = <int>(xf + 0.5)
    y = <int>(yf + 0.5)
    z = <int>(zf + 0.5)
    
    colors[j]   = vcolors[x,y,z,0]
    colors[j+1] = vcolors[x,y,z,1]
    colors[j+2] = vcolors[x,y,z,2]
    colors[j+3] = vcolors[x,y,z,3]
    
    normals[i]   = vnorms[x,y,z,0]
    normals[i+1] = vnorms[x,y,z,1]
    normals[i+2] = vnorms[x,y,z,2]
    
    i += 3
    j += 4
 
  t3 = time()  
  
  #print "  Sub times:%.3f %.3f %.3f" % (t1-t0, t2-t1, t3-t2)
  
  return vertices, colors, normals


def calcCoordDepth(ndarray[double, ndim=2] coords, radius=5.0, nVoxels=100, chromoSizes=None, chromoDepth=False):
  
  cdef int i, n = len(coords)
  cdef int v = nVoxels
  cdef double r = radius
  cdef double dx, dy, dz, dmin = BIG
  
  cdef ndarray[double, ndim=2] surface
  cdef ndarray[double, ndim=2] nearest
  cdef ndarray[double, ndim=1] depths = empty(n, float)
  cdef ndarray[int, ndim=1] chromo_sizes
  
  if chromoSizes is not None:
    chromo_sizes = array(chromoSizes, int32)
    surface = calcTransInterface(coords, chromo_sizes, r, v, int32(chromoDepth))
    
  else:
    surface = calcSurface(coords, r, v)
  
  try:
    from scipy.spatial import cKDTree
    kt = cKDTree(surface, 10)
    dists, idx = kt.query(coords)
    nearest = surface[idx]
    
  except:  
    kt = KdTree(surface)
    nearest = kt.getNearest(coords)
  
  for i in range(n):
    dx = coords[i,0] - nearest[i,0]
    dy = coords[i,1] - nearest[i,1]
    dz = coords[i,2] - nearest[i,2]
    
    depths[i] = sqrt(dx*dx + dy*dy + dz*dz)
    
    if depths[i] < dmin:
      dmin = depths[i]
  
  for i in range(n):
    depths[i] -= dmin
  
  return depths


def calcBoundingBox(ndarray[double, ndim=2] coords):

  cdef int i, n = len(coords)
  cdef ndarray[double, ndim=1] bbMin = zeros(3, float)
  cdef ndarray[double, ndim=1] bbMax = zeros(3, float)
  cdef ndarray[double, ndim=1] bbSize = zeros(3, float)
  
  for k in range(3):
    bbMin[k] = BIG
    bbMax[k] = -BIG
  
  # Bounding bb
  for i in range(n):
  
    for k in range(3):
    
      if coords[i,k] > bbMax[k]:
        bbMax[k] = coords[i,k]
        
      elif coords[i,k] < bbMin[k]:
        bbMin[k] = coords[i,k]
 
  for k in range(3):
    bbSize[k] = bbMax[k]-bbMin[k]
  
  return bbMin, bbMax, bbSize
  

def _getSphereVolumeCoords(int r):

  cdef int x, y, z, d2, r2, ne
  cdef ndarray[int, ndim=2] volume = empty((8*r*r*r, 3), int32)
  
  r2 = r * r
  ne = 0
  
  for x in range(-r, r+1):
    for y in range(-r, r+1):
      for z in range(-r, r+1):
        d2 = x*x + y*y + z*z

        if d2 <= r2:
          volume[ne, 0] = x
          volume[ne, 1] = y
          volume[ne, 2] = z
          ne += 1
  
  return ne, volume
  
  
def calcSurfaceBuried(ndarray[double, ndim=2] surfCoords, # Those that define a surface, e.g. one chromsome
                      ndarray[double, ndim=2] testCoords, # Those to test, e.g. other chromosomes
                      radius=5.0, nVoxels=100, nIntersect=4):
  """
  Determines which points within an array of test oordinates are within the surface, as
  defined by a radius, of another set of coordinates.
  Returns buried indices of test coords
  """
  
  cdef int i, j, k, x, y, z, x1, y1, z1, w
  cdef int nv, nb, n = len(surfCoords)
  cdef int nt = len(testCoords)
  cdef int ni = nIntersect
  cdef ndarray[int, ndim=2] volume
  cdef ndarray[double, ndim=1] bbMin = zeros(3, float)
  cdef ndarray[double, ndim=1] bbMax = zeros(3, float)
  cdef ndarray[double, ndim=1] bbSize = zeros(3, float)
  cdef double size, dx, dy, dz, d2, r2 = 2*radius
  cdef ndarray[int, ndim=3] voxels = zeros((nVoxels, nVoxels, nVoxels), int32)
  cdef ndarray[double, ndim=2] testBuried = empty((nt, 3), float)
  cdef ndarray[int, ndim=1] buriedIdx = empty(nt, int32)
  
  bbMin, bbMax, bbSize = calcBoundingBox(surfCoords)
    
  # calc voxel size and voxel radius
  size = (2*r2+max([bbSize[0], bbSize[1], bbSize[2]])) / <double>nVoxels
  w = <int>ceil(radius/size) + 1
  
  # setup of volume of influence for each point
  nv, volume = _getSphereVolumeCoords(w)  

  # Fill voxels with volumes of spheres around each surface generating point  
  for i in range(n):
    # get voxel bin for each coord
    x = <int>( (surfCoords[i,0]-bbMin[0]+r2)/size )
    y = <int>( (surfCoords[i,1]-bbMin[1]+r2)/size )
    z = <int>( (surfCoords[i,2]-bbMin[2]+r2)/size )
 
    # add chomo volume values
    for j in range(nv):
      x1 = x+volume[j,0]
      y1 = y+volume[j,1]
      z1 = z+volume[j,2]
      voxels[x1, y1, z1] += 1
  
  # Dected test points that are in the filled volume
  nb = 0
  for i in range(nt):
    # get voxel bin for each coord
    
    x = <int>( (testCoords[i,0]-bbMin[0]+r2)/size )
    if 0 <= x < nVoxels:
      
      y = <int>( (testCoords[i,1]-bbMin[1]+r2)/size )
      if 0 <= y < nVoxels:
      
        z = <int>( (testCoords[i,2]-bbMin[2]+r2)/size )
        if 0 <= z < nVoxels:
        
          if voxels[x, y, z] >= ni:
            testBuried[nb, 0] = testCoords[i, 0]
            testBuried[nb, 1] = testCoords[i, 1]
            testBuried[nb, 2] = testCoords[i, 2]
            buriedIdx[nb] = i
            nb += 1

  testBuried = testBuried[:nb]
  buriedIdx = buriedIdx[:nb]
  
  # For the buried short-list of points test to see whether any are within radius
  # of surface coords without any binning using a KD-tree
  
  try:
    from scipy.spatial import cKDTree
    kt = cKDTree(surfCoords, 10)
    dists, idx = kt.query(testBuried)
    nearest = surfCoords[idx]
    
  except:  
    kt = KdTree(surfCoords)
    nearest = kt.getNearest(testBuried)
  
  j = 0
  for i in range(nb):
    k = buriedIdx[i]
    dx = testBuried[i, 0] - nearest[i, 0]
    dy = testBuried[i, 1] - nearest[i, 1]
    dz = testBuried[i, 2] - nearest[i, 2]
    
    d2 = dx*dx + dy*dy + dz*dz
    
    if d2 < r2: # Point really is within surface raduis: burried
      buriedIdx[j] = k
      j += 1
 
  return buriedIdx[:j]
  
  
def calcVoxels(ndarray[double, ndim=2] coords, ndarray[double, ndim=2] colors, double radius, nVoxels=100):

  cdef int i, j, k, x, y, z, x1, y1, z1, w, w2, d2
  cdef int nv, n = len(coords)
  cdef ndarray[double, ndim=1] bbMin = zeros(3, float)
  cdef ndarray[double, ndim=1] bbMax = zeros(3, float)
  cdef ndarray[double, ndim=1] bbSize = zeros(3, float)
  cdef double v1[3]
  cdef double s, size, rw, r2 = 2*radius
  cdef ndarray[double, ndim=3] voxels = zeros((nVoxels, nVoxels, nVoxels), float)
  cdef ndarray[double, ndim=4] colorVox = zeros((nVoxels, nVoxels, nVoxels, 4), float)
  cdef ndarray[double, ndim=4] normVox  = zeros((nVoxels, nVoxels, nVoxels, 3), float)
  cdef ndarray[double, ndim=3] voxCount = zeros((nVoxels, nVoxels, nVoxels), float)
  cdef ndarray[double, ndim=1] offset

  bbMin, bbMax, bbSize = calcBoundingBox(coords)
 
  for k in range(3):
    bbMin[k] -= r2
    
  # calc voxel size and voxel radius
  size = (2*r2+max([bbSize[0], bbSize[1], bbSize[2]])) / <double>nVoxels
  rw = ceil(radius/size)
  w = <int>rw
  w2 = w*w

  cdef ndarray[int, ndim=2] volume = empty((8*w2*w, 3), int32)
  cdef ndarray[double, ndim=1] intensity = empty((8*w2*w), float)
  cdef ndarray[double, ndim=2] uNorms = empty((8*w2*w, 3), float)
  
  nv = 0
  for x in range(-w, w+1):
    for y in range(-w, w+1):
      for z in range(-w, w+1):
        d2 = x*x + y*y + z*z

        if d2 < w2:
          volume[nv, 0] = x
          volume[nv, 1] = y
          volume[nv, 2] = z
          
          v1[0] = <double>x
          v1[1] = <double>y
          v1[2] = <double>z
           
          # unit vector normal
          s = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]
          if s != 0.0:
            s = sqrt(s)
            v1[0] = v1[0] / s
            v1[1] = v1[1] / s
            v1[2] = v1[2] / s
          
          uNorms[nv, 0] = v1[0]
          uNorms[nv, 1] = v1[1]
          uNorms[nv, 2] = v1[2]
          
          if rw-s > 1.0:
            intensity[nv] = 1.0
          else:
            intensity[nv] = rw-s
          
          nv += 1
  
  #t0 = time()
  
  # fill volumes
  for i in range(n):
    # get voxel bin for each coord
    x = <int>( 0.5+(coords[i,0]-bbMin[0])/size )
    y = <int>( 0.5+(coords[i,1]-bbMin[1])/size )
    z = <int>( 0.5+(coords[i,2]-bbMin[2])/size )
 
    # add values over volume
    for j in range(nv):
      x1 = x+volume[j,0]
      y1 = y+volume[j,1]
      z1 = z+volume[j,2]
      voxels[x1, y1, z1] += intensity[j]
      
      # Set colour voxels and normals based on central point
      normVox[x1, y1, z1, 0] += uNorms[j, 0]
      normVox[x1, y1, z1, 1] += uNorms[j, 1]
      normVox[x1, y1, z1, 2] += uNorms[j, 2]
      colorVox[x1, y1, z1, 0] += colors[i, 0]
      colorVox[x1, y1, z1, 1] += colors[i, 1]
      colorVox[x1, y1, z1, 2] += colors[i, 2]
      colorVox[x1, y1, z1, 3] += colors[i, 3]
      voxCount[x1, y1, z1] += 1.0
  
  #t1 = time()
  
  # calc colour averages
  for x in range(nVoxels):
    for y in range(nVoxels):
      for z in range(nVoxels):
        nc = voxCount[x, y, z]

        if nc > 1.0:
          normVox[x, y, z, 0] /= nc
          normVox[x, y, z, 1] /= nc
          normVox[x, y, z, 2] /= nc
          colorVox[x, y, z, 0] /= nc
          colorVox[x, y, z, 1] /= nc
          colorVox[x, y, z, 2] /= nc
          colorVox[x, y, z, 3] /= nc
          
          if voxels[x, y, z] > 1.0:
            voxels[x, y, z] = 1.0
          
          #voxels[x, y, z] /= nc
          
  
  #t2 = time()
  
  #print '     : %.3f %.3f' % (t1-t0, t2-t1)

  offset = array([bbMin[0], bbMin[1], bbMin[2]])
  return voxels, colorVox, normVox, offset, size
  
  
def calcTransInterface(ndarray[double, ndim=2] coords, ndarray[int, ndim=1] chromoSizes,
                       double radius, int nVoxels=100, int includeExterior=0):

  cdef int i, j, x, y, z, x1, y1, z1, x2, y2, z2, w, d2
  cdef int ne, n = len(coords)
  cdef int c, d, nc = len(chromoSizes)
  
  cdef ndarray[double, ndim=1] bbMin = zeros(3, float)
  cdef ndarray[double, ndim=1] bbMax = zeros(3, float)
  cdef ndarray[double, ndim=1] bbSize = zeros(3, float)

  cdef double v1[3]
  cdef double size, r2 = 2*radius
  
  cdef ndarray[int, ndim=4] voxels = zeros((nVoxels, nVoxels, nVoxels, 2), int32)
  cdef ndarray[int, ndim=1] chromoIds = zeros(n, int32)
  
  cdef ndarray[double, ndim=2] surface
  cdef ndarray[int, ndim=2] volume
  
  j = 0
  for c in range(nc):
    d = c+1
    for i in range(j, j+chromoSizes[c]):
      chromoIds[i] = d
    
    j += chromoSizes[c]
    
  bbMin, bbMax, bbSize = calcBoundingBox(coords)
    
  # calc voxel size and voxel radius
  size = (2*r2+max([bbSize[0], bbSize[1], bbSize[2]])) / <double>nVoxels
  w = <int>ceil(radius/size)

  # setup of volume of influence for each point
  ne, volume = _getSphereVolumeCoords(w)  
  
  # fill volumes
  ns = 0
  d = -1 # index indicating different chromosome overlap
  
  for i in range(n):
    # get voxel bin for each coord
    x = <int>( (coords[i,0]-bbMin[0]+r2)/size )
    y = <int>( (coords[i,1]-bbMin[1]+r2)/size )
    z = <int>( (coords[i,2]-bbMin[2]+r2)/size )
    c = chromoIds[i]
 
    # add chomo volume values
    for j in range(ne):
      x1 = x+volume[j,0]
      y1 = y+volume[j,1]
      z1 = z+volume[j,2]
      ns += 1
      
      x2 = (x-x1) * (x-x1)
      y2 = (y-y1) * (y-y1)
      z2 = (z-z1) * (z-z1)
      d2 = x2 + y2 + z2
      
      if voxels[x1, y1, z1, 0] == 0:
        voxels[x1, y1, z1, 0] = c
        voxels[x1, y1, z1, 1] = d2
 
      elif voxels[x1, y1, z1, 0] == c: # same chromo, different coord
        voxels[x1, y1, z1, 1] = min(d2, voxels[x1, y1, z1, 1])
 
      elif voxels[x1, y1, z1, 1] == d2: # different chromo, equal dist
        voxels[x1, y1, z1, 0] = d # mark overlap
      
      elif d2 < voxels[x1, y1, z1, 1]: # different chromo is closest
        voxels[x1, y1, z1, 0] = c
        voxels[x1, y1, z1, 1] = d2
      
  surface = empty((ns, 3), float) # array of surface coords
  ns = 0
  for x in range(nVoxels):
    for y in range(nVoxels):
      for z in range(nVoxels):
        if voxels[x,y,z,0] == d:  #  Where chromosomes mismatched at equal distance
          surface[ns, 0] = (size * <double>x)-r2+bbMin[0]
          surface[ns, 1] = (size * <double>y)-r2+bbMin[1]
          surface[ns, 2] = (size * <double>z)-r2+bbMin[2]
          ns += 1
        
        elif includeExterior and voxels[x,y,z,0] == 0:  #  Outside chromosomes
          surface[ns, 0] = (size * <double>x)-r2+bbMin[0]
          surface[ns, 1] = (size * <double>y)-r2+bbMin[1]
          surface[ns, 2] = (size * <double>z)-r2+bbMin[2]
          ns += 1
        
  return array(surface[:ns])


def calcSurface(ndarray[double, ndim=2] coords, double radius, int nVoxels=100):

  cdef int i, j, x, y, z, x1, y1, z1, w, w2, wa, d2
  cdef int ne, nm, n = len(coords)
  cdef ndarray[double, ndim=1] bbMin = zeros(3, float)
  cdef ndarray[double, ndim=1] bbMax = zeros(3, float)
  cdef ndarray[double, ndim=1] bbSize = zeros(3, float)

  cdef double v1[3], size,r2 = 2*radius
  cdef ndarray[int, ndim=3] voxels = zeros((nVoxels, nVoxels, nVoxels), int32)
  cdef ndarray[double, ndim=2] surface

  bbMin, bbMax, bbSize = calcBoundingBox(coords)
    
  # calc voxel size and voxel radius
  size = (2*r2+max([bbSize[0], bbSize[1], bbSize[2]])) / <double>nVoxels
  w = <int>ceil(radius/size)
  w2 = w*w
  wa = w2-1

  cdef ndarray[int, ndim=2] edge = empty((8*w2*w, 3), int32)
  cdef ndarray[int, ndim=2] middle = empty((8*w2*w, 3), int32)
  
  ne = 0
  nm = 0
  for x in range(-w, w+1):
    for y in range(-w, w+1):
      for z in range(-w, w+1):
        d2 = x*x + y*y + z*z

        if d2 <= wa:
          middle[nm, 0] = x
          middle[nm, 1] = y
          middle[nm, 2] = z
          nm += 1
        
        elif d2 <= w2:
          edge[ne, 0] = x
          edge[ne, 1] = y
          edge[ne, 2] = z
          ne += 1
            
  # fill volumes
  ns = 0
  for i in range(n):
   # get voxel bin for each coord
   x = <int>( (coords[i,0]-bbMin[0]+r2)/size )
   y = <int>( (coords[i,1]-bbMin[1]+r2)/size )
   z = <int>( (coords[i,2]-bbMin[2]+r2)/size )
 
   # add edge values
   for j in range(ne):
     x1 = x+edge[j,0]
     y1 = y+edge[j,1]
     z1 = z+edge[j,2]
     ns += 1

     if voxels[x1, y1, z1] == 0:
       voxels[x1, y1, z1] = 1
 
   # add middle values
   for j in range(nm):
     x1 = x+middle[j,0]
     y1 = y+middle[j,1]
     z1 = z+middle[j,2]
     voxels[x1, y1, z1] = 2

  surface = empty((ns, 3), float) # array of surface coords
  ns = 0
  for x in range(nVoxels):
    for y in range(nVoxels):
      for z in range(nVoxels):
        if voxels[x,y,z] == 1:  #  An unsullied edge
          surface[ns, 0] = (size * <double>x)-r2+bbMin[0]
          surface[ns, 1] = (size * <double>y)-r2+bbMin[1]
          surface[ns, 2] = (size * <double>z)-r2+bbMin[2]
          ns += 1
 
  return array(surface[:ns])
  

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

  #print(posDict, prevPosDict)
  
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
        p2 = min(p1+1, m-1)
      
      else: #new pos is below p1
        p2 = p1
        p1 = max(0, p1-1)
        dMin = positions[i] - prevPositions[p1]
      
      #calculate coordinates
      if prevPositions[p2] == prevPositions[p1]:
        newCoords[i0+i, 0] = coords[j0+p1, 0]
        newCoords[i0+i, 1] = coords[j0+p1, 1]
        newCoords[i0+i, 2] = coords[j0+p1, 2]
 
      else: #interpolate
        f = <float>dMin/<float>(prevPositions[p2]-prevPositions[p1])
        g = 1.0 - f
        
        newCoords[i0+i, 0] = g * coords[j0+p1, 0] + f * coords[j0+p2, 0]
        newCoords[i0+i, 1] = g * coords[j0+p1, 1] + f * coords[j0+p2, 1]
        newCoords[i0+i, 2] = g * coords[j0+p1, 2] + f * coords[j0+p2, 2]
  
  return newCoords
  
  

def concatenateRestraints(restraintDict, posDict, ambigDict, seqScale, backboneLower=0.1, backboneUpper=1.1):
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
      ambiguity = ambigDict[chrA][chrB]
      
      n = restraints.shape[1]
      
      for i in range(n):
        indicesArray[m,0] = <int>restraints[0,i] + startA
        indicesArray[m,1] = <int>restraints[1,i] + startB
 
        boundsArray[m,0] = restraints[4,i] # lower
        boundsArray[m,1] = restraints[5,i] # upper
 
        ambigArray[m] = ambiguity[i] # Group size for next n ambigous restraints
        
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
  
  return -1


def calcRestraints(chromosomes, ambigChromos, contactDict, pyIsSingleCell,
                   pyBboneSep=[500000, 500000], domLstDict={},
                   float scale=1.0, float exponent=-0.33,
                   float lower=0.8, float upper=1.2,
                   pyMinCount=2, float maxPopDist=5.0,
                   pySubset=None):
    
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
  cdef int minCount = pyMinCount
     
  cdef int i, j, k, a, b, c, n, na, nb, n_ambig
  cdef int cFilter, cFilterLocal, binned
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
  cdef ndarray[int, ndim=2] limits   # shape: (chromoId, 2:[start, end])
  cdef ndarray[int, ndim=1] bboneSep # Min and max binSizes, given domain list
  cdef ndarray[int, ndim=1] subsets  # TO restrict contacts if they have a subset/model label
  
  if pySubset is None:
    cFilter = 0    
  else:
    cFilter = 1
    subsets = array(pySubset, int32)
  
  if pyBboneSep is None:
    binned = 0
  else:
    binned = 1
    bboneSep = array(pyBboneSep, int32)

  binLstDict = {}
  nContDict  = {}
  posDict    = {}
  ambigDict  = {}
  bboneDict  = {}
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
    
    for chrB in contactDict[chrA]:
      if chrB not in chromos:
        continue
        
      b = chrIdx[chrB]
      
      positionsA = posDict[chrA]
      positionsB = posDict[chrB]
      contacts = contactDict[chrA][chrB]
      n = len(contacts[0])
      
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
      
        # Note: sometimes a == b
        positionsA[counts[a]] = contacts[0,i]
        counts[a] += 1
        
        positionsB[counts[b]] = contacts[1,i]
        counts[b] += 1

  if (domLstDict is None) and binned: # Regular backbone bins, no domains
    avSep = (bboneSep[0] + bboneSep[1]) / 2
    for a in range(c):
      limits[a,0] = avSep * (limits[a,0]/avSep)
      limits[a,1] = avSep * <int>(ceil(<float>limits[a,1]/avSep))

  elif (domLstDict is not None) and binned: # Domain backbone bins
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

    if (domLstDict is not None) and binned: # Set regularly spaced seq positions only
      if len(binLstDict[chrA]) < 2:
        posDict[chrA] = array([limits[a,0],], int32)
      else:
        #spacer = ((limits[a,1] + bboneSep[1]) - limits[a,0]) / (len(binLstDict[chrA]) - 1)
        #posDict[chrA] = arange(limits[a,0], limits[a,1] + spacer, spacer, int32)
        posDict[chrA] = arange(limits[a,0], limits[a,1] + bboneSep[0], bboneSep[0], int32)
        
      bboneDict[chrA] = ones(len(posDict[chrA])+1, int32)
      
    else:
      positionsA = posDict[chrA][:counts[a]] # Chop unused allocation
      
      
      if binned and (domLstDict is None): # inject backbone spacers
        sep = (bboneSep[0] + bboneSep[1]) / 2
        positionsA = positionsA + arange(limits[a,0], limits[a,1]+sep, sep, int32)
      
      positionsA = sort(positionsA)
      
      n = len(positionsA)
 
      if n ==0:
        continue
 
      positionsB = empty(n, int32)
      positionsB[0] = positionsA[0]
      
      if binned and (domLstDict is None):
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
    
    if chrA in ambigChromos:
      n_ambig = 2
    else:
      n_ambig = 1  
      
    a = chrIdx[chrA]
    backboneA = bboneDict[chrA]
    restraintDict[chrA] = {}
    ambigDict[chrA] = {}
    
    for chrB in contactDict[chrA]:
      if chrB not in chromos:
        continue
    
      if chrA in ambigChromos:
        n_ambig *= 2  
        
      b = chrIdx[chrB]
      backboneB = bboneDict[chrB]
      
      contacts = contactDict[chrA][chrB]
      n = len(contacts[0])  
      ambuguity = ones(n, int32)
      restraints = empty((6, n), float)
      
      if cFilter and len(contacts) == 4:
        cFilterLocal = 1
      else:
        cFilterLocal = 0
 
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
      
        for i in range(n):
          if cFilterLocal:
            for j in subsets:
              if contacts[3,i] == j: # This subset included
                break
            
            else:
              continue
            
          j = bpToBin(contacts[0,i], binLstDict[chrA], na)
          k = bpToBin(contacts[1,i], binLstDict[chrB], nb)
          
          if j < 0:  # Not in any of the bin definitions...
            continue
          if k < 0:
            continue
          
          binMatrix[j,k] += contacts[2,i]
        
        #loop over all binned contacts, and calculate the constraint target distance
        #using a powerlaw function and the number of observations
        k = 0
        for i in range(na):
          for j in range(nb):
            if binMatrix[i,j] > 0:
              if binMatrix[i,j] < minCount:
                continue
              
              restraints[0,k] = <double>i #binA
              restraints[1,k] = <double>j #binB
              restraints[2,k] = 1.0 # Weighting not currently used
                
              if isSingleCell:
                v = scale * binMatrix[i,j] ** exponent
                restraints[3,k] = v #target value
                restraints[4,k] = v * lower #constraint lower bound
                restraints[5,k] = v * upper #constraint upper bound
              
              backboneA[i] = 0
              backboneB[j] = 0
              
              k += 1
        
        restraints = restraints[:,:k]
        
      else:
        positionsA = posDict[chrA]
        positionsB = posDict[chrB]
        
        na = len(positionsA)
        nb = len(positionsB)
        indices = empty((n,2), int32)
        
        # find which positional index each contact corresponds to
        for i in range(n):
          for j in range(na):
            if positionsA[j] == contacts[0,i]:
              indices[i,0] = j
              break
          
          for j in range(nb):
            if positionsB[j] == contacts[1,i]:
              indices[i,1] = j
              break        
        
        k = 0
        for i in range(n):
          if cFilterLocal:
            for j in subsets:
              if contacts[3,i] == j: # This subset included
                break
            
            else:
              continue
            
          if contacts[2,i] < minCount:
            continue
           
          restraints[0,k] = <double>indices[i,0]
          restraints[1,k] = <double>indices[i,1]
          restraints[2,k] = 1.0 # Weighting, not currently used
          
          #print i, indices[i,0], indices[i,1]
          
            
          if isSingleCell:
            v = scale * contacts[2,i] ** exponent
            restraints[3,k] = v
            restraints[4,k] = v * lower
            restraints[5,k] = v * upper
          
          k += 1
          
        if k != n:
          restraints = restraints[:,:k]
      
      restraintDict[chrA][chrB] = restraints
      ambigDict[chrA][chrB] = ambuguity
      
      """
      else:
        # Assume normalised for restriction site density
        minObs = obsArray.min()
        deltaObs = obsArray.max() - minObs
        th = 1.0/maxPopDist
 
        adjObs = (obsArray - minObs) / deltaObs # 0 .. 1
        adjObs *= 100.0-th                      # 0 .. ~100
        adjObs += th                            # th .. 100
 
        values = scale * power(adjObs, exponent)
        lowers = values * lower
        uppers = values * upper
 
        values += scale
        uppers += scale
        lowers += scale
      
        restraints = vstack([indicesA, indicesB, weights, values, lowers, uppers])
      """
      
  # restraints as long array - separate function to convert
  
  return restraintDict, posDict, ambigDict, bboneDict, binLstDict


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def getInterpolatedColors(ndarray[double, ndim=2] colors, nPoints, alpha=None, blend=True):

  cdef int i, j, k, c, p1, p2, delta
  cdef int bl = int(blend)
  cdef int np = nPoints
  cdef int nc = len(colors)
  cdef double blockSize, f, g, flixedAlpha 
  cdef ndarray[double, ndim=2] colorsOut = empty((np, 4))
  
  if alpha is None:
    c = 4    
  else:
    flixedAlpha = alpha
    c = 3

  if (nc < 2) or (np < 2):
    return colors[:1]
  
  if nc == np:
    for i in range(np):
      for j in range(c):
        colorsOut[i,j] = colors[i,j]
  
  else:
    for j in range(c):
      colorsOut[np-1,j] = colors[nc-1,j]    
 
    blockSize = <double>(np-1)/(nc-1)   # Num new points per input color region
    for k in range(nc-1):
      p1 = <int>(k * blockSize)         # First point in this region
      p2 = <int>((k+1) * blockSize)     # First point in next region
      delta = p2-p1
      
      if bl > 0: 
        for i in range(delta):
          f = <double>i/<double>delta
          g = 1.0 - f
 
          for j in range(c):
            colorsOut[p1+i,j] = g * colors[k,j] + f * colors[k+1,j]
      
      else:
        for i in range(delta):
          for j in range(c):
            colorsOut[p1+i,j] = colors[k,j]
        
      
  if c == 3:
    for i in range(np):
      colorsOut[i,3] = flixedAlpha
  
  return colorsOut

   
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def getInterpolatedCoords(ndarray[int, ndim=1] seqPos, ndarray[double, ndim=2] coords,
                          ndarray[int, ndim=1] queryPos, clipEnds=False):
  
  cdef int i, j, k, a, b, k2 = 0
  cdef int ne = 0
  cdef int ce = 1 if clipEnds else 0
  cdef int nPos = len(seqPos)
  cdef int nQuery = len(queryPos)
  cdef int d, dMin, dMax=seqPos.max()
  cdef double f
  cdef ndarray[double, ndim=2] coords2 = zeros((nQuery,3))  
  
  if len(seqPos) != len(coords):
    msg = 'Number of sequence positions must match number of coordinates'    
    raise Exception(msg)

  for k in range(nQuery):
    dMin = dMax
  
    for i in range(nPos):
      d = seqPos[i]-queryPos[k]
 
      if abs(d) < abs(dMin):
        a = i
        dMin = d
        
      else:
        break  

    if dMin == 0: # At exact coords pos
      b = a
      
    elif dMin > 0: 
      if a-1 < 0:
        if ce:
          continue
        else:
          ne += 1
        
      b = a
      a = max(0, a-1)
      dMin = seqPos[b]-seqPos[a] - dMin
   
    else: 
      if a+1 > nPos-1:
        if ce:
          continue
        else:
          ne += 1
        
      b = min(a+1, nPos-1)
      dMin = - dMin
      
    if a == b:  
      for j in range(3):
        coords2[k2,j] = coords[a,j]
    
    else:
      f = <float>dMin/<float>(seqPos[b]-seqPos[a])
      for j in range(3):
        coords2[k2,j] = (1.0-f) * coords[a,j] + f * coords[b,j] 
    
    k2 += 1
    
  if ne > 0:
    msg = "Warning: %d out-of-bounds positions were mapped to coordinate ends"
    print msg % ne
  
  else:
    coords2 = coords2[:k2]
  
  return coords2
  

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
                int offsetA, int offsetB, int binSize=2000000,
                int symm=0, int transpose=0, int split_side=0):
  
  cdef int i, a, b
  cdef int n, m, nCont = len(contacts[0])
  
  n = len(binMatrix)
  m = len(binMatrix[0])
  
  # clauses outside for speed
  
  if split_side < 0:
    if transpose:
      for i in range(nCont):
        b = (contacts[0,i]-offsetA)/binSize
        a = (contacts[1,i]-offsetB)/binSize
 
        if (0 <= a < n) and (0 <= b < m):
          if contacts[1,i] <= contacts[0,i]:
            binMatrix[a,b] += contacts[2,i]
 
          if symm and (a != b):
            if contacts[0,i] <= contacts[1,i]:
              binMatrix[b,a] += contacts[2,i]
 
    else:
      for i in range(nCont):
        a = (contacts[0,i]-offsetA)/binSize
        b = (contacts[1,i]-offsetB)/binSize
 
        if (0 <= a < n) and (0 <= b < m):
          if contacts[0,i] <= contacts[1,i]:
            binMatrix[a,b] += contacts[2,i]
 
          if symm and (a != b):
            if contacts[1,i] <= contacts[0,i]:
              binMatrix[b,a] += contacts[2,i]
  
  elif split_side > 0:
    if transpose:
      for i in range(nCont):
        b = (contacts[0,i]-offsetA)/binSize
        a = (contacts[1,i]-offsetB)/binSize
 
        if (0 <= a < n) and (0 <= b < m):
          if contacts[1,i] >= contacts[0,i]:
            binMatrix[a,b] += contacts[2,i]
 
          if symm and (a != b):
            if contacts[0,i] >= contacts[1,i]:
              binMatrix[b,a] += contacts[2,i]
 
    else:
      for i in range(nCont):
        a = (contacts[0,i]-offsetA)/binSize
        b = (contacts[1,i]-offsetB)/binSize
 
        if (0 <= a < n) and (0 <= b < m):
          if contacts[0,i] >= contacts[1,i]:
            binMatrix[a,b] += contacts[2,i]
 
          if symm and (a != b):
            if contacts[1,i] >= contacts[0,i]:
              binMatrix[b,a] += contacts[2,i]
 
  
  else:
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
def outerProdRelEntropy(ndarray[double, ndim=2] matrix,
                        ndarray[double, ndim=1] pExpectA,
                        ndarray[double, ndim=1] pExpectB):
   
  cdef int i, j, n, m
  cdef double o, b, e, s = 0.0
  
  n = len(matrix)
  m = len(matrix[0])
  
  if len(pExpectA) != n:
    raise Exception('pExpectA array length (%d) does not match matrix rows (%d)' % (len(pExpectA), n))

  if len(pExpectB) != m:
    raise Exception('pExpectB array length (%d) does not match matrix columns (%d)' % (len(pExpectB), m))
  
  for i in range(n):
    for j in range(m):
      s += matrix[i,j]
  
  if s == 0.0:
    return s, matrix
  
  b = 1.0 / s
  
  for i in range(n):
    for j in range(m):
      if matrix[i,j] != 0.0:
        o = matrix[i,j] / s
        e = pExpectA[i]*pExpectA[j]
        
        if e < b:
          e = b
        
        matrix[i,j] = o * log(o/e)
  
  return s, matrix
    

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

    
