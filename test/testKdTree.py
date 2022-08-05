import sys, os, time, random
from os.path import dirname, abspath, exists
from scipy.spatial import cKDTree

# Setup import path

thisDir = dirname(abspath(__file__))
sys.path.remove(thisDir)

nucDir = dirname(thisDir)
sys.path.append(nucDir)

# 
from cUtil.kdTree import KdTree
from numpy import array, dot, random, sqrt, empty
coords = random.uniform(0.0, 10.0, (10000, 3))
myPoints = random.uniform(0.0, 10.0, (1000, 3))

print "Running..."

t0 = time.time()
t = cKDTree(coords)
dists, idx = t.query(myPoints, 2)
nearest = coords[idx]
t1 = time.time()
print '%.3f' % (t1-t0)
#print '%.3f %.3f' % (t1-t0, t2-t1)
known = []
for pt in myPoints:
  d2 = empty(len(coords))
  for i, v in enumerate(coords):
    d = v-pt
    d2[i] = dot(d,d)

  j = d2.argmin()
  known.append(j)
t2 = time.time()

for i, k in enumerate(known):
  delta = coords[k]-nearest[i]
  
  pt = myPoints[i]
  delta1 = pt-coords[k]
  delta2 = pt-nearest[i]
  
  diff = sqrt(dot(delta, delta))
  dist1 = sqrt(dot(delta1, delta1))
  dist2 = sqrt(dot(delta2, delta2))
  
  if diff == 0.0:
    same = 'OK'
    
  else:
    same =  '%4.2f' %diff
  
  print '%4.2f %4.2f %4.4s' % (dist1, dist2, same), coords[k], nearest[i]

print "Timing. KdTree:%.3f  Orig:%.3f" % (t1-t0, t2-t1)
"""
"""
from cUtil.apiUtil import calcSurface, calcCoordDepth

print calcSurface(coords, 0.2, 100)

print calcCoordDepth(coords, 0.2, 100)
