from . import Image

from math import ceil, log

from numpy import cov, linalg, array, uint8, empty, sqrt, exp, eye, dot
from numpy import ones, vstack, cross, mgrid, zeros, concatenate, sin, cos


TAU = 2.0 * 3.14159265358979323846

ROOT_3OVER2 = 0.75**0.5

HEX_GRID = array( [(0, 0), (1, 0), (0.5, ROOT_3OVER2), (-0.5, ROOT_3OVER2),
                   (-1, 0), (-0.5, -ROOT_3OVER2), (0.5, -ROOT_3OVER2),
                   (-1, 2*ROOT_3OVER2), (0,  2*ROOT_3OVER2), ( 1,  2*ROOT_3OVER2),
                   ( 1.5,  ROOT_3OVER2), ( 2, 0), ( 1.5, -ROOT_3OVER2),
                   (1, -2*ROOT_3OVER2), (0, -2*ROOT_3OVER2), (-1, -2*ROOT_3OVER2),
                   (-1.5, -ROOT_3OVER2), (-2, 0), (-1.5,  ROOT_3OVER2)] )
                   
                      
def indexToHilbert3D(index):  
  
  def grey(start, end, i):
    g = i^(i/2) * ((start^end) * 2)
    return ( (g|g/8) & 7 ) ^ start
  
  nChunks = int( ceil(log(index+1,8)) ) 
  nChunks = max(1, nChunks)
  iChunks = [0]*nChunks
  
  for j in range(nChunks-1, -1, -1):
    iChunks[j] = index % 8
    index /= 8
    
  start = 0
  end = 2**( (-nChunks-1) % 3 )
  chunks = [0]*nChunks
  
  for j, i in enumerate(iChunks):
    chunks[j] = grey(start, end, i)
    start = grey(start, end, max(0, (i-1) & ~1 ))
    end   = grey(start, end, min(7, (i+1) |  1 ))
  
  srcs = list(chunks)
  coord = [0,0,0]
  
  for j in (2,1,0):
    pos = 0

    for k in range(nChunks):
      pos = pos * 2 + srcs[k] % 2
      srcs[k] /= 2
        
    coord[j] = pos
      
  return coord  


def getRotationMatrix(axis, angle):
    
    axis = array(axis)
    axis /= sqrt((axis*axis).sum())
    
    x, y, z = axis
    c = cos(angle)
    d = 1-c
    s = sin(angle)
    R = array([[c+d*x*x,   d*x*y-s*z, d*x*z+s*y],
               [d*y*x+s*z, c+d*y*y,   d*y*z-s*x],
               [d*z*x-s*y, d*z*y+s*x, c+d*z*z  ]])
               
    return R
    
    
def getCoiledCoilCoords(baseCoord, numPoints, dAngle1, dAngle2, rad1, rad2, contact):

  indices = array(range(-1, numPoints+1))
  angles1 = indices * dAngle1 % TAU
  angles2 = indices * dAngle2 % TAU
  
  dx = rad1 * cos(angles1)
  dy = rad1 * sin(angles1)
  dz = (2.0 * rad2 + contact) * (indices * dAngle1) / TAU
  coords1 = vstack([dx, dy, dz]).T

  coords = coords1[1:-1]
  tangents = coords1[2:] - coords1[:-2]
  
  spokes = array(coords)
  spokes[:,:2] = (0.0, 0.0)
  
  spokes = coords - spokes
  
  perps = cross(tangents, spokes, axis=1)
  vLens = sqrt((perps*perps).sum(axis=1))
  
  for i in range(3):
    perps[:,i] /= vLens
  
  perps *= rad2
  
  for i, coord in enumerate(coords):
    rMat = getRotationMatrix(tangents[i], angles2[i])
    coords[i] += dot(perps[i], rMat)
  
  coords += array(baseCoord)

  return coords
  

def getArcChromoCoords(seqPos, seqOffset, totalSeqLen, scale):

  genomeFracs = (seqPos + seqOffset) / float(totalSeqLen)
  angles = genomeFracs * TAU     
  coords = zeros((len(seqPos), 3), float)
  coords[:,0] = scale * cos(angles)
  coords[:,1] = scale * sin(angles)
 
  angle = angles[int(len(angles)/2)]
 
  x = 1.04 * scale * cos(angle)
  y = 1.04 * scale * sin(angle)
  centre = array([x, y, 0.0])
 
  return coords, centre
 
   
def getCircleChromoCoords(seqPos, rank, nChromos, maxLen, scale):
   
  angle = TAU * float(rank)/nChromos
  r = 2.0 * nChromos / TAU
 
  cx = r * cos(angle)
  cy = r * sin(angle)
 
  seqMax = seqPos[-1]
  seqMin = seqPos[0]
  f = sqrt(seqMax/float(maxLen))

  angles = TAU * (seqPos-seqMin)/seqMax
  coords = zeros((len(seqPos), 3), float)
  coords[:,0] = cx + 0.5 * f * cos(angles)
  coords[:,1] = cy + 0.5 * f * sin(angles)

  centre = array([cx, cy, 0.0])
 
  return coords, centre
 
 
def getLinearChromoCoords(seqScale, seqPos, rank, nChromos,  scale):
  
  offset = seqPos[0]
  x = scale * (seqPos-offset) / seqScale
  x -= scale / 2.0
 
  y = scale * ((rank / float(nChromos)) - 0.5)
  z = 0.0
  
  coords = empty((len(seqPos), 3), float)
  coords[:,0] = x
  coords[:,1] = y
  coords[:,2] = z
  
  centre = array([x[0], y, z])
  
  return coords, centre


def gaussMatrix(r=2, sigma=1.4):

  x, y, z = mgrid[-r:r+1, -r:r+1, -r:r+1]

  s2 = 2.0 * sigma * sigma
  
  x2 = x * x / s2
  y2 = y * y / s2
  z2 = z * z / s2
  
  matrix = exp( -(x2 + y2 + z2))
  matrix /= matrix.sum()
  
  return matrix
  
def alignCoords(coordsA, coordsB, weights=None):

  n = len(coordsA)
  if weights is None:
    weights = ones(n)

  rMat = dot(coordsB.transpose()*weights, coordsA)

  rMat1, scales, rMat2 = linalg.svd(rMat)

  sign = linalg.det(rMat1) * linalg.det(rMat2)

  if sign < 0:
    rMat1[:,2] *= -1
  
  rotation = dot(rMat1, rMat2)
  
  coordsB = dot(coordsB, rotation)
  
  return rotation, coordsB
  
def calcRmsds(refCoords, allCoords, weights):

  rmsds = []
  totalWeight = sum(weights)
  totalDeltas = zeros(refCoords.shape)
  
  for coords in allCoords:
    delta = coords-refCoords
    delta *= delta
    totalDeltas += delta
    rmsds.append( sqrt(sum(weights*delta.sum(axis=1))/totalWeight) )
  
  atomRmsds = sqrt(totalDeltas.sum(axis=1)/len(rmsds))

  return rmsds, atomRmsds
  
def centerCoords(coords, weights):

  wCoords = coords.transpose() * weights

  xyzTotals = wCoords.sum(axis=1)

  center = xyzTotals/sum(weights)

  coords -= center
  
  return coords, center

def superimposeCoordPair(coords1, coords2, threshold=1.0, extraRefine=True):

  n = len(coords1)
  weights = ones(n, float)
  
  coords1, cen1 = centerCoords(coords1, weights)
  coords2, cen2 = centerCoords(coords2, weights)

  rotationA, coords2a = alignCoords(coords1, coords2, weights)
  rmsdA, atomRmsdsA = calcRmsds(coords1, [coords2a], weights)

  rotationB, coords2b = alignCoords(coords1, -coords2, weights)
  rmsdB, atomRmsdsB = calcRmsds(coords1, [coords2b], weights)
  
  if rmsdA[0] < rmsdB[0]:
    coords2 = coords2a
    atomRmsds = atomRmsdsA
    rotation = rotationA
    
  else: # Mirror is best
    coords2 = coords2b
    atomRmsds = atomRmsdsB
    rotation = dot(rotationB, -eye(3))
        
  # Refine
  weightScale = atomRmsds/threshold
  weightsExp = ones(n, float)
  weightsExp *= exp(-weightScale*weightScale)

  coords1, cen1 = centerCoords(coords1, weightsExp)
  coords2, cen2 = centerCoords(coords2, weightsExp)

  rotation2, coords2 = alignCoords(coords1, coords2, weightsExp)
  rotation = dot(rotation, rotation2)
  
  if extraRefine:
    rmsd, atomRmsds = calcRmsds(coords1, [coords2], weights) # Not exp
    
    weightScale = atomRmsds/threshold
    weightsExp = ones(n, float)
    weightsExp *= exp(-weightScale*weightScale)
 
    coords1, cen1 = centerCoords(coords1, weightsExp)
    coords2, cen2 = centerCoords(coords2, weightsExp)
  
    rotation2, coords2 = alignCoords(coords1, coords2, weightsExp)
    rotation = dot(rotation, rotation2)
   
  rmsds, atomRmsds = calcRmsds(coords1, [coords2], weights) # Not exp
  
  return coords2, rotation, rmsds[0], atomRmsds
  
  
def superimposeCoordArray(ensemble, threshold=1.0, extraRefine=True):

  # Align to first

  refCoords = ensemble[0]
  rotations = [eye(3)]
  
  for i, coords in enumerate(ensemble[1:], 1):
    coords2, rotation, rmsd, null = superimposeCoordPair(refCoords, coords, threshold, extraRefine)
    ensemble[i] = coords2
    rotations.append(rotation)

  # Align to mean

  ensemble = array(ensemble)
  refCoords = ensemble.mean(axis=0)
  rmsdList = []
  atomRmsds = []
  
  for i, coords in enumerate(ensemble):
    coords2, rotation, rmsd, atomRmsds2 = superimposeCoordPair(refCoords, coords, threshold, extraRefine)
    ensemble[i] = coords2
    rotations[i] = dot(rotations[i], rotation)
    rmsdList.append(rmsd)
    atomRmsds.append(atomRmsds2)
  
  atomRmsds = array(atomRmsds).mean(axis=0) # Could do RMS
  
  return ensemble, rotations, array(rmsdList), atomRmsds


def coordsToDepths(coords, nPoints=100, sigma=0.8, mapThreshold=0.0):

  # Redo in Cython

  nCoords = len(coords)
  dists = empty(nCoords)
  first = coords.min(axis=0)
  last = coords.max(axis=0)
  extent = last - first
  middle = (last + first) / 2.0
  binSize = 2.5
  m = int(5*sigma / binSize) + 1
  n2 = nPoints * 2
  center = array([int(n2//2.0)] * 3)
  
  indexDict = {}
  influence = gaussMatrix(m, sigma)
  voxelArray = zeros((n2, n2 ,n2), float)
  
  for c, coord in enumerate(coords):
    coordRel = coord - middle
    i, j, k = array(coordRel // binSize, int) + center
 
    voxelArray[i-m:i+m+1,
               j-m:j+m+1,
               k-m:k+m+1] += influence
 
    indexDict[c] = array((i, j, k))
 
  border = []
  boolArray1 = voxelArray > 0
  boolArray2 = voxelArray < influence[m/2,m/2,1]
  nonZeros = vstack((boolArray1 & boolArray2).nonzero()).T

  for i, j, k in nonZeros:
 
    for dx in (-1,0,1):
      i2 = i + dx
      for dy in (-1,0,1):
        j2 = j + dy
        for dz in (-1,0,1):
          k2 = k + dz
          value = voxelArray[i2,j2,k2]
 
          if not value:
            border.append((i,j,k))
            break
 
        else:
          continue
        break
 
      else:
        continue
      break
  
  border = array(border, float)
  distDict = {}

  for c in indexDict:
    indices = indexDict[c]
    i, j, k = indices
    key = (i, j, k)
    
    if key in distDict:
      minSq =  distDict[key]
      
    else:  
      diff = border-indices
      sumSq = (diff*diff).sum(axis=1)
      minDist = sumSq.min()
      distDict[key] = minDist
      
    dists[c] = minDist
    
    #xyz = (indices - center) * binSize
    #xyz += middle
    #surface[c].append(xyz.tolist())
  
  dArray = sqrt(dists) * binSize
  dArray -= dArray.min()
  
  return dArray


def coordsToVoxels(coords, colors, binSize, sigma=1.0, gamma=1.0):
  
  # Redo in Cython
  # Better algorithm?
  # Non-cube bbox
  
  first = coords.min(axis=0)
  last = coords.max(axis=0)
  extent = last - first
  middle = (last + first) / 2.0
  width = extent.max() * 1.4
  
  n = int(width/binSize)
  m = int(5*sigma / binSize) + 2
  n2 = n + m + m
  c = array([int(n2//2.0)] * 3)
  
  influence = gaussMatrix(m, sigma)
  m2 = len(influence)
  influence3 = influence.reshape(m2,m2,m2,1)
  influence3 = concatenate([influence3] * 3, axis=3)
  voxelArray = zeros((n2, n2 ,n2), float)
  colorArray = zeros((n2, n2, n2, 3), float)
  
  cubeDict = {}
  nC = len(coords)
  for x, coord in enumerate(coords):
    coordRel = coord - middle
    i, j, k = array(coordRel // binSize, int) + c
    key = (i,j,k)
    
    if key in cubeDict:
      cubeDict[key].append(colors[x])
    else:
      cubeDict[key] = [colors[x],]

  for key in cubeDict:
    i, j, k = key
    colors = array(cubeDict[key]).mean(axis=0)
  
    voxelArray[i-m:i+m+1,
               j-m:j+m+1,
               k-m:k+m+1] += influence
    
    colorArray[i-m:i+m+1,
               j-m:j+m+1,
               k-m:k+m+1] += influence3 * colors

  divisor = array(voxelArray).reshape(n2, n2, n2, 1)
  divisor[(divisor == 0.0).nonzero()] = 1.0
  colorArray /= divisor

  gammaArray = gamma * ones((n2, n2, n2, 1), float)   
  colorArray = concatenate([colorArray, gammaArray], axis=3)

  return voxelArray, colorArray





