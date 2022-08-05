import numpy, cython
from numpy import empty, zeros, array, int32
from numpy cimport ndarray
from libc.math cimport sqrt, sin, cos, exp, log, acos, abs, pow

cdef void rotMat3d(double m[3][3], double axis[3], double angle):

  cdef double c = cos(angle)
  cdef double s = sin(angle)
  cdef double d = 1.0 - c
  
  unit3d(axis)
  
  m[0][0] = c + d*axis[0]*axis[0]
  m[0][1] = d*axis[0]*axis[1] - s*axis[2]
  m[0][2] = d*axis[0]*axis[2] + s*axis[1]
  
  m[1][0] = d*axis[1]*axis[0] + s*axis[2]
  m[1][1] = c + d*axis[1]*axis[1]
  m[1][2] = d*axis[1]*axis[2] - s*axis[0]
  
  m[2][0] = d*axis[2]*axis[0] - s*axis[1]
  m[2][1] = d*axis[2]*axis[1] + s*axis[0]
  m[2][2] = c + d*axis[2]*axis[2]  


cdef void matMultVec3d(double v1[3], double v2[3], double m[3][3]):

  cdef double x = v2[0]
  cdef double y = v2[1]
  cdef double z = v2[2]
  
  v1[0] = m[0][0]*x + m[0][1]*y + m[0][2]*z
  v1[1] = m[1][0]*x + m[1][1]*y + m[1][2]*z
  v1[2] = m[2][0]*x + m[2][1]*y + m[2][2]*z
  

cdef void cross3d(double v1[3], double v2[3], double v3[3]):

  v1[0] = v2[1]*v3[2] - v2[2]*v3[1]
  v1[1] = v2[2]*v3[0] - v2[0]*v3[2]
  v1[2] = v2[0]*v3[1] - v2[1]*v3[0]


cdef double dot3d(double v1[3], double v2[3]):
  
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
   
   
cdef void copy3d(double v1[3], double v2[3]):
  
  v1[0] = v2[0]
  v1[1] = v2[1]
  v1[2] = v2[2]


cdef int equal3d(double v1[3], double v2[3]):

  #TODO float comparison < EPSILON?  
  if (v1[0] != v2[0]):
    return 0
  if (v1[1] != v2[1]):
    return 0
  if (v1[2] != v2[2]):
    return 0
  return 1


cdef double dist3d(double v1[3]):

  return sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])


cdef void sum3d(double v1[3], double v2[3], double v3[3]):
  
  v1[0] = v2[0] + v3[0]
  v1[1] = v2[1] + v3[1]
  v1[2] = v2[2] + v3[2]


cdef void mean3d(double v1[3], double v2[3], double v3[3]):
  
  v1[0] = 0.5 * (v2[0] + v3[0])
  v1[1] = 0.5 * (v2[1] + v3[1])
  v1[2] = 0.5 * (v2[2] + v3[2])


cdef void diff3d(double v1[3], double v2[3], double v3[3]):
  
  v1[0] = v2[0] - v3[0]
  v1[1] = v2[1] - v3[1]
  v1[2] = v2[2] - v3[2]


cdef void scale3d(double v1[3], double factor):
  
  v1[0] = v1[0] * factor
  v1[1] = v1[1] * factor
  v1[2] = v1[2] * factor
  
  
cdef void unit3d(double v1[3]):
  
  cdef double size = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])
  
  if size != 0.0:
    v1[0] = v1[0] / size
    v1[1] = v1[1] / size
    v1[2] = v1[2] / size
 
cdef void perp3d(double v1[3]):

  cdef double x, y, z
  
  x = v1[0]
  y = v1[1]
  z = v1[2]
  
  if x == 0.0:
    v1[0] = 1.0
    v1[1] = 0.0
    v1[2] = 0.0
  
  elif y == 0.0:
    v1[0] = 0.0
    v1[1] = 1.0
    v1[2] = 0.0
  
  elif z == 0.0:
    v1[0] = 0.0
    v1[1] = 0.0
    v1[2] = 1.0
  
  else:
    v1[0] = y
    v1[1] = -x
    v1[2] = 0.0


cdef double TAU = 2.0 * 3.14159265358979323846
X = -0.525731112119133606
Z = -0.850650808352039932

icoVertexCoords = array([(-X, 0.0, Z), (X, 0.0, Z), (-X, 0.0, -Z), (X, 0.0, -Z),
                         (0.0, Z, X), (0.0, Z, -X), (0.0, -Z, X), (0.0, -Z, -X),
                         (Z, X, 0.0), (-Z, X, 0.0), (Z, -X, 0.0), (-Z, -X, 0.0)], float)

icoTriIndices = array([(0,4,1),  (0,9,4),  (9,5,4),  (4,5,8),  (4,8,1),
                       (8,10,1), (8,3,10), (5,3,8),  (5,2,3),  (2,7,3),
                       (7,10,3), (7,6,10), (7,11,6), (11,0,6), (0,1,6),
                       (6,1,10), (9,0,11), (9,11,2), (9,2,5), (7,2,11)], int)

squareCoords = array([( 0.0, 1.0, 0.0),
                      ( 1.0, 0.0, 0.0),
                      ( 0.0,-1.0, 0.0),
                      (-1.0, 0.0, 0.0)], float)  

def interpolateCircle(coords):
  
  n = len(coords)
  coords2 = empty((n*2, 3))
  
  for i in range(n):
    j = 2*i
    v1 = coords[i]
    v2 = coords[(i+1) % n]
    
    v12 = v1 + v2
    v12 /= numpy.sqrt(numpy.dot(v12, v12))
    
    coords2[j]   = v1
    coords2[j+1] = v12    
  
  return coords2
 
def icosahedron(coords):
  
  for i in range(20):
    for j in range(3):
      k = 3*i + j
      coords[k,0] = icoVertexCoords[icoTriIndices[i,j],0]
      coords[k,1] = icoVertexCoords[icoTriIndices[i,j],1]
      coords[k,2] = icoVertexCoords[icoTriIndices[i,j],2]

def interpolateSphere(coords):

  n = len(coords)
  coords2 = empty((n*4, 3))
  
  for i in range(n/3):
    j = 3*i
    k = 4*j
    v1 = coords[j]
    v2 = coords[j+1]
    v3 = coords[j+2]
    
    v12 = v1+v2
    v23 = v2+v3  
    v31 = v3+v1
  
    coords2[k:k+3]    = [ v1,v31,v12]
    coords2[k+3:k+6]  = [v23, v2,v12]
    coords2[k+6:k+9]  = [v12,v31,v23]
    coords2[k+9:k+12] = [v23,v31, v3]
    
  
  sz = numpy.sqrt((coords2 * coords2).sum(axis=1))
  coords2[:,0] /= sz
  coords2[:,1] /= sz
  coords2[:,2] /= sz
  
  return coords2
    
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def createBallAndStick(ndarray[double, ndim=2] coords,
                       ndarray[double, ndim=2] colors,
                       double radius, double radius2,
                       int sphDetail):
  
  cdef int i, j, k, k2, n, m, nv, nCoords = len(coords)
  cdef double x, y, z, x2, y2, z2 = 0.0
  cdef double dx, dy, dz, bLen, angle, rLim = 4.0 * radius * radius
  cdef double vert[3]
  cdef double mat[3][3]
  cdef double norm[3]
  cdef double axis[3]
  
  cdef ndarray[double, ndim=2] ico = empty((60, 3))
  cdef ndarray[double, ndim=2] circ1
  cdef ndarray[double, ndim=2] circ2

  icosahedron(ico)
  
  for i in range(sphDetail):
    ico = interpolateSphere(-ico)
  
  circ1 = array(squareCoords)
  for i in range(sphDetail):
    circ1 = interpolateCircle(circ1)
  
  n = len(circ1)
  m = n * 6
  cdef ndarray[double, ndim=2] tube = empty((m, 3))
  cdef ndarray[double, ndim=2] tubeNorms = empty((m, 3))
  
  circ2 = circ1 + array([0.0, 0.0, 1.0])
  for i in range(n):
    j = i * 6
    k = (i+1) % n    
    tube[j]   = circ1[i]
    tube[j+1] = circ1[k]
    tube[j+2] = circ2[i]
    tube[j+3] = circ2[i]
    tube[j+4] = circ1[k]
    tube[j+5] = circ2[k]

    tubeNorms[j]   = circ1[i]
    tubeNorms[j+1] = circ1[k]
    tubeNorms[j+2] = circ1[i]
    tubeNorms[j+3] = circ1[i]
    tubeNorms[j+4] = circ1[k]
    tubeNorms[j+5] = circ1[k]

  n = len(ico)
  cdef ndarray[double, ndim=2] icoR = empty((n, 3))
  for j in range(n):
    icoR[j,0] = ico[j,0]*radius
    icoR[j,1] = ico[j,1]*radius
    icoR[j,2] = ico[j,2]*radius  
  
  nv = (n + m) * nCoords # Max possible sphere vertices + tube vertices
  cdef ndarray[double, ndim=1] pVerts = empty(nv*3)
  cdef ndarray[double, ndim=1] pNorms = empty(nv*3)
  cdef ndarray[double, ndim=1] pColors = empty(nv*4)
  
  k = 0
  k2 = 0
  for i in range(nCoords): 
    i2 = i * 3  
    x = coords[i,0]
    y = coords[i,1]
    z = coords[i,2]
    
    for j in range(n):
      pNorms[k]   = ico[j,0]
      pNorms[k+1] = ico[j,1]
      pNorms[k+2] = ico[j,2]
    
      pVerts[k]   = icoR[j,0]+x
      pVerts[k+1] = icoR[j,1]+y
      pVerts[k+2] = icoR[j,2]+z

      pColors[k2]   = colors[i,0]
      pColors[k2+1] = colors[i,1]
      pColors[k2+2] = colors[i,2]
      pColors[k2+3] = colors[i,3]
      
      k += 3
      k2 += 4

    if (radius2 > 0.01) and (i+1 < nCoords):
      j = i+1
      x2 = coords[j,0]
      y2 = coords[j,1]
      z2 = coords[j,2]
      
      dx = x - x2
      dy = y - y2
      dz = z - z2
     
      bLen = dx*dx + dy*dy + dz*dz
      
      if bLen > rLim:
        bLen = sqrt(bLen)
        dx /= bLen
        dy /= bLen
        dz /= bLen
        
        angle = acos(-dz)
        axis[0] = dy
        axis[1] = -dx
        axis[2] = 0.0        
        rotMat3d(mat, axis, angle)
        
        # copy tube, scale x,y to radius, stretch z to length, rotate
        for j in range(m):
          norm[0] = tubeNorms[j,0]
          norm[1] = tubeNorms[j,1]
          norm[2] = tubeNorms[j,2]
          
          vert[0] = tube[j,0] * radius2
          vert[1] = tube[j,1] * radius2
          vert[2] = tube[j,2] * bLen
        
          # rotate
          matMultVec3d(norm, norm, mat)
          matMultVec3d(vert, vert, mat)
                     
          pNorms[k]   = norm[0]
          pNorms[k+1] = norm[1]
          pNorms[k+2] = norm[2]
 
          pVerts[k]   = vert[0]+x
          pVerts[k+1] = vert[1]+y
          pVerts[k+2] = vert[2]+z

          pColors[k2]   = colors[i,0]
          pColors[k2+1] = colors[i,1]
          pColors[k2+2] = colors[i,2]
          pColors[k2+3] = colors[i,3]
          
          k += 3
          k2 += 4         
  
  pVerts = pVerts[:k]
  pNorms = pNorms[:k]
  pColors = pColors[:k2]
  
  return pVerts, pColors, pNorms


def smoothPath(ndarray[double, ndim=2] coords):
  
  cdef int i, j, k = 0
  cdef int n = len(coords)
  cdef int m = 2 * n -1
  cdef double v1[3], v2[3], v3[3]
  cdef double p0[3], p1[3], p2[3], p3[3]
  cdef double s1, s2, s3,
  
  cdef ndarray[double, ndim=2] coordsOut = empty((m, 3))

  for j in range(3):
    coordsOut[0,j] = coords[0,j]
  
  for i in range(1, n-1): 
    
    for j in range(3):
      p0[j] = coords[i,j]
      v1[j] = coords[i,j] - coords[i-1,j]
      v2[j] = coords[i+1,j] - coords[i,j]
    
    s1 = dist3d(v1)
    
    mean3d(v3, v1, v2)
    
    s3 = dist3d(v3)
    if s3 != 0.0:
      scale3d(v3, s1 * 0.5/s3)

    scale3d(v1, 0.5)
    
    diff3d(p1, p0, v3)
    diff3d(p2, p0, v1)

    mean3d(p3, p1, p2)

    k = 2 * i -1
    for j in range(3):
      coordsOut[k,j] = p3[j]
          
    k += 1
    for j in range(3):
      coordsOut[k,j] = coords[i,j]
   
  
  s2 = dist3d(v2)
  s3 = dist3d(v3)
  
  mean3d(v3, v1, v2)
  
  if s3 != 0.0:
    # scale3d(v3, s2 * 0.5/s3)
    scale3d(v3, 0.5)
  
  scale3d(v2, 0.5)
  
  sum3d(p1, p0, v3)
  sum3d(p2, p0, v2)
  
  mean3d(p3, p1, p2)
    
  k = m-2
  for j in range(3):
    coordsOut[k,j] = p3[j]
    
  k = m-1
  for j in range(3):
    coordsOut[k,j] = coords[n-1,j]
      
  return coordsOut
  

# Below needs testing
def meshNormals(ndarray[double, ndim=1] coords):
  
  cdef int i, j, k, l, d
  cdef int n = len(coords)
  cdef int idx[10]
  cdef int m = n/3
  cdef double va[3], vb[3], vc[3], norm[3]
  cdef double c, pa[4], pb[4]
  cdef double th=0.01
  cdef ndarray[double, ndim=1] normals = empty((n,))
  cdef ndarray[double, ndim=2] sortCoords = empty((m,4))
        
  for i in range(0, n, 9):
    
    for j in range(3):
      va[j] = coords[i+j]
      vb[j] = coords[i+3+j]
      vc[j] = coords[i+6+j]

    diff3d(va, vb, va)
    diff3d(vb, vc, va)
    cross3d(norm, va, vb)
    unit3d(norm)
    
    for j in range(0, 9, 3):
      normals[i+j] = norm[0]
      normals[i+j+1] = norm[1]
      normals[i+j+2] = norm[2]
      
  for i in range(0, n, 3):
    j = i/3
    sortCoords[j,0] = coords[i]
    sortCoords[j,1] = coords[i+1]
    sortCoords[j,2] = coords[i+2]
    sortCoords[j,3] = <double>i
  
  # TBD yuck
  sortCoords2 = array(sortCoords, dtype=[('x',float), ('y',float), ('z',float), ('w',float)])
  sortCoords2.sort(order=('x','y','z'))
  sortCoords = array(sortCoords2, float)
  
  for j in range(4):
    pa[j] = sortCoords[0,j]
  
  i = <int>pa[4]
  idx[0] = i
  for j in range(3):
    norm[j] = normals[i+j]
  
  d = 1
  c = 1.0
  for i in range(1, m):
    for j in range(4):
      pb[j] = sortCoords[i,j]
    
    if (pb[0]-pa[0] < th) and (pb[1]-pa[1] < th) and (pb[2]-pa[2] < th):
      k = <int>pb[4]
      idx[d] = k
      for j in range(3):
        norm[j] += normals[k+j]
      
      d += 1
      c += 1.0
      
    else:
      
      for j in range(3):
        norm[j] /= c
      
      for j in range(d):
        l = <int>idx[j]
        for k in range(3):
          normals[l+k] = norm[k]
      
      for j in range(4):
        pa[j] = pb[j]
      
      k = <int>pa[4]
      idx[0] = k
      for j in range(3):
        norm[j] = normals[k+j]
        
      d = 1
      c = 1.0
  
  # Last 
     
  for j in range(3):
    norm[j] /= c
  
  for j in range(d):
    for k in range(3):
      normals[idx[d]+k] = norm[k]
      
  return normals
  
  
#For an array of coordinates (and colours for each coordinate),
#calculate the tube (of given radius) surface coordinates,
#at numStep points per perimeter resolution.
#Return the coordinates, the surface normals and colours for each coordinate.
#At the ends, to close off the tube, the first/last coordinates are used as an extra coordinate with zero radius.
def createTube(ndarray[double, ndim=2] coords,
                       ndarray[double, ndim=2] colors,
                       double radius, int numSteps):
  
  cdef int i, i2, j, j2, j3, j4, k, nCoords = len(coords)
  cdef int n = nCoords*numSteps
  
  cdef double v0[3], v1[3], v2[3], v3[3]
  cdef double axis[3], perp[3]
  cdef double this[3], prev[3], next[3]
  cdef double o0[3], o1[3], o2[3]
  cdef double c0[4], c1[4], c2[4]
  cdef double n0[3], n1[3], n2[3], n3[3]
  cdef double proj, proj2, angle = TAU/<double>numSteps
  cdef double rMat[3][3]
  
  cdef ndarray[double, ndim=2] spokes  = empty((n, 3))
  cdef ndarray[double, ndim=1] tubeCoords  = empty(n*6*3)
  cdef ndarray[double, ndim=1] tubeNormals = empty(n*6*3)
  cdef ndarray[double, ndim=1] tubeColors  = empty(n*6*4)

  # To calculate the tube coordinates, let's build numSteps spokes around each coordinate,
  # (numSteps * n_coord) in total

  for i in range(nCoords):   

    # calculate directional vectors of chain v1, v2, v3:
    # v1_i = x_i - x_{i-1}    except for the first one: v_1_0 = x_1 - x_0
    # v2_i = x_{i+1} - x_i    except for the last one: v2_N = x_{N} - x_{N-1}
    
    for k from 0 <= k < 3:
      this[k] = coords[i,k] #this
      #next
      if i != nCoords-1:
        next[k] = coords[i+1,k]
      else:
        next[k] = coords[i-1,k] #for last, use prev
      #prev
      if i != 0:
        prev[k] = coords[i-1,k]
      else:
        prev[k] = coords[i+1,k] #for first, use next

    diff3d(v1,this,prev)
    diff3d(v2,next,this)

    # calculate the rotation matrix to create points on the tube surface
    # with x_{i-1} -> x_{i+1} direction axis
    # using a (3 PI / n_points) angle
    
    sum3d(axis, v1, v2)
    unit3d(axis) 

    rotMat3d(rMat, axis, angle)
      
    # calculate normal of the v1-v2 plane with a length of the tube radius
    if equal3d(v1, v2): #avoid zero normal vector
      if equal3d(v1,[0,0,1]): #use (arbitrarily) z axis if we can
        cross3d(perp, v1, [0,1,0]) #use y axis instead
      else:
        cross3d(perp, v1, [0,0,1])
    else:
      cross3d(perp, v1, v2)

    unit3d(perp)
    scale3d(perp, radius)

    #build the spokes around each coordinate, (numSteps * n_coord) in total

    if i == 0: # first

      #add spokes starting with an arbitrary one (perp)
      i2 = 0
      for k from 0 <= k < 3:
        spokes[i2,k] = perp[k]
 
      copy3d(n1, perp)
      for j from 1 <= j < numSteps:
        matMultVec3d(n1, n1, rMat)
 
        for k from 0 <= k < 3:
          spokes[i2+j,k] = n1[k]
    
    else: # all other coords

      #find the spoke for this coordinate which is most aligned with the previous coordinate's first spoke
      i2 = (i-1) * numSteps
      for k from 0 <= k < 3:
        v3[k] = spokes[i2,k]
        
      proj = dot3d(perp, v3)
      copy3d(v2, perp)
      
      copy3d(v1, perp)
      for j from 1 <= j < numSteps:
        matMultVec3d(v1, v1, rMat)
        proj2 = dot3d(v1, v3)
        
        if proj2 > proj:
          proj = proj2
          copy3d(v2, v1)

      #add spokes starting with this most aligned spoke
      i2 = i * numSteps
      for k from 0 <= k < 3:
        spokes[i2,k] = v2[k]
        
      copy3d(v1, v2)
      for j from 1 <= j < numSteps:
        matMultVec3d(v1, v1, rMat)
 
        for k from 0 <= k < 3:
          spokes[i2+j,k] = v1[k]


  #now that we've got the spokes, let's build the tube surface coordinates
  #for each coordinate, we need this:
  #  this coord this spoke
  #  prev coord this spoke
  #  this coord next spoke (using rewrapping)
  #  this coord this spoke
  #  next coord next spoke
  #  this coord next spoke
         
  for i from 0 <= i < nCoords: # all coordinates

    i2 = i * numSteps

    #coordinates
    for k from 0 <= k < 3:
      o0[k] = coords[i,k]
      o1[k] = coords[i-1,k] if i != 0 else coords[i,k]
      o2[k] = coords[i+1,k] if i != nCoords-1 else coords[i,k]

    #turn spoke vectors into surface coordinates
    for j from 0 <= j < numSteps:
      j2 = i2+j
      j3 = j2 * 18
      j4 = j2 * 24
      
      for k from 0 <= k < 3:
        v0[k] = spokes[j2, k] #this spoke vector on this coord
        if i == 1: #first coord
          v1[k] = 0 #collapse all spoke vectors of the prev coord into the centre
        else:
          v1[k] = spokes[j2-numSteps, k] #this spoke vector on prev coord
      
      #next spoke vectors on this/prev/next coordinates
      if j == numSteps-1: #rewrapping
        for k from 0 <= k < 3:
          v2[k] = spokes[i2, k] #this coord
          if i == nCoords - 1: #last coord
            v3[k] = 0 #collapse all spoke vectors of the next coord into the centre
          else:
            v3[k] = spokes[i2+numSteps, k] #next coord
      
      else:
        for k from 0 <= k < 3:
          v2[k] = spokes[j2+1, k] #next spoke vector on this coord
          if i == nCoords - 1: #last coord
            v3[k] = 0
          else: #collapse all spoke vectors of the next coord into the centre
            v3[k] = spokes[j2+1+numSteps, k]
          
      #for tube surface normals, use spoke vectors
      copy3d(n0, v0)
      if i != 1:
        copy3d(n1, v1)
      else: #first: use the axis vector to close it off
        diff3d(n1,this,next)
      copy3d(n2, v2)
      if i != nCoords - 1:
        copy3d(n3, v3)
      else: #last: use the axis vector to close it off
        diff3d(n1,this,prev)

      unit3d(n0)
      unit3d(n1)
      unit3d(n2)
      unit3d(n3)

      #actual coordinates on the tube surface (coord + spoke vector)
      sum3d(v0, o0, v0) #this spoke vector on this coord
      sum3d(v1, o1, v1) #this spoke vector on prev coord
      sum3d(v2, o0, v2) #next spoke vector on same coord
      sum3d(v3, o2, v3) #same spoke vector on next coord       
      
      for k from 0 <= k < 3:
        
        tubeCoords[j3+k]    = v0[k] #this spoke on this coord
        tubeCoords[j3+3+k]  = v1[k] #this spoke on prev coord
        tubeCoords[j3+6+k]  = v2[k] #next spoke on this coord
        tubeCoords[j3+9+k]  = v0[k] #this spoke on this coord
        tubeCoords[j3+12+k] = v3[k] #next spoke on next coord
        tubeCoords[j3+15+k] = v2[k] #next spoke on this coord

      for k from 0 <= k < 3:
        tubeNormals[j3+6+k]  = n2[k] #this spoke on this coord
        tubeNormals[j3+k]    = n0[k] #this spoke on prev coord
        tubeNormals[j3+3+k]  = n1[k] #next spoke on this coord
        tubeNormals[j3+9+k]  = n0[k] #this spoke on this coord
        tubeNormals[j3+12+k] = n3[k] #next spoke on next coord
        tubeNormals[j3+15+k] = n2[k] #next spoke on this coord
     
      for k from 0 <= k < 4:
        tubeColors[j4+k]    = colors[i,k] #this spoke on this coord
        #for first: to close tube off, use the first coordinate as the previous coord
        tubeColors[j4+4+k]  = colors[i-1,k] if i != 0 else colors[i,k]
        tubeColors[j4+8+k]  = colors[i,k] #next spoke on this coord
        tubeColors[j4+12+k] = colors[i,k] #this spoke on this coord
        #for last: to close tube off, use the last coordinate as the next coord
        tubeColors[j4+16+k] = colors[i+1,k] if i != nCoords -1 else colors[i,k]
        tubeColors[j4+20+k] = colors[i,k] #next spoke on this coord
      
  return tubeCoords, tubeColors, tubeNormals

cdef gaussMatrix(ndarray[double, ndim=3] matrix, double sigma):
  
  cdef int i, j, k, r, n = len(matrix)
  cdef double x, y, z, v, t=0.0
  cdef double s2 = 2.0 * sigma * sigma
  
  r =  (n-1)/2
  
  i = 0
  for x in range(-r, r+1):
    j = 0
    for y in range(-r, r+1):
      k = 0
      for z in range(-r, r+1):
        v = exp( -(x*x + y*y + z*z)/s2)
        matrix[i,j,k] = v
        t += v
        
        k += 1
      j += 1
    i += 1    
  
  
  for i in range(n):
    for j in range(n):
      for k in range(n):
        matrix[i,j,k]  /= t
  
  return matrix


def getContactMapRegions(ndarray[int, ndim=2] regions,
                         ndarray[double, ndim=1] values,
                         seqStartA, seqStartB,
                         pixStartA, pixStartB,
                         pixEndA, pixEndB, binSize):
 
  cdef int i, n = len(regions)
  cdef int px1, px2, py1, py2, m = 0
  
  cdef int sa = seqStartA
  cdef int sb = seqStartB
  cdef int isa = pixStartA
  cdef int isb = pixStartB
  cdef int iea = pixEndA
  cdef int ieb = pixEndB
 
  cdef double b = binSize
  cdef double ea = sa + <int>(b * <float>(iea-isa))
  cdef double eb = sb + <int>(b * <float>(ieb-isb))
  
  cdef ndarray[int, ndim=3] pixRegions = zeros((n,2,2), int32)
  cdef ndarray[double, ndim=1] pixValues = zeros(n, float)
    
  for i in range(n):
    
    if sa < regions[i,0] < ea:
      if sa < regions[i,1] < ea:
        py1 = isa + <int>(<double>(regions[i,0]-sa) / b)
        py2 = isa + <int>(<double>(regions[i,1]-sa) / b)
        
      else:
        py1 = isa + <int>(<double>(regions[i,0]-sa) / b)
        py2 = iea
    
    elif sa < regions[i,1] < ea:
      py1 = isa
      py2 = isa + <int>(<double>(regions[i,1]-sa) / b)
    
    else:
      continue
    
    if sb < regions[i,0] < eb:
      if sb < regions[i,1] < eb:
        px1 = isb + <int>(<double>(regions[i,0]-sb) / b)
        px2 = isb + <int>(<double>(regions[i,1]-sb) / b)
        
      else:
        px1 = isb + <int>(<double>(regions[i,0]-sb) / b)
        px2 = ieb
    
    elif sb < regions[i,1] < eb:
      px1 = isb
      px2 = isb + <int>(<double>(regions[i,1]-sb) / b)
    
    else:
      continue
    
    if (px2 - px1) < 1:
      px2 = px1 + 1

    if (py2 - py1) < 1:
      py2 = py1 + 1
    
    pixRegions[m,0,0] = px1
    pixRegions[m,0,1] = py1
    pixRegions[m,1,0] = px2
    pixRegions[m,1,1] = py2
    pixValues[m] = values[i]
    
    m += 1
  
  return pixRegions[:m], pixValues[:m]
  

def addContactMapRegions(ndarray[int, ndim=3] pixmap,
                         ndarray[int, ndim=3] regions, # Array of box (start(x,y),end(x,y)) points
                         ndarray[double, ndim=1] values,  # Value for each box
                         int r, int g, int b):
                      
  cdef int i, j, sx, sy, ex, ey, n, w, h, rScaled, gScaled, bScaled
  cdef float rFloat, gFloat, bFloat, f, maxVal = 0.0
  
  n = len(regions)
  w = len(pixmap[0])
  h = len(pixmap)
  
  if len(values) != n:
    data = (n, len(values))
    raise Exception('Number of regions (%d) does not match number of values (%d)' % data)
  
  rFloat = <double>r
  gFloat = <double>g
  bFloat = <double>b

  for i in range(n):
    if abs(values[i]) > maxVal:
      maxVal = abs(values[i])
  
  if maxVal == 0.0:
    maxVal = 1.0
  
  for i in range(n):
    sx = regions[i,0,0] 
    sy = regions[i,0,1] 
    ex = regions[i,1,0]
    ey = regions[i,1,1]
    f = abs(values[i]) / maxVal
    rScaled = <int>(f * rFloat)
    gScaled = <int>(f * gFloat)
    bScaled = <int>(f * bFloat)
      
    for j in range(sx,ex+1):
      if 0 <= j < w:
        if 0 <= sy < h:
          if bScaled > pixmap[sy,j,0]:
            pixmap[sy,j,0] = bScaled
          
          if gScaled > pixmap[sy,j,1]:
            pixmap[sy,j,1] = gScaled
          
          if rScaled > pixmap[sy,j,2]:
            pixmap[sy,j,2] = rScaled
 
        if 0 <= ey < h:
          if bScaled > pixmap[ey,j,0]:
            pixmap[ey,j,0] = bScaled
          
          if gScaled > pixmap[ey,j,1]:
            pixmap[ey,j,1] = gScaled
          
          if rScaled > pixmap[ey,j,2]:
            pixmap[ey,j,2] = rScaled
      
    for j in range(sy,ey+1):
      if 0 <= j < h:
        if 0 <= sx < w:
          if bScaled > pixmap[j,sx,0]:
            pixmap[j,sx,0] = bScaled
          
          if gScaled > pixmap[j,sx,1]:
            pixmap[j,sx,1] = gScaled
          
          if rScaled > pixmap[j,sx,2]:
            pixmap[j,sx,2] = rScaled
 
        if 0 <= ex < w:
          if bScaled > pixmap[j,ex,0]:
            pixmap[j,ex,0] = bScaled
          
          if gScaled > pixmap[j,ex,1]:
            pixmap[j,ex,1] = gScaled
          
          if rScaled > pixmap[j,ex,2]:
            pixmap[j,ex,2] = rScaled
 
  return pixmap


def addContactMapMesh(ndarray[int, ndim=3] pixmap,
                      ndarray[int, ndim=1] xGrid,
                      ndarray[int, ndim=1] yGrid,
                      int r, int g, int b):
                      
  cdef int i, j, n, m
  
  n = len(pixmap)
  m = len(pixmap[0])
  
  for i in yGrid:

    for j in range(m):
      pixmap[i,j,0] = b
      pixmap[i,j,1] = g
      pixmap[i,j,2] = r
      
  for j in xGrid:

    for i in range(n):
      pixmap[i,j,0] = b
      pixmap[i,j,1] = g
      pixmap[i,j,2] = r
  
  return pixmap
      
def addContactMapPoints(ndarray[int, ndim=3] pixmap,
                        ndarray[int, ndim=2] matrix,
                        int r, int g, int b):

  cdef int i, j, n, m, v, r2, g2, b2
  
  r2 = (4*r)/3
  g2 = (4*g)/3
  b2 = (4*b)/3
  
  n = min(len(matrix), len(pixmap))
  m = min(len(matrix[0]), len(pixmap[0]))
  
  for i in range(n):
    for j in range(m):
      if matrix[i,j] != 0:
        pixmap[i,j,0] += b
        pixmap[i,j,1] += g
        pixmap[i,j,2] += r

        if pixmap[i,j,0] > 255:
          pixmap[i,j,0] = 255

        if pixmap[i,j,1] > 255:
          pixmap[i,j,1] = 255

        if pixmap[i,j,2] > 255:
          pixmap[i,j,2] = 255

  """
        for p in range(-4, 5):
          for q in range(-4, 5):
            
            if (p*p + q*q) > 16:
              continue
             
            if i+p < 0:
              continue
            if j+q < 0:
              continue
            if i+p >= n:
              continue
            if j+q >= m:
              continue
 
            pixmap[i+p,j+q,0] += b2
            pixmap[i+p,j+q,1] += g2
            pixmap[i+p,j+q,2] += r2

            if pixmap[i+p,j+q,0] > 255:
              pixmap[i+p,j+q,0] = 255

            if pixmap[i+p,j+q,1] > 255:
              pixmap[i+p,j+q,1] = 255

            if pixmap[i+p,j+q,2] > 255:
              pixmap[i+p,j+q,2] = 255
         

  """
        
  return pixmap
                      
                      
def addContactMapBlocks(ndarray[int, ndim=3] pixmap,
                        ndarray[int, ndim=2] matrix,
                        int r, int g, int b,
                        int nPix, int nMax,
                        double gamma=1.0):

  cdef int i, j, i2, j2, r2, g2, b2, n, m
  cdef double v, vMax
  
  n = min(len(matrix), len(pixmap))
  m = min(len(matrix[0]), len(pixmap[0]))
  
  vMax = <float>nMax
  
  if vMax == 0.0:  # Nothing to add
    return pixmap
  
  vMax = log(vMax)

  for i in range(n):
    for j in range(m):
      if matrix[i,j] != 0:
        v = log(matrix[i,j])/vMax
        if v > 1.0:
          v = 1.0

        v = pow(v, gamma)

        r2 = <int>(r*v)
        g2 = <int>(g*v)
        b2 = <int>(b*v)
      
        for i2 in range(i, i+nPix):
          if i2 == n:
            break
        
          for j2 in range(j, j+nPix):
            if j2 == m:
              break
            
            pixmap[i2,j2,0] = b2
            pixmap[i2,j2,1] = g2
            pixmap[i2,j2,2] = r2
 
            if pixmap[i2,j2,0] > 255:
              pixmap[i2,j2,0] = 255
 
            if pixmap[i2,j2,1] > 255:
              pixmap[i2,j2,1] = 255
 
            if pixmap[i2,j2,2] > 255:
              pixmap[i2,j2,2] = 255
 
  return pixmap


def getFractionPickColors(ndarray[double, ndim=1] fractions):

  cdef int i, n = len(fractions)
  cdef long j, x, r, g, b
  cdef ndarray[double, ndim=2] colors = empty((n, 4), float)
  
  for i in range(n):
    j = <long>(4228250625*fractions[i]) # index wrt max index
    
    r = j / 16581375
    x = j % 16581375
    g = x / 65025
    x = x % 65025
    b = x / 255
    a = x % 255
    
    colors[i,0] = <double>(r)/255.0
    colors[i,1] = <double>(g)/255.0
    colors[i,2] = <double>(b)/255.0
    colors[i,3] = <double>(a)/255.0

  return colors


def coordsToVoxels(ndarray[double, ndim=2] coords,
                   ndarray[double, ndim=2] colors,
                   double binSize, double sigma=1.0):
  
  start = coords.min(axis=0)
  end = coords.max(axis=0)
  extent = end-start
  
  cdef int i, j, k, u, v, w, c
  cdef int n1, n2, n3, c1, c2, c3
  cdef int m, m2, n = len(coords)
  cdef double inf, x, y, z, r, g, b, a
  
  m = 5 * <int>(sigma/binSize)
  m2 = (m-1)/2
  
  n1 = m + m + <int>(extent[0]/binSize)
  n2 = m + m + <int>(extent[1]/binSize)
  n3 = m + m + <int>(extent[2]/binSize)
  
  c1 = n1/2
  c2 = n2/2
  c3 = n3/2
  
  cdef ndarray[double, ndim=1] middle = (start + end) / 2.0
  cdef ndarray[double, ndim=3] influence =  empty((m, m ,m))
  cdef ndarray[double, ndim=3] voxelArray = zeros((n1, n2 ,n3))
  cdef ndarray[double, ndim=4] colorArray = zeros((n1, n2, n3, 4))
  
  gaussMatrix(influence, sigma)
  
  for c in range(n):    
    u = <int>( (coords[c,0]-middle[0]) / binSize ) + c1 - m2
    v = <int>( (coords[c,1]-middle[1]) / binSize ) + c2 - m2
    w = <int>( (coords[c,2]-middle[2]) / binSize ) + c3 - m2
    
    r = colors[c,0]
    g = colors[c,1]
    b = colors[c,2]
    a = colors[c,3]
    
    for i in range(m):
      x = u+i
      
      for j in range(m):
        y = v+j
        
        for k in range(m):
          z = v+k
          
          inf = influence[i,j,k]
          voxelArray[x,y,z] += inf
          
          colorArray[x,y,z,0] += r * inf
          colorArray[x,y,z,1] += g * inf
          colorArray[x,y,z,2] += b * inf
          colorArray[x,y,z,3] += a * inf
  
  for x in range(n1):
    for y in range(n2):
      for z in range(n3):
        inf = voxelArray[x,y,z]
        
        if inf:
          colorArray[x,y,z,0] /= inf
          colorArray[x,y,z,1] /= inf
          colorArray[x,y,z,2] /= inf
          colorArray[x,y,z,3] /= inf

  return voxelArray, colorArray

