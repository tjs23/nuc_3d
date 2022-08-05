from numpy cimport ndarray
from numpy import empty, array, zeros, int32

from libc.math cimport sqrt, sin, cos, exp, ceil, floor

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
    
   
cdef void copy3d(double v1[3], double v2[3]):
  
  v1[0] = v2[0]
  v1[1] = v2[1]
  v1[2] = v2[2]

cdef void sum3d(double v1[3], double v2[3], double v3[3]):
  
  v1[0] = v2[0] + v3[0]
  v1[1] = v2[1] + v3[1]
  v1[2] = v2[2] + v3[2]


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


#   0                                 1                               2
#   |                |                |                |              |
#   0           s    1                2           e    3              4
#               s                                 e

cdef double TAU = 2.0 * 3.14159265358979323846

def getSymbolCoords(ndarray[int, ndim=1] pos,
                    ndarray[double, ndim=2] coords,
                    ndarray[int, ndim=2] regions,
                    ndarray[double, ndim=2] values,
                    double dAngle=TAU/6.0, int nPoints=6,
                    int lineSmooth=0, double scale=1.0,
                    double minVal=0.0, double minRadius=0.1,
                    double maxRadius=10.0):
  
  cdef int i, i2, j, k, k2
  cdef int aPos, bPos, cPos, delta, jOrd
  cdef int n = len(pos)
  cdef int m = len(regions)
  cdef int nCoords = len(coords)
  cdef int step = 2**lineSmooth
  cdef int stride = 2*nPoints
  
  cdef double fstep = float(step)  
  cdef double aFrac, bFrac
  cdef double vec1[3], vec2[3], vec3[3]
  cdef double centre[3], perp[3], axis[3]
  cdef double rMat[3][3]
  
  cdef ndarray[int, ndim=1] posOrder = array(pos.argsort(), int32)
  cdef ndarray[int, ndim=1] regOrder = array(regions[:,0].argsort(), int32)
  
  cdef ndarray[int, ndim=1] idxStart = zeros(m, int32)
  cdef ndarray[int, ndim=1] idxEnd   = zeros(m, int32)
  
  cdef ndarray[double, ndim=1] fracStart = empty(m, float)
  cdef ndarray[double, ndim=1] fracEnd   = empty(m, float)
  cdef ndarray[double, ndim=2] symbolCoords = empty((m*4*nPoints, 3), float)
  
  j = 0
  for i in range(1, n):
    aPos = pos[posOrder[i-1]]
    bPos = pos[posOrder[i]]
    delta = bPos-aPos
    jOrd = regOrder[j]
    
    while regions[jOrd,0] < bPos:
      
      if values[jOrd,1] < minVal:
        j += 1
        if j >= m:
          break
 
        jOrd = regOrder[j]
        continue
      
      k = i
      while (k < n) and (regions[jOrd,1] < pos[posOrder[k]]):
        k += 1
      
      k -= 1
      cPos = pos[posOrder[k-1]]
      idxStart[j] = i-1
      fracStart[j] = (regions[jOrd,0]-aPos) / <double>delta
      
      idxEnd[j]   = k-1
      fracEnd[j] = (regions[jOrd,1]-cPos) / (pos[posOrder[k]]-cPos)
      
      j += 1
      if j >= m:
        break
      
      jOrd = regOrder[j]
      
    
    if j >= m:
      break
  
  i2 = 0  
  for i in range(m):
    # Starts
    
    aFrac = fracStart[i] * fstep
    j = idxStart[i] * step + <int>floor(aFrac)    
    bFrac = aFrac - floor(aFrac)
    aFrac = 1.0 - bFrac
    
    if j >= nCoords-1:
      continue
      
    elif j < 0:
      continue
    
    vec1[0] = coords[j,0]
    vec1[1] = coords[j,1]
    vec1[2] = coords[j,2]
    
    vec2[0] = coords[j+1,0]
    vec2[1] = coords[j+1,1]
    vec2[2] = coords[j+1,2]
    
    centre[0] = (aFrac * vec1[0]) + (bFrac * vec2[0])
    centre[1] = (aFrac * vec1[1]) + (bFrac * vec2[1])
    centre[2] = (aFrac * vec1[2]) + (bFrac * vec2[2])
    
    diff3d(axis, vec2, vec1)
    unit3d(axis)
    
    copy3d(perp, axis)
    perp3d(perp)
    unit3d(perp)
    scale3d(perp, min(maxRadius, max(values[regOrder[i],1] * scale, minRadius)))
    rotMat3d(rMat, axis, dAngle)
           
    for k in range(nPoints):
      sum3d(vec3, centre, perp)
      
      k2 = k*2
      symbolCoords[i2+k2,0] = vec3[0]
      symbolCoords[i2+k2,1] = vec3[1]
      symbolCoords[i2+k2,2] = vec3[2]
    
      matMultVec3d(perp, perp, rMat)
      sum3d(vec3, centre, perp)
          
      k2 += 1
      symbolCoords[i2+k2,0] = vec3[0]
      symbolCoords[i2+k2,1] = vec3[1]
      symbolCoords[i2+k2,2] = vec3[2]
    
    i2 += stride

    
    # Ends
    
    if idxEnd[i] == idxStart[i]:
      continue
      
    aFrac = fracEnd[i] * fstep
    j = idxEnd[i] * step + <int>floor(aFrac)    
    
    bFrac = aFrac - floor(aFrac)
    aFrac = 1.0 - bFrac
    
    if j >= nCoords-1:
      continue
      
    elif j < 0:
      continue
    
    vec1[0] = coords[j,0]
    vec1[1] = coords[j,1]
    vec1[2] = coords[j,2]
    
    vec2[0] = coords[j+1,0]
    vec2[1] = coords[j+1,1]
    vec2[2] = coords[j+1,2]
    
    centre[0] = (aFrac * vec1[0]) + (bFrac * vec2[0])
    centre[1] = (aFrac * vec1[1]) + (bFrac * vec2[1])
    centre[2] = (aFrac * vec1[2]) + (bFrac * vec2[2])
    
    diff3d(axis, vec2, vec1)
    unit3d(axis)
    
    copy3d(perp, axis)
    perp3d(perp)
    unit3d(perp)
    scale3d(perp, min(maxRadius, max(values[regOrder[i],1] * scale, minRadius)))
    rotMat3d(rMat, axis, dAngle)
        
    for k in range(nPoints):
      sum3d(vec3, centre, perp)
      
      k2 = k*2
      symbolCoords[i2+k2,0] = vec3[0]
      symbolCoords[i2+k2,1] = vec3[1]
      symbolCoords[i2+k2,2] = vec3[2]
    
      matMultVec3d(perp, perp, rMat)
      sum3d(vec3, centre, perp)
          
      k2 += 1
      symbolCoords[i2+k2,0] = vec3[0]
      symbolCoords[i2+k2,1] = vec3[1]
      symbolCoords[i2+k2,2] = vec3[2]
    
    i2 += stride
  
  return symbolCoords[:i2-stride]
  
  
def addLayerColor(ndarray[double, ndim=2] colors,
                  ndarray[double, ndim=1] signals,
                  ndarray[int, ndim=1] positions,
                  ndarray[int, ndim=2] regions,
                  ndarray[double, ndim=2] values,
                  ndarray[double, ndim=1] sColor,
                  double threshold):

  cdef int i, j, k, k2, n = len(positions)
  
  cdef double div, div2, divFirst, divLast, pStart, pEnd
  cdef double s, s_max, frac, val, val2
  cdef double r, g, b, a

  cdef ndarray[double, ndim=1] signalsNew = zeros(n)
  cdef ndarray[double, ndim=1] divPoints  = empty(n)
  cdef ndarray[int, ndim=1] starts = regions[:,0]
  cdef ndarray[int, ndim=1] ends = regions[:,1]
  cdef ndarray[int, ndim=1] indices
  
  for i in range(n-1):
    divPoints[i] = <double>(positions[i] + positions[i+1])/2.0
  
  divPoints[n-1] = positions[n-1]
  
  r = sColor[0]
  g = sColor[1]
  b = sColor[2]
  a = sColor[3]
  
  divFirst = divPoints[0]
  divLast = divPoints[n-1]
  
  k = 0
  div = divPoints[k]
  
  indices = array(starts.argsort(), int32)
  
  for j in indices:
    val = values[j,1]
    if val < threshold:
      continue
    
    val -= threshold
    
    pStart = starts[j]
    pEnd = ends[j]
    
    if pStart > divLast:
      break

    if pEnd < divFirst:
      continue
    
    while (div < pStart) and (k+1 < n):
      k += 1
      div = divPoints[k]
    
    if pEnd <= div:
      signalsNew[k] += val
      
    else:
      frac = 1.0
      #frac = (div-pStart) / (pEnd-pStart)
      val2 = frac * val
      signalsNew[k] += val2
     
      k2 = k + 1
      
      while k2 < n:
        div2 = divPoints[k2]
        
        if pEnd < div2:
          frac = 1.0
          #frac = (pEnd-divPoints[k2-1]) / (pEnd-pStart)
          val2 = frac * val
          signalsNew[k2] += val2
          break
        
        else:
          frac = 1.0
          #frac = (div2-divPoints[k2-1]) / (pEnd-pStart)
          val2 = frac * val
          signalsNew[k2] += val2
           
        k2 += 1
      
  s_max = 2.0 * signalsNew.mean()
  
  if not s_max:
    s_max = max(signalsNew.max(), 1.0)
  
  for i in range(n):
    
    if signalsNew[i] > s_max:
      s = 1.0
    else:
      s = signalsNew[i]/s_max
    
    colors[i,0] += s * r
    colors[i,1] += s * g
    colors[i,2] += s * b
    colors[i,3] += s
    signals[i] += min(1.0, s)
        
  return s_max
  
  
def regionBinValues(ndarray[int, ndim=2] regions, ndarray[double, ndim=1] values,
                    int binSize=1000, int start=0, int end=-1, double dataMax=0.0,
                    double scale=1.0, double threshold=0.0):
                    
  cdef int i, p1, p2, b1, b2, b3, s, e
  cdef int nBins, n = len(values)
  cdef double f, r, v, vMin, vMax
  
  if len(regions) != n:
    data = (len(regions), n)
    raise Exception('Number of regions (%d) does not match number of values (%d)' % data) 
  
  if end < 0:
    end = binSize * int32(regions.max() / binSize)
  
  s = start/binSize
  e = end/binSize
  nBins = 1+e-s
  
  cdef ndarray[double, ndim=1] hist = zeros(nBins, float)
  
  for i in range(n):
    
    v = values[i]
    if abs(v) < threshold:
      continue
    
    if regions[i,0] > regions[i,1]:
      p1 = regions[i,1] 
      p2 = regions[i,0]
    
    else:
      p1 = regions[i,0]
      p2 = regions[i,1]
    
    if end < p1:
      continue
    
    if start > p2:
      continue  
    
    b1 = p1 / binSize
    b2 = p2 / binSize
    
    if b1 == b2:
      if b1 < s:
        continue
      
      if b1 > e:
        continue
        
      hist[b1-s] += v

    else:
      r = <double> (p2-p1)
    
      for b3 in range(b1, b2+1):
        if b3 < s:
          continue
        
        if b3 >= e:
          break  
      
        if b3 * binSize < p1:
          f = <double> ((b3+1)*binSize - p1) / r 
        
        elif (b3+1) * binSize > p2:
          f = <double> (p2 - b3*binSize) / r 
        
        else:
          f = 1.0
      
        hist[b3-s] += v * f
  
  if dataMax != 0.0:
    vMin = hist[0]
    vMax = hist[0]
    
    for i in range(1, nBins):
      if hist[i] < vMin:
        vMin = hist[i]
      
      elif hist[i] > vMax:
        vMax = hist[i]
    
    vMax = max(abs(vMin), vMax, dataMax)

    if vMax > 0.0:
      for i in range(0, nBins):
        hist[i] = hist[i]/vMax
  
  for i in range(0, nBins):
    hist[i] = hist[i] * scale  
  
  return hist
  
  
def addPixmapHistogram(ndarray[char, ndim=3] pixmap,
                       ndarray[double, ndim=1] values,
                       ndarray[int, ndim=1] colors,
                       asDensity=False, stranded=False):
                    
  cdef int i, j, a, b, h, w, n = len(values)
  cdef int cPos[4], cZero[4], cNeg[4]
  cdef double th, v
  h = pixmap.shape[0]
  w = pixmap.shape[1]
  
  if w != n:
    raise Exception('Pixmap width (%d) does not match number of values (%d)' % (w,n)) 
  
  cPos[0]  = colors[2]
  cPos[1]  = colors[1]
  cPos[2]  = colors[0]
  cPos[3]  = colors[3]
  cZero[0] = colors[6]
  cZero[1] = colors[5]
  cZero[2] = colors[4]
  cZero[3] = colors[7]
  cNeg[0]  = colors[10]
  cNeg[1]  = colors[9]
  cNeg[2]  = colors[8]
  cNeg[3]  = colors[11]

  if asDensity:
    
    if stranded:
      for i in range(n):
        v = min(1.0, max(-1.0, values[i]))
      
        if v == 0.0:
          continue
        
        elif v > 0.0:
          a = min(255, <int>(255.0 * v))
          b = 255 - a
 
          for j in range(0, h/2):
            pixmap[j,i,0] = (a * cPos[0] + b * cZero[0]) / 255
            pixmap[j,i,1] = (a * cPos[1] + b * cZero[1]) / 255
            pixmap[j,i,2] = (a * cPos[2] + b * cZero[2]) / 255
            pixmap[j,i,3] = (a * cPos[3] + b * cZero[3]) / 255
        
        else:
          a = min(255, <int>(-255.0 * v))
          b = 255 - a
 
          for j in range(h/2, h):
            pixmap[j,i,0] = (a * cNeg[0] + b * cZero[0]) / 255
            pixmap[j,i,1] = (a * cNeg[1] + b * cZero[1]) / 255
            pixmap[j,i,2] = (a * cNeg[2] + b * cZero[2]) / 255
            pixmap[j,i,3] = (a * cNeg[3] + b * cZero[3]) / 255
 
    
    else:    
      for i in range(n):
        v = min(1.0, max(-1.0, values[i]))
        a = min(255, <int>abs(255.0 * v))
        b = 255 - a
        
        for j in range(0, h):
          pixmap[j,i,0] = (a * cPos[0] + b * cZero[0]) / 255
          pixmap[j,i,1] = (a * cPos[1] + b * cZero[1]) / 255
          pixmap[j,i,2] = (a * cPos[2] + b * cZero[2]) / 255
          pixmap[j,i,3] = (a * cPos[3] + b * cZero[3]) / 255
   
    
  else:
    
    if stranded:
      a = h/2
      th = <float>a
      
      for i in range(n):
        v = min(1.0, max(-1.0, values[i]))
        
        if values[i] == 0.0:
          continue
          
        elif values[i] > 0.0:
          b = a - <int>ceil(v*th)
          if b < 0:
            b = 0
        
          for j in range(b,a):
            pixmap[j,i,0] = cPos[0]
            pixmap[j,i,1] = cPos[1]
            pixmap[j,i,2] = cPos[2]
            pixmap[j,i,3] = cPos[3]
            
        else:
          b = a + <int>floor(v*th)
          if b > h:
            b = h
        
          for j in range(a,b):
            pixmap[j,i,0] = cNeg[0]
            pixmap[j,i,1] = cNeg[1]
            pixmap[j,i,2] = cNeg[2]
            pixmap[j,i,3] = cNeg[3]
        
    else:
      a = h
      th = <float>a
    
      for i in range(n):
        v = min(1.0, max(-1.0, values[i]))
        b = a - <int>ceil(abs(v)*th)
        if b < 0:
          b = 0
        
        for j in range(b,a):
          pixmap[j,i,0] = cPos[0]
          pixmap[j,i,1] = cPos[1]
          pixmap[j,i,2] = cPos[2]
          pixmap[j,i,3] = cPos[3]
      
  
  return pixmap
  
  

  
