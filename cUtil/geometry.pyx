from libc.math cimport sqrt, sin, cos, exp

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

