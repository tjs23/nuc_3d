from libc.math cimport sqrt, acos, sin, cos
from numpy cimport ndarray
import cython

from numpy import array, empty
import numpy

TAU = 2.0 * 3.14159265358979323846

cdef extern from "GL/gl.h":

  ctypedef double GLdouble
  ctypedef int GLint

  cdef void glTranslated(GLdouble, GLdouble, GLdouble)
  cdef void glRotated(GLdouble, GLdouble, GLdouble, GLdouble)
  cdef void glColor4d(GLdouble, GLdouble, GLdouble, GLdouble)
  cdef void glVertex3d(GLdouble, GLdouble, GLdouble)
  cdef void glNormal3d(GLdouble, GLdouble, GLdouble)
    
    
cdef extern from "GL/glu.h":
  
  ctypedef struct GLUquadric:
    pass
  ctypedef GLUquadric GLUquadricObj
  
  cdef int GLU_SMOOTH, GL_TRIANGLES, GL_TRIANGLE_STRIP, GL_LINE_STRIP
  cdef GLUquadricObj* gluNewQuadric( )
  cdef void gluSphere( GLUquadricObj*, GLdouble, GLint, GLint )
  cdef void gluCylinder( GLUquadricObj*, GLdouble, GLdouble, GLdouble, GLint, GLint )
  cdef void gluQuadricNormals( GLUquadricObj*, GLint)
  cdef void glPushMatrix()
  cdef void glPopMatrix()
  cdef void glBegin(GLint)
  cdef void glEnd()
      
      
@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def createArc(double startAngle, double dAngle, int nAngles,
              double radius, color):
  
  cdef int i
  cdef double angle, x, y, z = 0.0
  cdef double r, g, b, a
  
  r = color[0]
  g = color[1]
  b = color[2]
  a = color[3]
  
  glBegin(GL_LINE_STRIP)
  glColor4d(r, g, b, a)
  
  angle = startAngle
  
  for i in range(nAngles):
    x = radius * cos(angle)
    y = radius * sin(angle)
    glVertex3d(x, y, z)
    angle += dAngle
  
  glEnd()
