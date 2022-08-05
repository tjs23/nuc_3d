from numpy import array, sqrt, dot, random, zeros, ones, arange
from numpy import cross, eye, outer, linalg, float32, vstack, empty, append

from random import shuffle, randint
from os import path
from colorsys import hsv_to_rgb

uniform = random.uniform

from math import exp, sin, cos, floor, atan2, degrees, radians, acos, tan

from PySide import QtCore, QtGui, QtOpenGL
from OpenGL import GL, GLU
from OpenGL import GLUT # Temporary fopr postional labels

Qt = QtCore.Qt
QGL = QtOpenGL.QGL

LeftButton = QtCore.Qt.LeftButton
MiddleButton = QtCore.Qt.MiddleButton
RightButton = QtCore.Qt.RightButton
QPoint = QtCore.QPoint

GL_COLOR_BUFFER_BIT = GL.GL_COLOR_BUFFER_BIT
GL_DEPTH_BUFFER_BIT = GL.GL_DEPTH_BUFFER_BIT
GL_MODELVIEW = GL.GL_MODELVIEW
GL_PROJECTION = GL.GL_PROJECTION
#GLUT_BITMAP_HELVETICA_12 = GLUT.GLUT_BITMAP_HELVETICA_12

glBegin = GL.glBegin
glCallList = GL.glCallList
glClear = GL.glClear
glColor3f = GL.glColor3f
glColor4f = GL.glColor4f
glDisable = GL.glDisable
glEnable = GL.glEnable
glEnd = GL.glEnd
glGetFloatv = GL.glGetFloatv
glLineWidth = GL.glLineWidth
glLoadIdentity = GL.glLoadIdentity
glMatrixMode = GL.glMatrixMode
glMultMatrixf = GL.glMultMatrixf
glNormal3f = GL.glNormal3f
glPopMatrix = GL.glPopMatrix
glPushMatrix = GL.glPushMatrix
glRasterPos2d = GL.glRasterPos2d
glRasterPos3d = GL.glRasterPos3d
glRotatef = GL.glRotatef
glScalef = GL.glScalef
glTranslatef = GL.glTranslatef
gluCylinder = GLU.gluCylinder
gluDeleteQuadric = GLU.gluDeleteQuadric
gluNewQuadric = GLU.gluNewQuadric
gluPerspective = GLU.gluPerspective
gluSphere = GLU.gluSphere
glutBitmapCharacter = GLUT.glutBitmapCharacter
glVertex2d = GL.glVertex2d
glVertex3d = GL.glVertex3d

    
START_ROTATION = array([[ 1.0,  0.0,  0.0,  0.0],
                        [ 0.0,  1.0,  0.0,  0.0],
                        [ 0.0,  0.0,  1.0,  0.0],
                        [ 0.0,  0.0,  0.0,  1.0]],
                        float32)


TAU = 2.0 * 3.14159265358979323846

VERTEX_SHADER_GLSL = """\
#version 120
 
attribute vec4 coord;
varying vec2 texcoord;
 
void main(void) {
  gl_Position = vec4(coord.xy, 0, 1);
  texcoord = coord.zw;
}
"""

FRAGMENT_SHADER_GLSL = """\
#version 120
 
varying vec2 texcoord;
uniform sampler2D tex;
uniform vec4 color;
 
void main(void) {
  gl_FragColor = vec4(1, 1, 1, texture2D(tex, texcoord).a) * color;
}
"""

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
   
class Gl3dPanel(QtOpenGL.QGLWidget):
  
  
  def __init__(self, parent, mainApp, openFunc, **kw):
   
    format = QtOpenGL.QGLFormat(QGL.DepthBuffer | QGL.DirectRendering | QGL.AlphaChannel | QGL.SampleBuffers)
    QtOpenGL.QGLWidget.__init__(self, format, parent)

    self.mainApp = mainApp
    self.openFunc = openFunc
    
    self.bgColor = [1.0, 1.0, 1.0, 1.0]
    #self.glGeometry = [400, 300, 8192, 60] # w, h, farPt, fov
    self.glGeometry = [400, 300, 8192, 45] # w, h, farPt, fov
    self.view = [0.0, 0.0, 1.0] # dx, dy, z scale
    self.refPos = [0.0, 0.0, 1.0] # dx, dy, z scale
    self.glLists = []
    self.clipPlaneFront = 1.0
    self.clipPlaneBack = 1.0
    self.clipLimit = 512.0
    self.rulerPos = None
    self.lighting = set()
    self.pickBuffer = set()
    self.reference = set()
    self.moveRef = False
    self.moveMain = True
    self.prevVec = self.currVec = array([0.0, 0.0, 0.0])
    self.rotation = START_ROTATION
    self.rotation2 = START_ROTATION
    self.projectMatrix = eye(4)
    self.modelMatrix = eye(4)
    self.origin = zeros(3)
 
    self.mouseLeftMoveFunc = self._mouseRotate
    self.mouseMiddleMoveFunc = self._mouseTranslate
    self.mouseRightMoveFunc = self._mouseRotate
    self.translateStep = 10.0
    self.zoomScale = 1.1
    self.setAcceptDrops(True) 

    self.prevPos = QPoint()
    
    
  def dragEnterEvent(self, event):

    if event.mimeData().urls():
      event.acceptProposedAction()
        
    else:
      event.ignore()


  def dragMoveEvent(self, event):
    
    self.dragEnterEvent(event)


  def dragLeaveEvent(self, event):
  
    event.ignore()

    
  def dropEvent(self, event):
     
    if event.mimeData().urls():
      mods = event.keyboardModifiers()
      haveCtrl = mods & QtCore.Qt.ControlModifier
      haveShift = mods & QtCore.Qt.ShiftModifier
      
      if haveShift or haveCtrl:
        replace = True
      else:
        replace = False 
      
      filePaths = [url.path() for url in event.mimeData().urls() if path.isfile(url.path())]
      
      print(filePaths, self.openFunc )
      
      if filePaths and self.openFunc:
        self.openFunc(filePaths, replace)
        event.acceptProposedAction()     
      else:
        event.ignore()
        
    else:
      event.ignore()


  def minimumSizeHint(self):
  
    return QtCore.QSize(50, 50)


  def sizeHint(self):
  
    return QtCore.QSize(400, 300)
     
  def mouseZoom(self, delta):
    
    if delta:
      if delta < 0:
        z = self.zoomScale
      else:
        z = 1.0/self.zoomScale
      
      self.zoom(z)
  
  def zoom(self, z):
    
    if self.moveRef:
      prev = self.refPos[2]
      self.refPos[2] = min(max(0.01, prev/z), 64.0)
      
    if self.moveMain:
      prev = self.view[2]
      self.view[2] = min(max(0.01, prev*z), 64.0)
    
    self.update()
    
  def resetZoom(self):
    
    if self.moveRef:
      self.refPos[2] = 1.0
    
    if self.moveMain:
      self.view[2] = 1.0
      
    self.update()
  
  def resetView(self):
    
    if self.moveRef:
      z =  self.refPos[2]
      self.refPos = [0.0, 0.0, z]
    
    if self.moveMain:
      z =  self.view[2]
      self.view = [0.0, 0.0, z]
    
    self.prevVec = self.currVec = self._pointToVector(0, 0)
    self.update()
    
  def translate(self, dx, dy):
     
    if dx or dy:
      if self.moveRef:
        x, y, z = self.refPos
        self.refPos = [x+dx, y+dy, z]
        
      if self.moveMain:
        x, y, z = self.view
        self.view = [x+dx, y+dy, z]
        
      self.update()
  
  def _rotate(self, dx, dy):
  
    w, h = self.glGeometry[:2]
    w /= 2
    h /= 2
    self.prevVec = self._pointToVector(w, h)
    self.currVec = self._pointToVector(w+dx, h+dy)
    self.update()
   
  def rotateRight(self):
  
    self._rotate(self.translateStep, 0)
  
  def rotateLeft(self):

    self._rotate(-self.translateStep, 0)
    
  def rotateUp(self):

    self._rotate(0, -self.translateStep)
  
  def rotateDown(self):
  
    self._rotate(0, self.translateStep)

  def rotateView(self, angle, axis):
    """
    rotate the view by a given angle around a given axis
    hint: [0,1,0] is the vertical axis pointing upwards
          [1,0,0] is the horizontal axis pointing to the right
          [1,0,0] is [0,1,0] x [1,0,0] i.e. the axis pointing out of the screen
    """

    glPushMatrix()
    glLoadIdentity()
    glRotatef(angle, axis[0], axis[1], axis[2])
    glMultMatrixf(self.rotation)
    self.rotation = glGetFloatv(GL.GL_MODELVIEW_MATRIX)
    glPopMatrix()
    glMultMatrixf(self.rotation)
    
    self.paintEvent(None)
        
  def moveLeft(self):

    self.translate(self.translateStep,0)
    
  def moveRight(self):
    
    self.translate(-self.translateStep,0)

  def moveUp(self):

    self.translate(0,self.translateStep)
     
  def moveDown(self):
  
    self.translate(0,-self.translateStep)
      
  def construct(self, nuc=None):
    
    for i in self.glLists:
      GL.glDeleteLists(i, 1)
    
    self.lighting = set()
    
    glList = 1

    GL.glNewList(glList, GL.GL_COMPILE)
    quadric = GLU.gluNewQuadric()
    
    GL.glBegin(GL.GL_LINES)
    glColor4f(1.0, 0.0, 0.0, 1.0)
    glVertex3d(100,0,0)
    glVertex3d(0,0,0)
    
    glColor4f(0.0, 1.0, 0.0, 1.0)
    glVertex3d(0,100,0)
    glVertex3d(0,0,0)

    glColor4f(0.0, 0.0, 1.0, 1.0)
    glVertex3d(0,0,100)
    glVertex3d(0,0,0)
    GL.glEnd()

    glColor4f(1.0, 1.0, 1.0, 0.5)
    gluSphere(quadric, 14.0, 32, 32)
    
    GL.glEndList()
     
    self.glLists = [glList,]
    
    return glList

  def initGlParams(self):
    
    white = (1.0, 1.0, 1.0, 1.0)
    black = (0.0, 0.0, 0.0, 1.0)
    
    GL.glEnable(GL.GL_DEPTH_TEST)
    GL.glEnable(GL.GL_BLEND)
    GL.glEnable(GL.GL_LINE_SMOOTH)
    GL.glEnable(GL.GL_CLIP_PLANE0)
    GL.glEnable(GL.GL_CLIP_PLANE1)
    #GL.glEnable(GL.GL_CULL_FACE)
    GL.glEnable(GL.GL_MULTISAMPLE)
    GL.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA)
    GL.glDepthFunc(GL.GL_LEQUAL)
    GL.glEnable(GL.GL_SAMPLE_ALPHA_TO_COVERAGE)

    GL.glEnable(GL.GL_LIGHTING)
    GL.glEnable(GL.GL_LIGHT1)
    GL.glEnable(GL.GL_COLOR_MATERIAL)
    GL.glColorMaterial(GL.GL_FRONT, GL.GL_AMBIENT_AND_DIFFUSE)
 
    GL.glMaterialfv(GL.GL_FRONT, GL.GL_SPECULAR, white)
    GL.glMaterialfv(GL.GL_FRONT, GL.GL_EMISSION, black)
    GL.glMaterialf(GL.GL_FRONT, GL.GL_SHININESS, 100.0)
    GL.glLightfv(GL.GL_LIGHT1, GL.GL_AMBIENT,  black)
    GL.glLightfv(GL.GL_LIGHT1, GL.GL_DIFFUSE,  white)
    GL.glLightfv(GL.GL_LIGHT1, GL.GL_SPECULAR, white)
    GL.glLightModeli(GL.GL_LIGHT_MODEL_LOCAL_VIEWER, GL.GL_FALSE)

    #GL.glEnable(GL.GL_FOG)
    #GL.glFogi(GL.GL_FOG_MODE, GL.GL_EXP)
    #GL.glFogfv(GL.GL_FOG_COLOR, (0.0, 0.0, 0.0, 1.0))
    #GL.glFogf(GL.GL_FOG_DENSITY, 0.005);
    #GL.glFogf(GL.GL_FOG_START, 128.0)
    #GL.glFogf(GL.GL_FOG_END, 2048.0)
    #GL.glFogi(GL.GL_FOG_COORD_SRC, GL.GL_FRAGMENT_DEPTH)


    # For text glyphs
    """
    vShader = GL.shaders.compileShader(VERTEX_SHADER_GLSL, GL.GL_VERTEX_SHADER)  
    fShader = GL.shaders.compileShader(FRAGMENT_SHADER_GLSL, GL.GL_FRAGMENT_SHADER)
   
    shaderProg = GL.shaders.compileProgram(vShader, fShader)
    GL.shaders.glUseProgram(shaderProg)
    
    progTex = GL.glGetUniformLocation(shaderProg, 'tex')
    progCoord = GL.glGetAttribLocation(shaderProg, 'coord')

    GL.glActiveTexture(GL.GL_TEXTURE0)
    
    tex = GL.glGenTextures(1)
    GL.glBindTexture(GL.GL_TEXTURE_2D, tex)
    GL.glUniform1i(progTex, 0)
    
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
    
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
    
    GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)

    vbObj = GL.glGenBuffers(1)
    GL.glEnableVertexAttribArray(progCoord)
    GL.glBindBuffer(GL.GL_ARRAY_BUFFER, vbObj)
    GL.glVertexAttribPointer(progCoord, 4, GL.GL_FLOAT, GL.GL_FALSE, 0, 0)
    """
    
    # Temporary for chromosome labels
    GLUT.glutInit()
     
  def initializeGL(self):
  
    print(("OpenGL version:\t\t%s" % GL.glGetString(GL.GL_VERSION) ))
     
    self.construct(self.mainApp.nuc)
    self.rotation = eye(4).astype(float32)
    self.initGlParams()


  def _pointToVector(self, x, y):

    w, h = self.glGeometry[:2]
    w = float(w)
    h = float(h)
    
    vx = (2.0 * x - w) / w
    vy = (h - 2.0 * y) / h
    vz = cos(1.570796326794896 * sqrt(vx * vx + vy * vy))
    a = 1.0 / sqrt(vx * vx + vy * vy + vz * vz)
    
    return a * array([vx, vy, vz])
  
    
  def paintEvent(self, event):
    
    self.makeCurrent()
    self.initGlParams()    
    self.updateGlLayer()
    
    painter = QtGui.QPainter()
    painter.begin(self)
    
    
    #painter.beginNativePainting()
    
    #painter.endNativePainting()
   
    self.updateQtLayer(painter)
    
    painter.end()
  
  
  def updateQtLayer(self, painter):
  
    pen = QtGui.QPen(QtGui.QColor(255, 255, 0, 255))
    pen.setWidth(5)
    pen.setStyle(QtCore.Qt.SolidLine)
   
    painter.setPen(pen)
    painter.setBrush(QtGui.QColor(0, 128, 128, 128))
    
    painter.drawRect(QtCore.QRect(50, 50, 150, 150))
    w, h = self.glGeometry[:2]
    painter.drawLine(10,10,100,100)
    painter.drawText(QtCore.QPoint(50, 50), "hello")
  

  def updateGlLayer(self, doPick=False):
    
    # Perspective and view
    
    w, h, far, fov = self.glGeometry
    x, y, z = self.view
    z *= -32
    w = float(w)
    h = float(h)
    
    
    # Clear buffers
    
    if doPick:
      GL.glClearColor(0.0, 0.0, 0.0, 0.0)    
    
    else:
      GL.glClearColor(*self.bgColor)    
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    
    
    # Projection matrix
    
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(fov, w/h, 0.1, far)
    
    # Model matrix
    glMatrixMode(GL_MODELVIEW)

    # Update latest rotation vector, get rotation axis
    if not doPick:
      delta = self.currVec - self.prevVec
      angle = 90.0 * sqrt(dot(delta, delta))
      axis = cross(self.prevVec, self.currVec)
      self.prevVec = self.currVec
 
      # Create rotation matrices
      
      if self.moveRef:
        glLoadIdentity()
        glRotatef(angle, axis[0], axis[1], axis[2])
        glMultMatrixf(self.rotation2)
        self.rotation2 = glGetFloatv(GL.GL_MODELVIEW_MATRIX)
      
      if self.moveMain:
        glLoadIdentity()
        glRotatef(angle, axis[0], axis[1], axis[2])
        glMultMatrixf(self.rotation)
        self.rotation = glGetFloatv(GL.GL_MODELVIEW_MATRIX)
     
    # Move model in view
    glLoadIdentity()
    glTranslatef(-z*x/h, z*y/h, z)
    
    # Standard light source
    GL.glLightfv(GL.GL_LIGHT1, GL.GL_POSITION, (1.0, 1.0, 1.0, 0.0) )
 
    # Apply clipping
    limit = self.clipLimit
    if limit:
      glEnable(GL.GL_CLIP_PLANE0)
      glEnable(GL.GL_CLIP_PLANE1)
      GL.glClipPlane(GL.GL_CLIP_PLANE0, (0.0, 0.0, -1.0, self.clipPlaneFront*limit))
      GL.glClipPlane(GL.GL_CLIP_PLANE1, (0.0, 0.0, 1.0, self.clipPlaneBack*limit))
    
    # Accumulated rotation
    glMultMatrixf(self.rotation)    
     
    # Move model to origin
    ox, oy, oz = self.origin
    glTranslatef(-ox, -oy, -oz)    
    
    # Store matrices
    self.projectMatrix = glGetFloatv(GL.GL_PROJECTION_MATRIX)
    self.modelMatrix = glGetFloatv(GL.GL_MODELVIEW_MATRIX)
    
    # Call model vertex lists
    glDisable(GL.GL_LIGHTING)
    GL.glDisable(GL.GL_LINE_SMOOTH)
    GL.glDisable(GL.GL_MULTISAMPLE)
    if doPick:
       for i in self.glLists:
        if i in self.pickBuffer:
          glCallList(i)
    
    else:
      #GL.glDisable(GL.GL_BLEND)
      GL.glEnable(GL.GL_LINE_SMOOTH)
      GL.glEnable(GL.GL_MULTISAMPLE)
      
      # Call standard vertex lists
      #GL.glDepthMask(False)
      for i in self.glLists:
        if i in self.pickBuffer:
          continue
        
        #if i in self.reference:
        #  continue
        
        if i in self.lighting:
          glDisable(GL.GL_BLEND)
          glEnable(GL.GL_LIGHTING)
          glCallList(i)
          glDisable(GL.GL_LIGHTING)
          glEnable(GL.GL_BLEND)
 
        else:
          glCallList(i)
        
      #GL.glDepthMask(True)
      """
      if self.reference:
      
        GL.glDisable(GL.GL_LINE_SMOOTH)
        GL.glDisable(GL.GL_MULTISAMPLE)
       # Call separate, reference vertex lists
        glDisable(GL.GL_CLIP_PLANE0)
        glDisable(GL.GL_CLIP_PLANE1)
        x2, y2, z2 = self.refPos
       
        #vec = array([x2, y2, 0.0, 1.0])
        #vec = dot(-vec, linalg.inv(self.rotation2))
        #x3, y3, z3 = vec[:3]
 
        glLoadIdentity()
        glTranslatef(-z*(x+x2)/h, z*(y+y2)/h, z+z2)
        glMultMatrixf(self.rotation2)
        #glTranslatef(0, 0, z2)    
        
        glScalef(z2, z2, z2)
        glTranslatef(-ox, -oy, -oz)    
 
        for i in self.reference:
          if i in self.lighting:
            glEnable(GL.GL_LIGHTING)
            glCallList(i)
            glDisable(GL.GL_LIGHTING)
          else:
            glCallList(i)
      """
       
      GL.glEnable(GL.GL_MULTISAMPLE)
      GL.glEnable(GL.GL_LINE_SMOOTH)
      #GL.glEnable(GL.GL_BLEND)
    
    
  def resizeGL(self, width, height):
    
    self.glGeometry[:2] = width, height
    
    # Whole thing gets smaller, no clipping
    GL.glViewport(0, 0, width, height)
    
    
  def wheelEvent(self, event):
    
    mods = event.modifiers()
    haveCtrl = bool(mods & Qt.CTRL)
    haveShift = bool(mods & Qt.SHIFT)
    self.moveRef = haveShift or haveCtrl
    self.moveMain = not haveShift
    
    delta = event.delta()
    self.prevPos = QPoint(event.pos())
    self.mouseZoom(delta)
 
 


 
  def keyPressEvent(self, event):
  
    key = event.key()
    
    d_ref = 2.0
    a_ref = 3.14159265/90.0
    z_ref = 1.1
    
    if key == Qt.Key_PageUp:
      self.mouseZoom(-100)
      
    elif key == Qt.Key_PageDown:
      self.mouseZoom(100)

    elif key == Qt.Key_Up:
      self.rotateUp()
 
    elif key == Qt.Key_Down:
      self.rotateDown()
 
    elif key == Qt.Key_Left:
      self.rotateLeft()
 
    elif key == Qt.Key_Right:
      self.rotateRight()
    
    elif key == Qt.Key_N:
      self.image_scale *= z_ref
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)

    elif key == Qt.Key_M:
      self.image_scale /= z_ref
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)
    
    elif key == Qt.Key_E:
      self.image_translate[0] -= d_ref
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)

    elif key == Qt.Key_D:
      self.image_translate[1] -= d_ref
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)

    elif key == Qt.Key_C:
      self.image_translate[2] -= d_ref
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)
   
    elif key == Qt.Key_R:
      self.image_translate[0] += d_ref
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)

    elif key == Qt.Key_F:
      self.image_translate[1] += d_ref
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)

    elif key == Qt.Key_V:
      self.image_translate[2] += d_ref
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)
    
    elif key == Qt.Key_Q:
      r = getRotationMatrix([1.0, 0.0, 0.0], -a_ref)
      self.image_rotate = dot(self.image_rotate, r)
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)
      
    elif key == Qt.Key_A:
      r = getRotationMatrix([0.0, 1.0, 0.0], -a_ref)
      self.image_rotate = dot(self.image_rotate, r)
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)
      
    elif key == Qt.Key_Z:
      r = getRotationMatrix([0.0, 0.0, 1.0], -a_ref)
      self.image_rotate = dot(self.image_rotate, r)
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)
      
    elif key == Qt.Key_W:
      r = getRotationMatrix([1.0, 0.0, 0.0], a_ref)
      self.image_rotate = dot(self.image_rotate, r)
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)
      
    elif key == Qt.Key_S:
      r = getRotationMatrix([0.0, 1.0, 0.0], a_ref)
      self.image_rotate = dot(self.image_rotate, r)
      self.construct(self.mainApp.nuc)
      self.paintEvent(None)

    elif key == Qt.Key_X:
      r = getRotationMatrix([0.0, 0.0, 1.0], a_ref)
      self.image_rotate = dot(self.image_rotate, r)
      self.construct(self.mainApp.nuc)
      self.update()
      #self.paintEvent(None)


  def destroy(self, *args):
  
    QtOpenGL.QGLWidget.destroy(self, args)
  
  
  def leaveEvent(self, event):
    
    pass
    
    
  def hideEvent(self, event):

    pass


  def enterEvent(self, event):
    
    self.setFocus(Qt.OtherFocusReason)
  
  
  def _mouseRotate(self, event):
    
    mods = event.modifiers()
    haveCtrl = bool(mods & Qt.CTRL)
    haveShift = bool(mods & Qt.SHIFT)
  
    pos = event.pos()
    self.currVec = self._pointToVector(pos.x(), pos.y())
    self.moveRef = not haveCtrl # haveShift
    self.moveMain = not haveShift
    
    self.update()
  
  
  def _mouseTranslate(self, event):
  
    dx = event.x() - self.prevPos.x()
    dy = event.y() - self.prevPos.y()
    
    mods = event.modifiers()
    haveCtrl = bool(mods & Qt.CTRL)
    haveShift = bool(mods & Qt.SHIFT)
    
    self.moveRef = not haveCtrl
    self.moveMain = not haveShift
    self.translate(dx, dy)
    
    self.prevVec = self.currVec
  
  
  def mouseReleaseEvent(self, event):
    
    button = event.button()
    
    if button == LeftButton:
      self.makeCurrent()
   
      GL.glDisable(GL.GL_BLEND)
      GL.glDisable(GL.GL_LINE_SMOOTH)
      GL.glDisable(GL.GL_MULTISAMPLE)
      GL.glEnable(GL.GL_DEPTH_TEST)
      
      self.updateGlLayer(doPick=True)  # Render only pick buffer
      self.paintGL()                   # Lowest level paint - doesn't actually update screen 
     
      pos = event.pos()
      x = float(pos.x())
      y = float(pos.y())
   
      viewport = GL.glGetIntegerv(GL.GL_VIEWPORT)
      pickColor = GL.glReadPixels(x, viewport[3] - y, 1, 1, GL.GL_RGBA, GL.GL_UNSIGNED_BYTE)
      
      r, g, b, a = ord(pickColor[0]), ord(pickColor[1]), ord(pickColor[2]), ord(pickColor[3])
      index = a + b*255 + g*65025 + r*16581375
      
      GL.glEnable(GL.GL_MULTISAMPLE)
      GL.glEnable(GL.GL_LINE_SMOOTH)
      GL.glEnable(GL.GL_BLEND)
    
    else:
      index = 0   
      
    self.prevVec = self.currVec
    self.moveRef = False
    self.moveMain = True
    
    return index
    
    
  def mousePressEvent(self, event):
  
    pos = event.pos()
    self.currVec = self.prevVec = self._pointToVector(pos.x(), pos.y())
    self.prevPos = QPoint(event.pos())

    mods = event.modifiers()
    haveCtrl = bool(mods & Qt.CTRL)
    haveShift = bool(mods & Qt.SHIFT)
    self.moveRef = not haveCtrl
    self.moveMain = not haveShift


  def mouseMoveEvent(self, event):
    
    buttons = event.buttons()
    
    if buttons & LeftButton:
      self.mouseLeftMoveFunc(event)
      
    elif buttons & MiddleButton:
      self.mouseMiddleMoveFunc(event)
      
    elif buttons & RightButton:
      self.mouseRightMoveFunc(event)

    self.prevPos = QPoint(event.pos())
