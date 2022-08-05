#!/usr/bin/python3.3
from __future__ import division

import itertools
import os.path
from math import radians
from ctypes import c_void_p as void_p

from PySide import QtGui, QtOpenGL, QtCore
from OpenGL import GL
from OpenGL.GL import shaders
import numpy

from gui.qtgui.FileSelect import FileType, selectFile

from util import XForm
from util.Arcball import ArcBall

try:
  import tifffile
except ImportError:
  import util.tifffile as tifffile

class GLTexture:
  def __init__(self, handle, sampler):
    self.handle = handle
    self.sampler = sampler

  @property
  def image_unit(self):
    return self.sampler + GL.GL_TEXTURE0

class ImagePanel(QtGui.QWidget):
  def __init__(self, parent=None):
    super(self.__class__, self).__init__(parent)

    # Widget parents are set when added to a layout
    settings_layout = QtGui.QHBoxLayout()
    open_button = QtGui.QPushButton("Open")
    open_button.clicked.connect(self.openFile)
    settings_layout.addWidget(open_button)
    gain_slider = QtGui.QSlider(QtCore.Qt.Horizontal)
    settings_layout.addWidget(gain_slider)
    bg_selector = QtGui.QComboBox()
    bg_selector.addItems(["Black", "White"])
    settings_layout.addWidget(bg_selector)
    blend_selector = QtGui.QComboBox()
    blend_selector.addItems(["Linear", "Multiplicative", "Maximum"])
    settings_layout.addWidget(blend_selector)

    layout = QtGui.QVBoxLayout(self)
    layout.addLayout(settings_layout)
    self.display = ImageDisplay()
    layout.addWidget(self.display)

    gain_slider.valueChanged.connect(self.display.setGain)
    blend_selector.activated[str].connect(self.display.setBlend)
    bg_selector.activated[str].connect(self.display.setBackground)
    self.setLayout(layout)

  @QtCore.Slot()
  def openFile(self):
    file_name = selectFile(self, fileTypes=ImageDisplay.valid_file_types)
    if file_name is not None:
      self.display.openFile(file_name)

class ImageDisplay(QtOpenGL.QGLWidget):
  valid_file_types = [FileType("TIFF Image", ["*.tif", "*.tiff"]),
                      FileType("Raw Metadata", ["*.dat"]),]
   
  def __init__(self, parent=None):
    super(self.__class__, self).__init__(parent)
    
    self.cam_distance = 3
    self.fov = radians(50)
    self.background = 'black'
    self.blend = 'linear'
    self.arcball = ArcBall()
    
    self.model_world_xform = numpy.identity(4)
    size = self.geometry().size()
    self.camera_clip_xform = XForm.perspective(self.fov, aspect=size.width()/size.height())

    self.render_stages = ['volume', 'post']

  def sizeHint(self):
    return QtCore.QSize(400, 300) 
  
  def initializeGL(self):
    self.gl_programs = {k: GL.glCreateProgram() for k in self.render_stages}

    shader_dir = os.path.join(os.path.dirname(__file__), 'ImagePanelShaders')
      
    self.attribs = {k: {'position': 0} for k in self.render_stages}
    self.programs = {k: GL.glCreateProgram() for k in self.render_stages}
    for render_stage in self.render_stages:
      shaders = {'vert': GL.glCreateShader(GL.GL_VERTEX_SHADER),
                 'frag': GL.glCreateShader(GL.GL_FRAGMENT_SHADER)}
      for shader_type, shader in shaders.items():
        shader_filename = os.path.join(shader_dir, '.'.join((render_stage, shader_type)))
        with open(shader_filename) as src:
          GL.glShaderSource(shader, src.read())
        GL.glCompileShader(shader)
        GL.glAttachShader(self.programs[render_stage], shader)
      for attrib, location in self.attribs[render_stage].items():
        GL.glBindAttribLocation(self.programs[render_stage], location, attrib)
      GL.glLinkProgram(self.programs[render_stage])

    uniform_names = {'volume': ['camera_clip_xform', 'world_camera_xform', 'model_world_xform',
                                'slices', 'axis', 'gain', 'volume_sampler'],
                     'post': ['fb_sampler', 'invert']}
    self.uniforms = {render_stage: {name: GL.glGetUniformLocation(self.programs[render_stage], name) 
                                    for name in uniform_names[render_stage]}
                     for render_stage in self.render_stages}
    
    # Might be able to use same sampler since used in different programs?
    self.volume_tex = GLTexture(handle=GL.glGenTextures(1), sampler=1)
    self.fb_tex = GLTexture(handle=GL.glGenTextures(1), sampler=2)

    # Set-up volume texture
    GL.glBindTexture(GL.GL_TEXTURE_3D, self.volume_tex.handle)
    GL.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1)
    GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
    GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
    GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_BORDER)
    GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_BORDER)
    GL.glTexParameteri(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_WRAP_R, GL.GL_CLAMP_TO_BORDER)
    GL.glTexParameterfv(GL.GL_TEXTURE_3D, GL.GL_TEXTURE_BORDER_COLOR, numpy.zeros(4))
    GL.glBindTexture(GL.GL_TEXTURE_3D, 0)
    
    # Set-up framebuffer
    self.fb = GL.glGenFramebuffers(1)
    GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fb)
    GL.glBindTexture(GL.GL_TEXTURE_2D, self.fb_tex.handle)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_BORDER)
    GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_BORDER)
    GL.glTexParameterfv(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_BORDER_COLOR, numpy.zeros(4))
    size = self.geometry().size()
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RED, size.width(), size.height(), 0, GL.GL_RED, GL.GL_FLOAT, None)
    GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0, GL.GL_TEXTURE_2D, self.fb_tex.handle, 0)
    GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)
    GL.glBindTexture(GL.GL_TEXTURE_2D, 0)
    
    # Both post and volume sampling works from unit square as geometry
    self.vao = {k: v for k, v in zip(self.render_stages, GL.glGenVertexArrays(2))}
    self.vbo = {k: v for k, v in zip(self.render_stages, GL.glGenBuffers(2))}
    unit_square = numpy.array(list(itertools.product((-1, 1), (-1, 1))), dtype='f')
    for render_stage in self.render_stages:
      GL.glBindVertexArray(self.vao[render_stage])
      GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self.vbo[render_stage])
      GL.glBufferData(GL.GL_ARRAY_BUFFER, unit_square.nbytes, unit_square, GL.GL_STATIC_DRAW)
      GL.glEnableVertexAttribArray(self.attribs[render_stage]['position'])
      GL.glVertexAttribPointer(self.attribs[render_stage]['position'], 2, GL.GL_FLOAT, False,
                               unit_square.itemsize * unit_square.shape[1], void_p(0))

      GL.glBindBuffer(GL.GL_ARRAY_BUFFER, 0)
      GL.glBindVertexArray(0)
      GL.glDisableVertexAttribArray(self.attribs[render_stage]['position'])
    
    GL.glUseProgram(self.programs['volume'])
    world_camera_xform = XForm.lookAt((0, 0, self.cam_distance)).dot(self.arcball.totalRotation())
    GL.glUniformMatrix4fv(self.uniforms['volume']['world_camera_xform'], 1, True, world_camera_xform.astype('f'))
    volume_xform = self.camera_clip_xform.dot(world_camera_xform).dot(self.model_world_xform)
    GL.glUniform1i(self.uniforms['volume']['axis'], XForm.viewAxis(volume_xform))
    GL.glUniform1i(self.uniforms['volume']['volume_sampler'], self.volume_tex.sampler)
    self.num_volume_slices = 100
    GL.glUniform1i(self.uniforms['volume']['slices'], self.num_volume_slices)
    
    GL.glUseProgram(self.programs['post'])
    GL.glUniform1i(self.uniforms['post']['fb_sampler'], self.fb_tex.sampler)
    GL.glUseProgram(0)
    
    GL.glActiveTexture(self.volume_tex.image_unit)
    GL.glBindTexture(GL.GL_TEXTURE_3D, self.volume_tex.handle)
    GL.glActiveTexture(self.fb_tex.image_unit)
    GL.glBindTexture(GL.GL_TEXTURE_2D, self.fb_tex.handle)
    GL.glActiveTexture(GL.GL_TEXTURE0)
    GL.glBindTexture(GL.GL_TEXTURE_3D, 0)
    GL.glBindTexture(GL.GL_TEXTURE_2D, 0)
    
    GL.glEnable(GL.GL_BLEND)
    self.setBlend(self.blend)
    self.setBackground(self.background)
    self.update()

  def openFile(self, filename):
    _, ext = os.path.splitext(filename)
    ext = ext.lower()
    if ext in (".tif", ".tiff"):
      tiff = tifffile.TiffFile(filename)
      axis_scales = itertools.repeat(1.0)
      if len(tiff.series) > 1:
        raise ValueError("{} series found in '{}' (expected 1).".format(len(tiff.series), filename))
      image_axes = tiff.series[0].axes
      # See `tifffile.TiffPage` and `tifffile.AXES_LABELS` for definitions
      spatial_axis_labels = set('XYZPI')
      num_spatial_axes = len(set(image_axes) & spatial_axis_labels)
      # TODO: Account for 2D images
      if num_spatial_axes != 3:
        raise ValueError("{} spatial axes found in '{}' (expected 3).".format(num_spatial_axes, filename))

      data = tiff.asarray()
      # TODO: Account for multiple data channels in one image
      data = data[[slice(None) if a in spatial_axis_labels else slice(1)
                   for a in image_axes]]
      data = numpy.squeeze(data)
    elif ext == ".dat":
      with open(filename) as f:
        dat = {k.strip(): v.strip() for k, v in (line.split(':') for line in f)}
      raw_path = os.path.join(os.path.dirname(filename), dat['ObjectFileName'])
      try:
        axis_scales = (float(s) for s in dat['SliceThickness'].split())
      except KeyError:
        axis_scales = itertools.repeat(1.0)
      resolution = [int(i) for i in dat['Resolution'].split()]
      with open(raw_path) as f:
        data = numpy.fromfile(f, dtype=dat['Format'].lower())
      data = data.reshape(*reversed(resolution))
    else:
      raise ValueError("Unrecognized file type.")

    axis_scales = [s * r for s, r in zip(axis_scales, reversed(data.shape))]
    axis_scales = [i / max(axis_scales) for i in axis_scales]
    self.setData(data, axis_scales)

  def setData(self, data, axis_scales):
    max_3D_size = GL.glGetInteger(GL.GL_MAX_3D_TEXTURE_SIZE)
    if any((s > max_3D_size for s in data.shape)):
      raise ValueError("Data set ({}) too large (max. {:d} allowed in any dimension)."\
        .format('x'.join(str(s) for s in data.shape), max_3D_size))
    
    data_types = {
      numpy.dtype(numpy.uint8): GL.GL_UNSIGNED_BYTE,
      numpy.dtype(numpy.uint16): GL.GL_UNSIGNED_SHORT,
      numpy.dtype(numpy.uint32): GL.GL_UNSIGNED_INT,
      numpy.dtype(numpy.float32): GL.GL_FLOAT
    }
    
    self.model_world_xform = XForm.scale(*axis_scales)
    GL.glUseProgram(self.programs['volume'])
    GL.glUniformMatrix4fv(self.uniforms['volume']['model_world_xform'], 1, True, self.model_world_xform.astype('f'))
    GL.glUseProgram(0)
    
    GL.glBindTexture(GL.GL_TEXTURE_3D, self.volume_tex.handle)
    GL.glTexImage3D(GL.GL_TEXTURE_3D, 0, GL.GL_RED, data.shape[2], data.shape[1], data.shape[0], 0, GL.GL_RED, data_types[data.dtype], data)
    GL.glBindTexture(GL.GL_TEXTURE_3D, 0)
  
  def resizeGL(self, w, h):
    self.camera_clip_xform = XForm.perspective(self.fov, aspect=w/h)
    
    GL.glUseProgram(self.programs['volume'])
    GL.glUniformMatrix4fv(self.uniforms['volume']['camera_clip_xform'], 1, True, self.camera_clip_xform.astype('f'))
    GL.glUseProgram(0)
    
    GL.glActiveTexture(self.fb_tex.image_unit)
    GL.glBindTexture(GL.GL_TEXTURE_2D, self.fb_tex.handle)
    GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RED, w, h, 0, GL.GL_RGB, GL.GL_FLOAT, None)
    GL.glActiveTexture(GL.GL_TEXTURE0)
    GL.glBindTexture(GL.GL_TEXTURE_2D, 0)
    
    GL.glViewport(0, 0, w, h)
  
  def paintGL(self):
    GL.glDisable(GL.GL_CULL_FACE)
    GL.glDepthMask(GL.GL_FALSE)
    
    GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.fb)
    GL.glEnable(GL.GL_BLEND)
    GL.glClear(GL.GL_COLOR_BUFFER_BIT)
    
    GL.glBindVertexArray(self.vao['volume'])
    GL.glUseProgram(self.programs['volume'])
    
    GL.glActiveTexture(self.volume_tex.image_unit) # FIXME: Don't think this is necessary
    GL.glDrawArraysInstanced(GL.GL_TRIANGLE_STRIP, 0, 4, self.num_volume_slices)
    
    GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, 0)
    GL.glDisable(GL.GL_BLEND)
    GL.glBindVertexArray(self.vao['post'])
    GL.glUseProgram(self.programs['post'])
    
    GL.glActiveTexture(self.fb_tex.image_unit) #FIXME: Don't think this is necessary
    GL.glDrawArrays(GL.GL_TRIANGLE_STRIP, 0, 4)
    
    GL.glBindVertexArray(0)
    GL.glUseProgram(0)
  
  def mousePressEvent(self, event):
    pos = numpy.array((event.posF().x() / self.geometry().size().width(),
                       event.posF().y() / self.geometry().size().height()))
    pos = pos * 2 - 1
    pos[1] *= -1
    self.arcball.startRotation(pos)
  
  def mouseMoveEvent(self, event):
    pos = numpy.array((event.posF().x() / self.geometry().size().width(),
                       event.posF().y() / self.geometry().size().height()))
    pos = pos * 2 - 1
    pos[1] *= -1
    self.arcball.updateRotation(pos)
    
    world_camera_xform = XForm.lookAt((0, 0, self.cam_distance)).dot(self.arcball.totalRotation())
    
    volume_xform = self.camera_clip_xform\
      .dot(world_camera_xform)\
      .dot(self.model_world_xform)
    volume_view_axis = XForm.viewAxis(volume_xform)
    
    GL.glUseProgram(self.programs['volume'])
    GL.glUniformMatrix4fv(self.uniforms['volume']['world_camera_xform'], 1, True,
      world_camera_xform.astype('f'))
    GL.glUniform1i(self.uniforms['volume']['axis'], volume_view_axis)
    GL.glUseProgram(0)
    
    self.update()

  def mouseReleaseEvent(self, event):
    self.arcball.finishRotation()
  
  @QtCore.Slot(str)
  def setBackground(self, bg):
    bg = bg.lower()
    self.background = bg

    GL.glUseProgram(self.programs['post'])
    if bg == 'white':
      GL.glUniform1i(self.uniforms['post']['invert'], 1)
    elif bg == 'black':
      GL.glUniform1i(self.uniforms['post']['invert'], 0)
    else:
      raise ValueError("Invalid background: {}".format(bg))
    GL.glUseProgram(0)
    
    self.update()
  
  @QtCore.Slot(str)
  def setBlend(self, blend):
    blend = blend.lower()
    self.blend = blend
    
    if blend == 'linear':
      GL.glBlendEquation(GL.GL_FUNC_ADD)
      GL.glBlendFunc(GL.GL_ONE, GL.GL_ONE)
    elif blend == 'multiplicative':
      GL.glBlendEquation(GL.GL_FUNC_ADD)
      GL.glBlendFunc(GL.GL_ONE_MINUS_DST_COLOR, GL.GL_ONE)
    elif blend == 'maximum':
      GL.glBlendEquation(GL.GL_MAX)
    else:
      raise ValueError("Invalid blend type: {}".format(blend))
    
    self.update()
 
  @QtCore.Slot(int)
  def setGain(self, value):
    gain = 10 ** (value / 100.0 * 6 - 3)
    # Applying in post program doesn't work (looks like it clips). Not sure why.
    GL.glUseProgram(self.programs['volume'])
    GL.glUniform1f(self.uniforms['volume']['gain'], gain)
    GL.glUseProgram(0)
    self.update()
  
  @QtCore.Slot(int, int, int)
  def setColor(self, r, g, b):
    pass
