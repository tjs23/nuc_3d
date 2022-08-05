import sys

def warn(msg, prefix='WARNING'):

  print('%8s : %s' % (prefix, msg))


def testImports(gui=False):
  
  critical = False
  cWarn    = False
  guiWarn  = False
    
  try:
    from h5py import File
  except ImportError:
    critical = True
    warn('Critical module "h5py" is not installed or accessible')
  
  try:
    from numpy import array
  except ImportError:
    critical = True
    warn('Critical module "numpy" is not installed or accessible')
  
  try:
    from cUtil.apiUtil import calcCoordDensity
  except ImportError:
    cWarn = True
    warn('Utility C/Cython code not compiled')
  
  try:
    from cUtil.samread import readPairedSam
  except ImportError as err:
    warn('Samtools code or its interface not compiled')
    warn('Importing BAM/SAM format files will not work', 'INFO')
    raise(err)
    
  try:
    from solve.SimAnneal import runDynamics
  except ImportError:
    cWarn = True
    warn('Structure calculation C/Cython code not compiled')
    
  try:
    from PySide import QtCore, QtGui
  except ImportError:
    guiWarn = True
    warn('GUI module "PySide" is not installed or accessible\n')
  
  try:
    from OpenGL import GL, GLU
  except ImportError:
    guiWarn = True
    warn('GUI module "OpenGL" is not installed or accessible\n')  
    
  if critical:
    warn('Nuc3d cannot proceed because critical Python modules are not available', 'ABORT')
  
  if cWarn:
    warn('Nuc3d cannot proceed because Cython code is not compiled', 'ABORT')
    warn('To compile the Cython code run the script "./compile.sh"', 'INFO')

  if gui and guiWarn:
    warn('Nuc3d cannot proceed in GUI mode because critical Python modules are not available', 'ABORT')
  
  if critical or cWarn or (gui and guiWarn):
    sys.exit(0)
    

from PySide import QtNetwork
from parallel import Controller
import loadLibs

testImports(gui=True)

from NucGui import main

main()
