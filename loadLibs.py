import os, sys
from ctypes import cdll

# Load C libraries directly via Python

if sys.platform == 'darwin':
  libbam = 'libbam.dylib'

else:
  libbam = 'libbam.so'


nucDir = os.path.dirname(os.path.abspath(__file__))
libbam = os.path.join(nucDir, 'cUtil', 'samtools', libbam)
cdll.LoadLibrary(libbam)
