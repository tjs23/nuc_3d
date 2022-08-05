from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

from numpy import get_include as numpy_get_include
numpy_includes = numpy_get_include()

ext_modules = [ Extension("SimAnneal", ["SimAnneal.pyx"], libraries=["m"],
                         include_dirs=[numpy_includes]) ]

setup(
  name = 'Nuc3D Simulated Annealing',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
