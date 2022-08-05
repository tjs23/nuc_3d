import os, sys
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

from numpy import get_include as numpy_get_include
numpy_includes = numpy_get_include()

libbam = os.path.join(os.path.dirname(__file__), "samtools", "libbam.a")

ext_modules = [Extension("drawing",    ["drawing.pyx"],    libraries=['m'],
                         include_dirs=['.','cUtil', numpy_includes]),
               Extension("dataLayer",  ["dataLayer.pyx"],  libraries=['m'],
                         include_dirs=['.','cUtil', numpy_includes]),
               Extension("kdTree",  ["kdTree.pyx"],  libraries=['m'],
                         include_dirs=['.','cUtil', numpy_includes]),
               Extension("contour",    ["contour.pyx"],    libraries=['m'],
                         include_dirs=[numpy_includes]),
               Extension("apiUtil",    ["apiUtil.pyx"],    libraries=['m'],
                         include_dirs=[numpy_includes]),
               Extension("getCoordsUtil",    ["getCoordsUtil.pyx"],  libraries=['m'],
                         include_dirs=[numpy_includes]),
               Extension("geometry",   ["geometry.pyx"],   libraries=['m']),
               Extension("overlapHelper",   ["overlapHelper.pyx"],   libraries=['m'],
                         include_dirs=[numpy_includes]),
               Extension("samread",    ["samread.pyx"], libraries=['z'],
                         extra_objects=[libbam], include_dirs=[numpy_includes])
]

if sys.platform == 'darwin':
  ext_modules.append(Extension("OpenGlUtil", ["OpenGlUtil.pyx"], libraries=['m', 'GL', 'GLU'],
                               include_dirs=[numpy_includes, "/usr/X11/include"],
                               library_dirs=["/usr/X11/lib/"]))
else:
  ext_modules.append(Extension("OpenGlUtil", ["OpenGlUtil.pyx"], libraries=['m', 'GL', 'GLU'],
                               include_dirs=[numpy_includes]))

setup(name = 'Nuc3D cUtil',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules)
