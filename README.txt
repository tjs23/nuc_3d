
REQUIRED PYTHON PACKAGES
 
 NumPy    - Numeric Python; generally bundled with SciPy
 H5py     - Hierarchical Data Format v5 bindings for Python
 Cython   - Python to C translator
 PyOpenGl - OpenGL bindings for Python
 PySide   - Qt GUI toolkit Python bindings - Will Change to PyQt5
 MPI4Py   - Message passing interface bindings for Python

 Example installation on Debian based systems:

   sudo apt-get install cython python-opengl python-pyside python-scipy python-h5py python-mpi4py

 Note: PySam is no longer a dependency as SAMtools is directly bundled
 with a bespoke Python interface.


COMPILING CYTHON/C CODE
 
 The following script runs Cython, to generate the C code,
 and then compiles the C code into .so Python modules:
 
 ./compile.sh
 


DIRECTORY CONTENTS

 nuc3d        - Shell script to start the main GUI or command line system.
                Sets environment variables etc.

 hdfView      - Shell script to start HdfViewer.py a graphical program to
                view the contents of HDF5 files (inlcuding .nuc files)
                
 nuc3d.py     - Python script to support the nuc3d executable.
                Handles command line arguments.

 NucApi.py    - Contains the Nucleus class, API and related non-GUI functions
 
 NucGui.py    - The startup script for the graphical interface

 analyses/    - Python code for high-level analysis of .nuc files
 
 cUtil/       - Cython and C code, except for structure calculation
 
 gui/         - Python graphical interface code
 gui/qtgui/   - Python convenience wrappings around PySide/Qt classes
 
 data/        - Data files for testing, research and reference
 
 solve/       - Cython and C code for structure calculation
 
 test/        - Code validation and testing scripts
 
 util/        - Python code with no dependency on the Nucleus class 
