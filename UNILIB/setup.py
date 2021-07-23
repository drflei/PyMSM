# setup.py file
from distutils.core import setup, Extension
from Cython.Build   import cythonize
#
#
# Note in the invocation below, and a quirk with Cythonize, the module name
# must equal to the prefix name on the Cython (*.pyx) file.
#
# For The compiled UNILIB library is assumed to be in "./lib/libunilib.so".
# PyUNILIB.pyx has now been added to ./src/, and a header file for the wrapper
# interface is in PyUNILIB_Interface.h.  The UNILIB shared object library
# should be built by the makefile before this Python script is invoked.
#
# Cython also needs the fortran library "libgfortran.so" in my installation
# located in /usr/lib64/
#
#
# This script is invoked as:
#
#   python3 setup.py build_ext --inplace
#
# and the resulting PyUNILIB.*.so file will be in ./lib64
#
setup(
  ext_modules = cythonize(
    Extension(
      name         = "lib.PyUNILIB",
      sources      = ["./src/PyUNILIB.pyx"],
      libraries    = ["unilib", "gfortran"],
      library_dirs = ["./lib/", "/usr/lib64/"],
      language="c",
    )
  )
)
