#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import sys
from sage.env import sage_include_directories

extensions = Extension("betti.bidegmap", ["betti/bidegmap.pyx"],
    extra_compile_args=["-march=native", "-funroll-loops"])

extensions=cythonize(extensions, include_path=sys.path)

# Run Distutils
setup(
    name="betti",
    version='0.1',
    packages=["betti"],
    ext_modules=extensions,
    include_dirs=sage_include_directories()
)
