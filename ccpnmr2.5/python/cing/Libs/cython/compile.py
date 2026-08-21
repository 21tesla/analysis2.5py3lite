# usage: python compile.py build_ext --inplace
# TODO automate this from CING setup script.
# Modernized for Python 3.13 (py2-era `distutils` / `Cython.Distutils` were removed);
# requires Cython >= 3 and setuptools.
from Cython.Build import cythonize
from setuptools import Extension, setup

setup(
  name = 'Superpose',
  ext_modules=cythonize([Extension("superpose", ["superpose.pyx"])], language_level=3)
)
