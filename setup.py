"""Build script for the ShapeFile Python extension (Python 3.13 C API spike)."""
import os
from setuptools import setup, Extension

CC_DIR = os.path.join("ccpnmr2.5", "c", "memops", "global")

shape_file_ext = Extension(
    "ShapeFile",
    sources=[
        os.path.join(CC_DIR, "py_shape_file.c"),
        os.path.join(CC_DIR, "shape_file.c"),
        os.path.join(CC_DIR, "python_util.c"),
        os.path.join(CC_DIR, "utility.c"),
    ],
    include_dirs=[CC_DIR],
    extra_compile_args=["-Wall", "-Wno-unused-function"],
)

setup(
    name="ShapeFile",
    version="2.5.2",
    ext_modules=[shape_file_ext],
)
