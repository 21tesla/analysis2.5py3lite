"""Build script for CCPNMR C extensions (Python 3.13)."""
import os
from setuptools import setup, Extension

CC = "ccpnmr2.5/c/memops/global"

def ext(name, sources, extra_sources=None):
    srcs = [os.path.join(CC, s) for s in sources]
    if extra_sources:
        srcs += [os.path.join(CC, s) for s in extra_sources]
    return Extension(
        name,
        sources=srcs,
        include_dirs=[CC],
        extra_compile_args=["-Wall", "-Wno-unused-function", "-Wno-unused-variable"],
    )

ext_modules = [
    ext("ShapeFile",
        ["py_shape_file.c", "shape_file.c", "python_util.c", "utility.c"]),

    ext("MemCache",
        ["py_mem_cache.c", "mem_cache.c",
         "hash_list.c", "hash_table.c", "int_array.c",
         "list.c", "mutex.c",
         "python_util.c", "utility.c"]),

    ext("BlockFile",
        ["py_block_file.c", "block_file.c",
         "py_mem_cache.c", "py_shape_file.c",
         "hash_list.c", "hash_table.c", "int_array.c",
         "list.c", "mutex.c", "mem_cache.c", "shape_file.c",
         "python_util.c", "utility.c"]),

    ext("FitMethod",
        ["py_fit.c", "fit.c", "fit1d.c", "nonlinear_model.c",
         "cpmg.c", "line_fit.c", "random.c", "gauss_jordan.c",
         "gamma.c",
         "python_util.c", "utility.c"]),

    ext("StoreFile",
        ["py_store_file.c", "store_file.c", "python_util.c", "utility.c"]),

    ext("StoreHandler",
        ["py_store_handler.c", "store_handler.c", "python_util.c", "utility.c"]),

    ext("PdfHandler",
        ["py_pdf_handler.c", "pdf_handler.c", "clipping.c",
         "python_util.c", "utility.c"]),

    ext("PsHandler",
        ["py_ps_handler.c", "ps_handler.c", "clipping.c",
         "python_util.c", "utility.c"]),
]

setup(
    name="ccpnmr-ext",
    version="2.5.2",
    ext_modules=ext_modules,
    zip_safe=False,
)
