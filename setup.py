#!/usr/bin/env python
import sys
from setuptools import setup, find_packages
from Cython.Build import cythonize
from Cython.Build.Dependencies import default_create_extension
import numpy


def create_extension(template, kwds: dict):
    define_macros = kwds.get("define_macros", [])
    
    #Use the new numpy API and remove all the compilation warnings
    define_macros.append(("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"))

    if sys.platform in ("darwin", "win32"):
        define_macros.append(("CYTHON_INLINE", ""))

    kwds["define_macros"] = define_macros
    return default_create_extension(template, kwds)


setup(
    ext_modules=cythonize("pyoptools/**/*.pyx", create_extension=create_extension,
                          language_level="3str"),
    include_dirs=[numpy.get_include()],   
    use_scm_version=True,
)
