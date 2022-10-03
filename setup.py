#!/usr/bin/env python
import sys

from setuptools import setup, find_packages
from Cython.Build import cythonize
from Cython.Build.Dependencies import default_create_extension

import numpy


def create_extension(template, kwds: dict):
    define_macros = kwds.get("define_macros", [])

    if sys.platform in ("darwin", "win32"):
        define_macros.append(("CYTHON_INLINE", ""))
    
    kwds["define_macros"] = define_macros
    return default_create_extension(template, kwds)


setup(
    name="pyoptools",
    version="0.1.1",
    packages=find_packages(exclude=["tests"]),
    scripts=["ipyoptools"],
    package_data={
        "pyoptools.raytrace.mat_lib": [
            "data/glass/*",
            "data/glass/*/*",
            "data/glass/*/*/*",
            "data/main*",
            "data/main/*/*",
            "data/aliases.cfg",
        ],
        "pyoptools.raytrace.library": ["Edmund/*.cmp", "Thorlabs/*.cmp"],
    },
    author="Ricardo Amezquita Orozco",
    author_email="ramezquitao@cihologramas.com",
    description="Optical ray tracing simulation system",
    license="GPLv3",
    url="https://github.com/cihologramas/pyoptools/",
    download_url="https://github.com/cihologramas/pyoptools/archive/v0.1.1.zip",
    ext_modules=cythonize("pyoptools/**/*.pyx", language_level="2", create_extension=create_extension),
    include_dirs=[numpy.get_include()],
    install_requires=['numpy', 'scipy', 'imageio', 'PyYAML', 'matplotlib'],
)
