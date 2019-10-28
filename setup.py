#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from setuptools import setup
from setuptools.extension import Extension
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    print("You don't seem to have Cython installed. Please get a")
    print("copy from www.cython.org and install it")
    sys.exit(1)


# Look for paths containing arrayobject.h
# Necessary for non-unix systems
def contains_arrayobject_h(path):
    """
    Returns True if the python path string contains the arrayobject.h
    include file where it is supposed to be.
    """
    f = False
    try:
        _s = os.stat(os.path.join(path, 'numpy', 'core', 'include',
                                  'numpy', 'arrayobject.h'))
        f = True
    except OSError:
        pass
    return f


# scan the directory for extension files, converting
# them to extension names in dotted notation
def scandir(dir_, files=[]):
    for file in os.listdir(dir_):
        path = os.path.join(dir_, file)
        if os.path.isfile(path) and path.endswith(".pyx"):
            files.append(path.replace(os.path.sep, ".")[2:-4])
        elif os.path.isdir(path):
            scandir(path, files)
    return files


def findpackages(dir_, files=[]):
    for file in os.listdir(dir_):
        if file != "build":
            path = os.path.join(dir_, file)
            if os.path.isdir(path):  # and path.endswith(".py"):
                for file1 in os.listdir(path):
                    if file1 == "__init__.py":
                        files.append(path.replace(os.path.sep, ".")[2:])
                findpackages(path, files)
    return files


# generate an Extension object from its dotted name
def makeExtension(extName):
    extPath = extName.replace(".", os.path.sep)+".pyx"
    return Extension(
        extName,
        [extPath],
        # adding the '.' to include_dirs is CRUCIAL!!
        include_dirs=[".", include_numpy_array],
        extra_compile_args=["-O3", "-Wall"],
        # extra_link_args = ['-g'],
        # libraries = ["dv",],
        )


# Check the availability of arrayobject.h
valid_paths = list(filter(contains_arrayobject_h, sys.path))
if len(valid_paths) == 0:
    print("No paths in the python path contain numpy/arrayobject.h")
    sys.exit(0)

# The base path is by default the first python path with arrayobject.h in it.
include_numpy_array = valid_paths[0]

if len(valid_paths) > 1:
    print("There are several valid include directories"
          "containing numpy/arrayobject.h")
    l = [('%d: %s' % (i+1, valid_paths[i])) for i in
         range(0, len(valid_paths))]
    s = -1
    print('\n'.join(l))
    # Prompt the user with a list of selections.

    # Fix input for Python 3-3
    try:
        input = raw_input
    except NameError:
        pass

    while not (s >= 1 and s <= len(valid_paths)):
        s = input('Selection [default=1]:')
        if s == '':
            s = 1
        else:
            s = int(s)
    include_numpy_array = valid_paths[s-1]

# Add the children directory path suffix to the base path.
include_numpy_array = os.path.join(include_numpy_array, 'numpy', 'core',
                                   'include')

# Need to create a pyOpTools package
extNames = scandir("./")

# and build up the set of Extension objects
extensions = [makeExtension(name) for name in extNames]

setup(name="pyoptools",
      version="0.1.1",
      packages=findpackages("./"),
      scripts=['ipyoptools'],
      # The names from pipy are used, not the deb package names
      requires=['numpy',
                'cython',
                'PyOpenGl',
                'ipython',
                'scipy',
                'six',
                ],
      package_data={
          'pyoptools.raytrace.mat_lib': ['data/*.mat'],
          'pyoptools.raytrace.library': ['Edmund/*.cmp','Thorlabs/*.cmp'],
          },
      author='Ricardo Amezquita Orozco',
      author_email='ramezquitao@cihologramas.com',
      description='Optical ray tracing simulation system',
      license='GPLv3',
      url='https://github.com/cihologramas/pyoptools/',
      download_url='https://github.com/cihologramas/pyoptools/archive/v0.1.1.zip',
      ext_modules=cythonize(extensions),
      cmdclass={'build_ext': build_ext},
      data_files=[("share/doc/pyoptools/examples/basic_course",
                   ["examples/basic_course/00-Introducción.ipynb",
                    "examples/basic_course/03-SimpleComponents.ipynb",
                    "examples/basic_course/05-Autocollimator.ipynb",
                    "examples/basic_course/01-IntroPython.ipynb",
                    "examples/basic_course/04-PredefinedComponents.ipynb",
                    "examples/basic_course/06-GeomWF.ipynb",
                    "examples/basic_course/02-Surfaces.ipynb",
                    "examples/basic_course/04-Simple RayTraces.ipynb",
                    "examples/basic_course/07-SimpleEODs.ipynb"])]
      )
