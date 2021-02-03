#!/usr/bin/env python
# -*- coding: utf-8-*-

# ------------------------------------------------------------------------------
# Copyright (c) 2007-2021, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Material definition helper class
# Symbols Defined: Material
#
#
# ------------------------------------------------------------------------------

'''
Material class definition, and helper functions used to load the
constants of the dispersion formula to be used in the calculation of
the refraction index.

It uses the database from  https://refractiveindex.info
'''

from os import walk
from os.path import join, expanduser, relpath
from pkg_resources import resource_filename
from .mat_eq import from_yml, ModelNotImplemented
from configparser import ConfigParser

mat_config = resource_filename("pyoptools.raytrace.mat_lib", 'data')

# mat_config = "../mat_lib/data/"

# Get library names from the system

# Get data from glass folder

libnames = []
libpath = join(mat_config, "glass")
for (dirpath, dirnames, filenames) in walk(libpath):
    library = relpath(dirpath, libpath)

    # Exclude some names that are not libraries
    if library in [".", ]:
        continue

    libnames.append((relpath(dirpath, libpath)).replace("/", "_"))

# Get data from main folder

mainlibnames = []
mainlibpath = join(mat_config, "main")
for (dirpath, dirnames, filenames) in walk(mainlibpath):
    library = relpath(dirpath, libpath)

    # Exclude some names that are not libraries
    if library in [".", ]:
        continue

    mainlibnames.append((relpath(dirpath, mainlibpath)).replace("/", "_"))


# Get library names from the user home

homelibpath = join(expanduser("~"), ".pyoptools", "material", "glass")

homelibnames = []
for (dirpath, dirnames, filenames) in walk(homelibpath):
    library = relpath(dirpath, homelibpath)

    if library in [".", ]:
        continue

    homelibnames.append((relpath(dirpath, homelibpath)).replace("/", "_"))

# Create the materials dictionary

# Note: If a home library has the same name as a system library, all the
# glasses defined will be merged in the same library

libset = list(set(libnames+mainlibnames+homelibnames))
liblist = []

for libname in libset:
    # Create the dictionaries where the materials will be saved. One dictionary
    # per library,
    globals()[libname] = {}
    liblist.append((libname, globals()[libname]))


# Fill the dictionaries with the current materials system wide, and then with
# the materials defined in the home of the user

for npath in [libpath, mainlibpath, homelibpath]:

    for (dirpath, dirnames, filenames) in walk(npath):
        library = (relpath(dirpath, npath)).replace("/","_")
        # Exclude some names that are not libraries
        if library in [".", ]:
            continue

        for name in filenames:
            try:
                matname = name.split(".")[0]
                globals()[library][matname] = from_yml(join(dirpath, name))
            except ModelNotImplemented:
                continue

# Create the aliases material library. It will read the information from the
# aliases.cfg file

aliases_path = join(mat_config,"aliases.cfg")


globals()["aliases"] = {}
config = ConfigParser()
config.read(aliases_path)

for i in config:
    if i == "DEFAULT": continue
    libr = config[i]["library"]
    mate = config[i]["material"]
    globals()["aliases"][i] =  globals()[libr][mate]

liblist.append(("aliases", globals()["aliases"]))

def find_material(material):
    """Search for a material in all the libraries

    This function prints all the libraries that contain the material

    Arguments:

    material
        String with the material name
    """
    retv=[]
    for libn, _ in liblist:
       if material in globals()[libn]:
          retv.append(libn)
    return retv


def get_material(material):
    """Search for a material in all the libraries

    This function search in all the material libraries, and return the
    first instance found, that matches the name of the material requested.
    If no material found, returns None
    Arguments:

    material
        String with the material name
    """
    for libn, _ in liblist:
        tdict=globals()[libn]
        if material in tdict:
            return tdict[material]
    print (material, " not found")
    raise KeyError



def mat_list():
    for libn, _ in liblist:
        tdict = globals()[libn]
        print(libn,  tdict.keys())
