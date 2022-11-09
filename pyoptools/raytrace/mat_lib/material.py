"""
Material class definition, and helper functions used to load the
constants of the dispersion formula to be used in the calculation of
the refraction index.

It uses the database from  https://refractiveindex.info
"""

from .mat_eq import from_yml, ModelNotImplemented

import sys
import json

from pathlib import Path
from importlib_resources import files

class MaterialLibrary:

    def __init__(self, prefix = None):

        self.prefix = prefix
        self._cache = {}

        self.dp = files('pyoptools.raytrace.mat_lib').joinpath('data')

        if self.prefix is None:
            self.glass_path = self.dp/'glass'
        else:
            self.glass_path = self.dp/'glass'/prefix

        with (self.dp/'aliases.json').open() as af:
            self.ailises = json.load(af)

    def _material_factory(self, name, mat_path):
        "Builds and caches a material given path to yml file"
        mat = from_yml(mat_path)
        self._cache[name] = mat
        return mat

    def _get_in_aliases(self, name):

        if name in self.ailises:
            a = self.ailises[name]
            if a['type'] == 'organic':
                return self.organic[a['material']]
            elif a['type'] == 'inorganic':
                return self.inorganic[a['material']]
            elif a['type'] == 'glass':
                ap = self.dp/'glass'/a['library']/f"{a['material']}.yml"
                return self._material_factory(name, ap)
        else:
            raise KeyError

    def __getitem__(self, name: str):

        if name in self._cache:
            return self._cache[name]

        # check in aliases, handling special case of a compond
        try:
            return self._get_in_aliases(name)
        except KeyError:
            pass

        # find in glasses
        matches = list(self.glass_path.glob(f"**/{name}.yml"))
        if len(matches) > 1:
            warning = (f"Multiple matches for glass type {name}. "
                       f"Use one of: ")
            for m in matches:
                warning += f"material.{m.parts[-2]}['{name}'] or "
            raise KeyError(warning[:-4])

        if matches:
            return self._material_factory(name, matches[0])
        else:
            raise KeyError(f"Material {name} not found.")

    def get_from(self, name: str, libs: str, check_aliases = True):
        """Finds a glass type located in a specific manufacturer library
        or space-seperate list of libraries,

        name : name of the glass to find
        libs : A name of a manufacturer library, e.g. 'schott', or a list
            of space-seperate manufacturer names, in which case they will be
            searched in the order given. Insensitive to capitalization.
        check_aliases : if True (default), will check within material
            aliases after checking in specified libraries.
            Note that aliases can include compounds in addition to glasses.

        Raise Exception if no material found in the given libraries.
        """

        if self.prefix is not None:
            raise Exception('get_from only availabe on base library.')

        for libname in libs.split(' '):
            libname = libname.lower()
            #print(f"checking {libname}")
            if (self.glass_path/libname).is_dir():
                try:
                    return MaterialLibrary(prefix = libname)[name]
                except AttributeError:
                    #print('No library')
                    pass
                except KeyError:
                    #print('Not in library')
                    pass

        warning = f"Material {name} not found in any of {libs.split()}."
        if check_aliases:
            warning = warning[:-1]+" or aliases."
            try:
                return self._get_in_aliases(name)
            except KeyError:
                pass
        else:
            warning = f"Material {name} not found in any of {libs.split()}."

        raise KeyError(warning)

    def __getattr__(self, name: str):
        # Guard for if instantiated as a sub-module
        if self.prefix is not None:
            raise AttributeError()

        if name == 'inorganic' or name == 'organic':
            return CompoundLibrary(self.dp/name)
        if (self.glass_path/name).is_dir():
            return MaterialLibrary(prefix = name)
        else:
            raise AttributeError(f"Material {name} not found.")

class CompoundLibrary:
    def __init__(self, lib_path):
        self.lib_path = lib_path

        # Populate dict of all available compound dirs
        # Organic compounds have long names with space so just use the
        # first part as the identifier
        self.compound_dirs = {}
        for i in lib_path.iterdir():
            if i.is_dir():
                self.compound_dirs[i.name.split(' ')[0]] = i

        self._cache = {}

    def _material_factory(self, name, mat_path):
        "Builds and caches a material given path to yml file"
        mat = from_yml(mat_path)
        self._cache[name] = mat
        return mat

    def __getitem__(self, identifier: str):

        if identifier in self._cache:
            return self._cache[identifier]

        tokens = identifier.split(':')
        tokens = [t for t in tokens if t != '']
        compound = tokens[0]
        try:
            reference = tokens[1]
        except IndexError:
            reference = None

        # Get the compound dir
        try:
            cd = self.compound_dirs[compound]
        except KeyError:
            raise KeyError(f"Compound {identifier} not found.")

        if reference is None:
            # if reference unspecified, default to the first item
            compound_file = list(cd.glob('*.yml'))[0]
        else:
            compound_file = cd/f"{reference}.yml"

        if compound_file.exists():
            print('cf ', compound_file)
            return self._material_factory(identifier, compound_file)
        else:
            raise KeyError(f"Compound reference not found for {identifier}.")

sys.modules[__name__] = MaterialLibrary()

