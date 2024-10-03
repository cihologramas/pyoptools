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


# This class overrides the module definition
class MaterialLibrary:
    """The material library provides access to refractive index data
    for common optical glasses and generic organic and inorganic compounds.

    All data is from www.refractiveindex.info. Please follow their guidelines
    for citation.

    `raytrace.mat_lib.material` instances for optical glasses can be retrieved
    directly using dictionary-style access. For example:

    `material['N-BK7']`

    In some cases, materials from different manufacturer catalogs may share
    the same identifier, in this case a KeyError will be raised.
    This can be resolved using attribute access given the name of
    the manufacturer catalog. For example:

    ```
    m = material['SF11']
        KeyError: "Multiple matches for glass type SF11. Use one of:
        material.hikari['SF11'] or material.schott['SF11']"
    ```

    Compounds are also available through the attributes `organic` and
    `inorganic`. Note that in refractiveindex.info, inorganic is named _main_.

    Often, multiple data sources are available for a given compound.
    These can be disambiguated using a trailing : followed by the name of the
    reference. For example:

    ```
    silver = material.inorganic['Ag:Yang']
    ```
    If the reference is omitted then the first available reference will be
    returned.

    Some glasses and compounds have common abbreviations. These can also be
    used, through the aliases system, which takes precedent in material
    retrieval. For example:

    ```
    m = material['SODA_LIME']
    ```

    the full dictionary of available in the attribute
    `material.aliases`.
    """

    def __init__(self, prefix=None):

        self.prefix = prefix
        self._cache = {}

        self.dp = files("pyoptools.raytrace.mat_lib").joinpath("data")

        if self.prefix is None:
            self.glass_path = self.dp / "glass"
        else:
            self.glass_path = self.dp / "glass" / prefix

        with (self.dp / "aliases.json").open() as af:
            self.aliases = json.load(af)

        self._compound_lib_names = ["organic", "inorganic"]

    def _material_factory(self, name, mat_path):
        "Builds and caches a material given path to yml file"
        mat = from_yml(mat_path)
        self._cache[name] = mat
        return mat

    def _get_in_aliases(self, name):
        if self.prefix is None:
            return self._get_in_aliases_base(name)
        else:
            return self._get_in_aliases_catalog(name)

    def _get_in_aliases_catalog(self, name):
        "Restricted version of get aliases for when accessing in a catalog"

        if name in self.aliases and self.aliases[name]["library"] == self.prefix:
            a = self.aliases[name]
            ap = self.dp / "glass" / a["library"] / f"{a['material']}.yml"
            return self._material_factory(name, ap)
        else:
            raise KeyError

    def _get_in_aliases_base(self, name):

        if name in self.aliases:
            a = self.aliases[name]
            if a["type"] == "organic":
                return self.organic[a["material"]]
            elif a["type"] == "inorganic":
                return self.inorganic[a["material"]]
            elif a["type"] == "glass":
                ap = self.dp / "glass" / a["library"] / f"{a['material']}.yml"
                return self._material_factory(name, ap)
        else:
            raise KeyError

    def __getitem__(self, name: str):

        if name in self._cache:
            return self._cache[name]

        # check in aliases, handling special case of a compound
        try:
            return self._get_in_aliases(name)
        except KeyError:
            pass

        # find in glasses
        matches = list(self.glass_path.glob(f"**/{name}.yml"))
        if len(matches) > 1:
            warning = f"Multiple matches for glass type {name}. " f"Use one of: "
            for m in matches:
                warning += f"material.{m.parts[-2]}['{name}'] or "
            raise KeyError(warning[:-4])

        if matches:
            return self._material_factory(name, matches[0])
        else:
            raise KeyError(f"Material {name} not found.")

    def __contains__(self, name):
        try:
            self[name]  # This will call __getitem__
            return True
        except KeyError:
            return False

    def find_material(self, search, printout=True, exact=False, unalias=False):
        """Search for a material where the string _search_ is found in the
        name or description. For example:

        If printout is True, prints a list of the commands used to access
        the matching materials.

        If exact is True, return only the items where the material name is an
        exact match of the search string.

        If unalias is True, the materials found in the alises list, will be
        translated to the real material.

        Returns a list of the keys which can be used to retrieve candidate
        materials.

        Example finding all SF-11 type glasses:
        ```
        material.find_material('SF11')
            material.hikari['J-SF11']
            material.hikari['E-SF11']
            material.hikari['SF11']
            material.schott['N-SF11']
            material.schott['SF11']
        ```

        Note that for compounds, the common name can be searched as well,
        for example:

        ```
        material.find_material('styrene')
            material.organic['C8H8:Sultanova']
            material.organic['C8H8:Myers']
        ```
        """

        results = []
        # Check in aliases:
        for k, v in self.aliases.items():
            if search in k or search in v["material"]:
                results.append((None, k, None))

        # Check in glasses:
        matches = list(self.glass_path.glob(f"**/*{search}*"))
        for m in matches:
            catalog = m.parts[-2]
            name = m.parts[-1].split(".")[0]
            results.append((catalog, name, None))

        # Check in compounds:
        for cln in self._compound_lib_names:
            r = getattr(self, cln).find_material(search, printout=False)
            if r is not None:
                results += r

        if exact:
            results = list(filter(lambda mat_item: mat_item[1] == search, results))
        if unalias:
            unaliased = []
            for r in results:
                if r[0] is None:
                    realmat = self.aliases[r[1]]

                    if realmat["type"] != "glass":
                        raise ValueError(
                            "Don't know how to handle {}".format(realmat["type"])
                        )

                    unaliased.append((realmat["library"], realmat["material"], None))
                else:
                    unaliased.append(r)
            results = unaliased

        if printout:
            for r in results:
                catalog, name, ref = r
                if catalog is None:
                    print(f"material['{name}']")
                if ref is None:
                    print(f"material.{catalog}['{name}']")
                else:
                    print(f"material.{catalog}['{name}:{ref}']")

        return results

    def get_from(self, name: str, libs: str, check_aliases=True):
        """Finds a glass type located in a specific manufacturer library,
        or space-separate list of libraries.

        This is mostly used internally because many optical component
        definitions specify the glass catalogs in this way.

        Note that this method is only available from the base material
        module.

        name : name of the glass to find
        libs : A name of a manufacturer library, e.g. 'schott', or a list
            of space-separate manufacturer names, in which case they will be
            searched in the order given. Insensitive to capitalization.
        check_aliases : if True (default), will check within material
            aliases after checking in specified libraries.
            Note that aliases can include compounds in addition to glasses.

        Raise KeyError if no material found in the given libraries.
        """

        if self.prefix is not None:
            raise Exception("get_from only available on base library.")

        for libname in libs.split(" "):
            libname = libname.lower()
            if (self.glass_path / libname).is_dir():
                try:
                    return MaterialLibrary(prefix=libname)[name]
                except AttributeError:
                    # print('No library')
                    pass
                except KeyError:
                    # print('Not in library')
                    pass

        warning = f"Material {name} not found in any of {libs.split()}."
        if check_aliases:
            warning = warning[:-1] + " or aliases."
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

        if name == "inorganic" or name == "organic":
            return CompoundLibrary(self.dp / name, aliases=self.aliases)
        if (self.glass_path / name).is_dir():
            return MaterialLibrary(prefix=name)
        else:
            raise AttributeError(f"Material {name} not found.")

    def get_glass_libraries(self):
        """Returns a list of strings containing the names of the glass
        libraries defined in PyOpTools"""

        libpaths = list(self.glass_path.glob("*"))

        return ["aliases"] + [lib_path.parts[-1] for lib_path in libpaths]

    def get_glass_materials_from_library(self, libname):
        """Returns a list of strings containing the names of the materials
        defined in the library 'libname'
        """

        if libname == "aliases":
            return list(self.aliases.keys())

        matfiles = self.glass_path.glob(f"{libname}/*.yml")
        return [mat.stem for mat in matfiles]


class CompoundLibrary:
    def __init__(self, lib_path, aliases=None):
        self.lib_path = lib_path
        self.aliases = aliases

        # Populate dict of all available compound dirs
        # Organic compounds have long names with space so just use the
        # first part as the identifier
        self.compound_dirs = {}
        for i in lib_path.iterdir():
            if i.is_dir():
                self.compound_dirs[i.name.split(" ")[0]] = i

        self._cache = {}

    def _material_factory(self, name, mat_path):
        "Builds and caches a material given path to yml file"
        mat = from_yml(mat_path)
        self._cache[name] = mat
        return mat

    def __getitem__(self, identifier: str):

        # Check in cache
        if identifier in self._cache:
            return self._cache[identifier]

        # Check in aliases
        if (
            self.aliases is not None
            and identifier in self.aliases
            and self.aliases[identifier]["type"] == self.lib_path.name
        ):
            identifier = self.aliases[identifier]["material"]

        tokens = identifier.split(":")
        tokens = [t for t in tokens if t != ""]
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
            compound_file = list(cd.glob("*.yml"))[0]
        else:
            compound_file = cd / f"{reference}.yml"

        if compound_file.exists():
            return self._material_factory(identifier, compound_file)
        else:
            raise KeyError(f"Compound reference not found for {identifier}.")

    def find_material(self, search, printout=True):

        results = []
        for i in self.lib_path.iterdir():
            if i.is_dir():
                if search in i.name:
                    for ref in i.glob("*.yml"):
                        results.append(
                            (
                                self.lib_path.name,
                                i.name.split(" ")[0],
                                ref.name.split(".")[0],
                            )
                        )

        return results


sys.modules[__name__] = MaterialLibrary()
