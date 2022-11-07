from pyoptools.raytrace._comp_lib.optic_factory import optic_factory

import ijson
from pathlib import Path
import sys

# This class overrides the module, to provide direct attribute and item access
class LibraryModule:

    def __init__(self):
        self._user_libraries = []

    def __getattr__(self, name: str):
        p = Path(__file__).parent/(name+'.json')
        if p.exists():
            return OpticCatalog(p)
        else:
            raise AttributeError(f"Catalog {name} not available in library.")

    def __getitem__(self, name: str):
        for jf in self._json_files():
            oc = OpticCatalog(jf)
            try:
                return oc[name]
            except KeyError:
                continue
        raise KeyError(f"{name} not found in any available optics catalog.")

    def items(self):
        """Generates a dictionary iterator which allows iteration through
        all optical components descriptors in all libraries.
        Note : to instantiate an optical component given a descriptor,
        use optic_factory(**descrioptor).
        """
        for jfp in self._json_files():
            oc = OpticCatalog(jfp)
            for k, v in oc.items():
                yield (k, v)

    def get(self, part_number: str):
        """Depricated method. Simply use dictionay style access to get
        an optic given the part number.
        """
        return optic_factory(**self[part_number])

    def add(self, userlib: str):
        """Add a user-defined .json file containing optical component
        descriptors to the library.
        userlib : string path to the .json file.
        """
        p = Path(userlib)

        if (p.suffix != '.json' or not p.exists()):
            raise FileNotFoundError(f"{userlib} not a valid .json catalog.")

        if p.name in [x.name for x in self._user_libraries]:
            raise ValueError(f"{p.name} already loaded.")

        self._user_libraries.append(p)

    def remove(self, libname: str):
        """Remove a user-defined optical component library."""
        for lib in self._user_libraries:
            if lib.name == libname:
                self._user_libraries.remove(lib)

    def _json_files(self):
        return list(Path(__file__).parent.glob('*.json')) + self._user_libraries

sys.modules[__name__] = LibraryModule()

class OpticCatalog:
    def __init__(self, catalog_path):
        self.catalog_path = catalog_path
        self.jf = catalog_path.open()

    def items(self):
        gen = ijson.kvitems(self.jf, '', use_float=True)
        self.jf.seek(0)
        return gen

    def descriptor(self, part):
        """Returns a dictionay-descriptor for a given optical component
        """
        results = list(ijson.items(self.jf, part, use_float=True))
        if results:
            self.jf.seek(0)
            return results[0]

        self.jf.seek(0)
        raise KeyError(f"{part} not in optic catalog {self.catalog_path.name}")

    def __getitem__(self, part):
        return optic_factory(**self.descriptor(part))

    def get(self, part_number):
        """Depricated method. Simply use dictionay style access to get
        an optic given the part number.
        """
        return optic_factory(**self.descriptor(part_number))

    def parts(self):
        """Deprecated method. Iterate using library.items() to get all
           the available optical components."""

        raise NotImplementedError("Method deprecated, it is now possible "
                                  "to simply iterate through all optics "
                                  "descriptors library.items().")

    def __del__(self):
        #print('closing ', self.jf)
        self.jf.close()

