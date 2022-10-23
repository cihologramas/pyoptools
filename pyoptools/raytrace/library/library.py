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
        raise KeyError(f"{item} not in optic any available optic catalog.")

    def items(self):
        for jfp in self._json_files():
            oc = OpticCatalog(jfp)
            for k, v in oc.items():
                yield (k, v)

    def get(self, part_number: str):
        return optic_factory(**self[part_number])

    def add(self, userlib: str):
        p = Path(userlib)

        if (p.suffix != '.json' or not p.exists()):
            raise FileNotFoundError(f"{userlib} not a valid .json catalog.")

        if p.name in [x.name for x in self._user_libraries]:
            raise ValueError(f"{p.name} already loaded.")

        self._user_libraries.append(p)

    def remove(self, libname: str):
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

    def __getitem__(self, item):
        results = list(ijson.items(self.jf, item, use_float=True))
        if results:
            self.jf.seek(0)
            return results[0]

        self.jf.seek(0)
        raise KeyError(f"{item} not in optic catalog {self.catalog_path.name}")

    def get(self, part_number):
        return optic_factory(**self[part_number])

    def parts(self):
        pass

    def __del__(self):
        #print('closing ', self.jf)
        self.jf.close()

