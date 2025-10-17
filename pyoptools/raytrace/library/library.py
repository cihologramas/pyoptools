from pyoptools.raytrace._comp_lib.optic_factory import optic_factory

import orjson
import warnings

from importlib_resources import files
from pathlib import Path


# This class overrides the module, to provide direct attribute and item access
class LibraryModule:
    """The optics library provides access to optical components from
    the catalogs of well-known optical suppliers or user-supplied libraries.

    Currently, catalogs from Thorlabs [thorlabs.com] and Edmund
    [edmundoptics.com] are available.

    Optical elements can be retrieved by part number using dictionary-style
    access from the library module i.e:

    `library['LB1862']` returns a SphericalLens based on an optic from
    Thorlabs.

    optionally, the specific supplier catalog can be specified via attribute
    access. This can be used in the rare case that different catalogs have
    parts with matching part number, or to speed up lookup time. I.e.

    `library.Thorlabs['LB1862']`

    It is also possible to retrieve a dictionary-descriptor for an optic
    without instantiating the optical component instance. This is done with
    the `.descriptor()` method. Descriptors can be used to instantiate
    optical components using the `component.optic_factory`.

    The library module is iterable using `.items()`. This dictionary-style
    iterator returns component descriptor. A common usage is to search for
    optical components matching specific requirements.

    User supplied catalogs can be added using the `.add()` method, which
    requires the path (string or Pathlib.Path) to a .json optic catalog file.
    The user supplied catalog can then be accessed by attribute using the
    given filename.
    """

    def __init__(self):
        self._user_libraries = []
        self.dp = files("pyoptools.raytrace.library").joinpath("catalogs")

    def __getattr__(self, name: str):
        p = self.dp / (name + ".json")
        if p.exists():
            return OpticCatalog(p)
        else:
            raise AttributeError(f"Catalog {name} not available in library.")

    def items(self):
        """Generates a dictionary iterator which allows iteration through
        all optical components descriptors in all libraries.
        Note : to instantiate an optical component given a descriptor,
        use component.optic_factory(**descriptor).
        """
        for jfp in self._json_files():
            oc = OpticCatalog(jfp)
            for k, v in oc.items():
                yield (k, v)

    def parts(self):
        """Returns a list of all available optical component by part number."""
        return [k for k, _ in self.items()]

    def descriptor(self, part):
        """Return the dictionary-descriptor for optic with given part number."""
        for jf in self._json_files():
            oc = OpticCatalog(jf)
            try:
                return oc.descriptor(part)
            except KeyError:
                continue
        raise KeyError(f"{part} not found in any available optics catalog.")

    def __getitem__(self, part: str):
        return optic_factory(**self.descriptor(part))

    def get(self, part: str):
        """Deprecated method. Simply use dictionary style access to get
        an optic given the part number.
        """
        warnings.warn(
            "This method is deprecated, you can use dictionary-style access " "instead",
            DeprecationWarning,
        )
        return optic_factory(**self.descriptor(part))

    def add(self, userlib):
        """Add a user-defined .json file containing optical component
        descriptors to the library.
        userlib : string path to the .json file or pathlib.Path
        """

        if isinstance(userlib, str):
            p = Path(userlib)
        elif isinstance(userlib, Path):
            p = userlib
        else:
            raise ValueError(f"{userlib} not a string or Pathlib.Path.")

        if p.suffix != ".json" or not p.exists():
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
        return list(self.dp.glob("*.json")) + self._user_libraries

    def catalogs(self):
        return [catalog.stem for catalog in self._json_files()]


class OpticCatalog:
    def __init__(self, catalog_path):
        self.catalog_path = catalog_path
        with open(catalog_path, "rb") as f:
            self.catalog_data = orjson.loads(f.read())

    def items(self):
        return self.catalog_data.items()

    def descriptor(self, part):
        """Return the dictionary-descriptor for optic with given part number."""
        try:
            return self.catalog_data[part]
        except KeyError:
            raise KeyError(f"{part} not in optic catalog {self.catalog_path.name}")

    def __getitem__(self, part):
        return optic_factory(**self.descriptor(part))

    def get(self, part):
        warnings.warn(
            "This method is deprecated, you can use dictionary-style access instead",
            DeprecationWarning,
        )
        return optic_factory(**self.descriptor(part))

    def parts(self):
        """Returns a list of all available optical component by part number."""
        return list(self.catalog_data.keys())
