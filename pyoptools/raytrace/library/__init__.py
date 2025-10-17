"""Optical component library module.

This module uses a module replacement pattern to provide a convenient API for
accessing optical components from vendor catalogs. The module imports the
LibraryModule class, creates an instance, and then replaces itself in
sys.modules with that instance.

This allows users to interact with the library as if it were an instance,
enabling dictionary-style access (library['part_number']), attribute access
(library.Thorlabs['part_number']), and iteration (library.items()).

The replacement must happen here in __init__.py after successfully importing
the classes from library.py. If the replacement were done in library.py itself,
the import would fail because Python wouldn't be able to find the class
definitions to import.

To support different import patterns used throughout the codebase, we expose:
- library: the LibraryModule instance (supports: from pyoptools.raytrace.library import library)
- LibraryModule: the class itself (supports: from pyoptools.raytrace.library import LibraryModule)
- OpticCatalog: the class itself (supports: from pyoptools.raytrace.library import OpticCatalog)

All three are set as attributes on the instance before the module replacement.
"""
from .library import LibraryModule, OpticCatalog
import sys

_library_instance = LibraryModule()

_library_instance.library = _library_instance
_library_instance.LibraryModule = LibraryModule
_library_instance.OpticCatalog = OpticCatalog

sys.modules[__name__] = _library_instance

__all__ = ["LibraryModule", "OpticCatalog", "library"]
