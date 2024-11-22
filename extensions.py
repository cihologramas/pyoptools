import sys
from Cython.Build.Dependencies import default_create_extension


def create_extension(template, kwds: dict):
    define_macros = kwds.get("define_macros", [])

    # Use the new numpy API and remove all the compilation warnings
    define_macros.append(("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"))

    if sys.platform in ("darwin", "win32"):
        define_macros.append(("CYTHON_INLINE", ""))

    kwds["define_macros"] = define_macros
    return default_create_extension(template, kwds)
