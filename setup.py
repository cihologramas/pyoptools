#!/usr/bin/env python
import sys
import subprocess
from setuptools import setup
from Cython.Build import cythonize
from Cython.Build.Dependencies import default_create_extension

def get_eigen_include():
    try:
        # Run pkg-config to get the include flags
        include_path = subprocess.check_output(["pkg-config", "--cflags-only-I", "eigen3"])
        include_path = include_path.decode("utf-8").strip().replace("-I", "")
        return include_path
    except subprocess.CalledProcessError:
        print("Could not find Eigen include path using pkg-config. Make sure Eigen is installed and pkg-config is configured correctly.")
        return None

eigen_include_path = get_eigen_include()


def create_extension(template, kwds: dict):
    define_macros = kwds.get("define_macros", [])

    #Use the new numpy API and remove all the compilation warnings
    define_macros.append(("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"))

    if sys.platform in ("darwin", "win32"):
        define_macros.append(("CYTHON_INLINE", ""))

    kwds["define_macros"] = define_macros
    return default_create_extension(template, kwds)


if __name__ == "__main__":
    # allow setup.py to run from another directory
    setup(
        ext_modules=cythonize("pyoptools/**/*.pyx", create_extension=create_extension,
                              language_level="3str"),
        include_dirs=[eigen_include_path],
        use_scm_version=True,
        include_package_data=True,
    )
