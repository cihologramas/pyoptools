#!/usr/bin/env python
import sys
import subprocess
import platform
from setuptools import setup, Command
from Cython.Build import cythonize
from Cython.Build.Dependencies import default_create_extension

def get_eigen_include():
    system = platform.system()

    if system == "Linux":
        return "/usr/include/eigen3"  # Example path for Ubuntu 24.04

    elif system == "Windows":
        return "C:\\path\\to\\eigen"  # Example path for Windows 2019

    elif system == "Darwin":
        # Check macOS version here if needed
        return "/usr/local/include/eigen3"  # Example path for macOS 11
    
    else:
        raise ValueError("Unsupported operating system")

eigen_include_path = get_eigen_include()


def create_extension(template, kwds: dict):
    define_macros = kwds.get("define_macros", [])

    # Use the new numpy API and remove all the compilation warnings
    define_macros.append(("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"))

    if sys.platform in ("darwin", "win32"):
        define_macros.append(("CYTHON_INLINE", ""))

    kwds["define_macros"] = define_macros
    return default_create_extension(template, kwds)


# Custom command to build extensions and run tests
class TestCommand(Command):
    description = "Build extensions and run tests."
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        # Build extensions
        self.run_command("build_ext")
        # Install the package in editable mode
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", ".[test]"])
        # Run tests using pytest
        errno = subprocess.call([sys.executable, "-m", "pytest"])
        raise SystemExit(errno)


if __name__ == "__main__":
    # allow setup.py to run from another directory
    setup(
        ext_modules=cythonize(
            "pyoptools/**/*.pyx",
            create_extension=create_extension,
            language_level="3str",
        ),
        include_dirs=[eigen_include_path],
        use_scm_version=True,
        include_package_data=True,
        cmdclass={"test": TestCommand},
        setup_requires=["setuptools_scm", "Cython"]
    )
