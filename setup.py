import os
import subprocess
import sys

from Cython.Build import cythonize
from Cython.Build.Dependencies import default_create_extension
from setuptools import Command, setup


def get_eigen_include_path():
    # 1. Environment variable
    env_path = os.environ.get("EIGEN3_INCLUDE_DIR")
    if env_path and (
        os.path.exists(os.path.join(env_path, "Eigen", "Dense"))
        or os.path.exists(env_path)
    ):
        return env_path

    # 2. Python package providing Eigen headers (managed via uv/pip)
    try:
        import eigency

        for inc in eigency.get_includes():
            if os.path.exists(os.path.join(inc, "Eigen", "Dense")):
                return inc
    except ImportError:
        pass

    # 3. Standard system locations
    for candidate in [
        "/usr/include/eigen3",
        "/usr/local/include/eigen3",
        "/opt/homebrew/include/eigen3",
    ]:
        if os.path.exists(os.path.join(candidate, "Eigen", "Dense")) or os.path.exists(
            candidate
        ):
            return candidate

    return "/usr/include/eigen3"


eigen_include_path = get_eigen_include_path()


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
    )
