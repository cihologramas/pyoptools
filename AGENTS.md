# Agent Guidelines for pyOpTools

## Build & Test Commands
- **Build**: `python setup.py build_ext -i` (builds Cython extensions in-place)
- **Install editable**: `pip install -e .` or `pip install -e .[test]`
- **Run all tests**: `pytest` or `python -m pytest`
- **Run single test**: `pytest tests/path/to/test_file.py::test_function_name`
- **Test with coverage**: `pytest --cov=pyoptools`
- **Linting**: `cython-lint pyoptools/**/*.pyx` (for Cython files)

## Pre-commit Hooks
- **Install hooks**: `pre-commit install` (run once to enable automatic checks)
- **Run manually**: `pre-commit run --all-files` (check all files without committing)
- **Run on staged**: `pre-commit run` (check only staged files)
- **Configured hooks**:
  - `cython-lint`: Lints all `.pyx` and `.pxd` files for style and potential issues
  - `double-quote-cython-strings`: Ensures Cython strings use double quotes for consistency
- Pre-commit hooks run automatically on `git commit` after installation
- All hooks must pass before commit is allowed; use `git commit --no-verify` to bypass (not recommended)
- Hooks are defined in [`.pre-commit-config.yaml`](.pre-commit-config.yaml)

## Build System Configuration
- **pyproject.toml**: PEP 517/518 build system declaration, project metadata, dependencies, and tool configuration
- **setup.py**: Required for custom Cython build logic including:
  - Custom `create_extension` function for numpy API macros and platform-specific defines
  - Runtime EIGEN3_INCLUDE_DIR environment variable handling
  - Custom `python setup.py test` command
  - Complex `cythonize()` configuration not expressible declaratively
- Both files must be maintained: pyproject.toml provides modern packaging interface, setup.py contains essential
compiled extension logic

## Code Style
- **Python/Cython**: Follow PEP8 with exceptions: E501 (line length - prefer 88 chars like Black)
- **Indentation**: 4 spaces for .py/.pyx/.pxd files, spaces not tabs
- **Imports**: Standard library first, then third-party (numpy as np, scipy, etc.), then local imports
- **Strings**: Prefer double quotes for consistency with Cython pre-commit hook
- **Type hints**: Encouraged for new code; use `from __future__ import annotations`
- **Docstrings**: Use NumPy-style (Parameters/Returns/Notes sections)
- **Exceptions**: Use specific types (ValueError, RuntimeError) not generic Exception
- **Cython**: Use language_level="3str", apply double-quote-cython-strings pre-commit hook
- **Line endings**: LF (Unix-style), files must end with newline
- **Encoding**: UTF-8 for all Python/Cython files (no encoding headers needed)
- **Naming**: lowercase_with_underscores for functions/variables, PascalCase for classes
- **NO COMMENTS**: Never add comments unless explicitly requested by the user

## Project Structure
- **Cython extensions**: In `pyoptools/**/*.pyx` with corresponding `.pxd` headers
- **Tests**: In `tests/` directory, mirroring package structure; ensure tests exist for new modules
- **Public API**: Use `__all__` to define exports in `__init__.py` files
- **Component library**: `_comp_lib/` contains implementations; `comp_lib.py` is public wrapper
- **Data files**: Material/glass data in `mat_lib/data/`, vendor catalogs in `library/catalogs/`
- **Notebooks**: Development notebooks go in `doc/notebooks/`, not repository root
- **Never commit**: Compiled artifacts (.c, .so, .cpp), __pycache__, .ipynb_checkpoints, Untitled*.ipynb
- **Import conventions**: Internal code may use `_comp_lib`, external users import from `comp_lib`

## Utility Scripts
Scripts in the `scripts/` directory provide maintenance and development tools:

- **clean.py**: Remove Cython-generated build artifacts
  - Deletes `.c`, `.cpp`, and `.so` files generated from `.pyx` sources
  - Usage: `python scripts/clean.py`
  - Run before committing to ensure no compiled artifacts are included
  - Useful when switching branches or cleaning build state

- **autopep8_cython.py**: Auto-format Cython source files
  - Applies autopep8 to all `.pyx` and `.pxd` files recursively
  - Fixes specific PEP8 violations compatible with Cython syntax
  - Usage: `python scripts/autopep8_cython.py`
  - Targets whitespace, indentation, and spacing issues (E114, E127, E201, E202, E211, E221, E228, E231, E251, E261, E262, E265, E271, E272, E301, E302, E303, E304, E701, E702, E703, E711, W291, W293, W391)
  - Note: Pre-commit hooks handle most formatting; use this for bulk cleanup or when hooks are bypassed
