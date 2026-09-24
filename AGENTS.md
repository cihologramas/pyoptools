# Agent Guidelines for pyOpTools

## Tooling & Execution Pattern
This project uses **uv** for dependency management, virtual environments, and command execution.

- **Fast in-place extension build**:
  ```bash
  uv run python setup.py build_ext --inplace
  ```
  (Compiles changed Cython `.pyx` files incrementally in 1–2 seconds. Avoid building binary wheels during development.)
- **Run all tests**:
  ```bash
  uv run pytest
  ```
- **Run single test**:
  ```bash
  uv run pytest tests/path/to/test_file.py::test_function_name
  ```
- **Linting & Formatting**:
  ```bash
  # Pure Python linting and formatting (Ruff)
  uv run ruff check .
  uv run ruff format .

  # Cython AST linting (.pyx and .pxd)
  uv run cython-lint pyoptools/**/*.pyx pyoptools/**/*.pxd
  ```

## Pre-commit Hooks
- **Install hooks**: `uv run pre-commit install`
- **Run manually**: `uv run pre-commit run --all-files`
- **Configured hooks** (`.pre-commit-config.yaml`):
  - `ruff`: Modern Python linting with automatic fixes
  - `ruff-format`: Fast, standard code formatting
  - `cython-lint`: Lints all `.pyx` and `.pxd` files for Cython style and potential issues
  - `double-quote-cython-strings`: Enforces consistent double quotes in Cython files

## Build System & Dependency Management
- **pyproject.toml**: Declares PEP 517/518 build dependencies (`setuptools`, `Cython`, `eigency`), project dependencies, optional test dependencies, and `[tool.ruff]` / `[tool.cython-lint]` configurations.
- **uv.lock**: Tracked in Git to ensure 100% reproducible environments across developer machines and CI.
- **setup.py**: Houses custom Cython extension setup, Eigen header discovery via `eigency`, and NumPy 2.x macro configuration (`NPY_NO_DEPRECATED_API`).

## Code Style & Guidelines
- **Python Formatting**: Target Python 3.10+, line length 88 (configured via Ruff in `pyproject.toml`).
- **Imports**: Standard library first, third-party second, local imports third. Re-export modules (`__init__.py`, `comp_lib.py`, `all.py`) use `__all__` or relative imports.
- **Cython**: Use `language_level="3str"`. Always prefer C-scalar comparisons and typed variables in inner ray-tracing loops (`System.propagate_ray`, `Component.distance`, `Component.propagate`) rather than allocating intermediate Python lists or invoking NumPy arrays in hot loops.
- **Documentation**: Use NumPy-style docstrings (`Parameters`, `Returns`, `Notes`). Retain existing documentation and add comments where algorithm rationale is non-obvious.

## Utility Scripts
Located in `scripts/`:
- **clean.py**: Removes Cython-generated build artifacts (`.c`, `.cpp`, `.so`).
  ```bash
  uv run python scripts/clean.py
  ```
