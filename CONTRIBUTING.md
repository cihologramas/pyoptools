# Contributing to pyOpTools

Thank you for your interest in contributing to pyOpTools! This document provides guidelines and instructions for setting up your development environment and contributing to the project.

## Development Setup

### Prerequisites

- Python 3.8 or higher
- C++ compiler (for Cython extensions)
- Git

### Installation

1. Clone the repository:
```bash
git clone https://github.com/cihologramas/pyoptools.git
cd pyoptools
```

2. Install the package in editable mode with development dependencies:
```bash
pip install -e .[test]
```

This command will automatically build the Cython extensions and install the package in development mode.

### Rebuilding Extensions

If you modify Cython source files (`.pyx` or `.pxd`), you can rebuild extensions without reinstalling:
```bash
python setup.py build_ext -i
```

Alternatively, re-run the editable install:
```bash
pip install -e .
```

## Pre-commit Hooks

We use pre-commit hooks to ensure code quality and consistency. These hooks automatically check your code before each commit.

### Setup

1. Install pre-commit hooks (one-time setup):
```bash
pre-commit install
```

### Usage

After installation, pre-commit hooks will run automatically on `git commit`. The hooks will:
- Lint all Cython files (`.pyx` and `.pxd`) for style and potential issues
- Ensure all Cython strings use double quotes for consistency

To run hooks manually without committing:
```bash
pre-commit run --all-files
```

To run hooks only on staged files:
```bash
pre-commit run
```

### Configured Hooks

The project uses the following hooks (defined in [`.pre-commit-config.yaml`](.pre-commit-config.yaml)):
- **cython-lint**: Lints all `.pyx` and `.pxd` files
- **double-quote-cython-strings**: Enforces double quotes in Cython strings

All hooks must pass before a commit is allowed. To bypass hooks (not recommended):
```bash
git commit --no-verify
```

## Running Tests

Run all tests:
```bash
pytest
```

Run tests with coverage:
```bash
pytest --cov=pyoptools
```

Run a specific test:
```bash
pytest tests/path/to/test_file.py::test_function_name
```

## Code Style

- Follow PEP8 with line length preference of 88 characters
- Use 4 spaces for indentation (no tabs)
- Prefer double quotes for strings (enforced by pre-commit hooks)
- Use NumPy-style docstrings
- Import order: standard library, third-party, local imports

For more detailed code style guidelines, see [`AGENTS.md`](AGENTS.md).

## Utility Scripts

The [`scripts/`](scripts/) directory contains helpful development tools:

- **clean.py**: Remove Cython-generated build artifacts
  ```bash
  python scripts/clean.py
  ```

- **autopep8_cython.py**: Auto-format Cython source files
  ```bash
  python scripts/autopep8_cython.py
  ```

## Making Changes

1. Create a new branch for your changes:
```bash
git checkout -b feature-or-fix-name
```

2. Make your changes and ensure tests pass
3. Commit your changes (pre-commit hooks will run automatically)
4. Push your branch and create a pull request

## Additional Resources

- **Build and test commands**: See [`AGENTS.md`](AGENTS.md) for comprehensive build instructions
- **Project structure**: See [`AGENTS.md`](AGENTS.md) for detailed project organization
- **Documentation**: Available at https://pyoptools.readthedocs.io/

## Questions?

If you have questions or need help, please:
- Open an issue on GitHub
- Check the existing documentation
- Contact the development team

Thank you for contributing to pyOpTools!