# Utility Scripts

This directory contains maintenance and development utility scripts for pyOpTools.

## clean.py

Remove Cython-generated build artifacts from the project.

**Purpose**: Deletes `.c`, `.cpp`, and `.so` files generated from `.pyx` sources during the build process.

**Usage**:
```bash
python scripts/clean.py
```

**When to use**:
- Before committing to ensure no compiled artifacts are included in version control
- When switching branches to avoid conflicts with different build states
- After encountering build errors that may be caused by stale artifacts
- When performing a clean rebuild of the project

**What it does**:
- Recursively scans the project for all `.pyx` files
- For each `.pyx` file, identifies and removes:
  - Corresponding `.c` and `.cpp` files
  - Associated `.so` (shared object) files with Python version identifiers

## autopep8_cython.py

Auto-format Cython source files to comply with PEP8 style guidelines.

**Purpose**: Applies autopep8 formatting to all `.pyx` and `.pxd` files in the project, fixing specific PEP8 violations that are compatible with Cython syntax.

**Usage**:
```bash
python scripts/autopep8_cython.py
```

**When to use**:
- For bulk cleanup of formatting issues across all Cython files
- When pre-commit hooks have been bypassed and files need formatting
- After manually editing multiple Cython files

**Targeted PEP8 violations**:
- E114, E127: Indentation issues
- E201, E202, E211: Whitespace around brackets
- E221, E228, E231, E251: Spacing around operators
- E261, E262, E265: Comment formatting
- E271, E272: Whitespace around keywords
- E301, E302, E303, E304: Blank line issues
- E701, E702, E703: Multiple statements on one line
- E711: Comparison to None
- W291, W293: Trailing whitespace
- W391: Blank line at end of file

**Note**: The project's pre-commit hooks handle most formatting automatically. This script is primarily useful for bulk operations or when hooks have been bypassed.

## Requirements

Both scripts require a properly configured pyOpTools development environment:
- Python 3.8+
- autopep8 (for autopep8_cython.py)
- Standard library modules (glob, os, re, subprocess)