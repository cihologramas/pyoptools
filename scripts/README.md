# Utility Scripts

This directory contains maintenance and development utility scripts for pyOpTools.

## clean.py

Remove Cython-generated build artifacts from the project.

**Purpose**: Deletes `.c`, `.cpp`, and `.so` files generated from `.pyx` sources during the build process.

**Usage**:
```bash
uv run python scripts/clean.py
```

**When to use**:
- Before committing to ensure no compiled artifacts are included in version control
- When switching branches to avoid conflicts with different build states
- After encountering build errors that may be caused by stale artifacts
- When performing a clean rebuild of the project
