Development Guide
=================

Development Setup
-----------------

We recommend using `uv <https://docs.astral.sh/uv/>`_ for fast, deterministic development environments:

.. code-block:: bash

    # Clone the repository
    git clone https://github.com/cihologramas/pyoptools.git
    cd pyoptools

    # Synchronize environment and install dependencies
    uv sync

    # Build Cython extensions in-place
    uv run python setup.py build_ext --inplace

Incremental Rebuilds
--------------------

When modifying Cython extensions (``.pyx`` or ``.pxd`` files), recompile incrementally with:

.. code-block:: bash

    uv run python setup.py build_ext --inplace

This compiles only the modified modules and typically finishes in 1-2 seconds.

Running Tests
-------------

Run the test suite with pytest:

.. code-block:: bash

    # Run all tests
    uv run pytest

    # Run specific test file
    uv run pytest tests/raytrace/test_system.py -v

Code Quality & Linting
----------------------

All code must conform to the project standards:

.. code-block:: bash

    # Python linting and formatting checks
    uv run ruff check .

    # Cython linting
    uv run cython-lint pyoptools

    # Automated coding standards test
    uv run pytest tests/test_coding_standards.py

Building Documentation Locally
------------------------------

To build the Sphinx documentation locally:

.. code-block:: bash

    uv run --with sphinx --with nbsphinx --with furo --with sphinxcontrib-apidoc sphinx-build -b html doc doc/_build/html

The generated HTML documentation will be located in ``doc/_build/html/index.html``.

Docstrings
----------

All Python and Cython documentation should follow the `NumPy Style Python Docstrings <https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard>`_.

Creating :class:`Surface` Subclasses
------------------------------------

The following rules must be followed when creating a :class:`Surface` subclass:

1. All surface subclasses must be picklable:
   * Inherit from ``pyoptools.misc.picklable.Picklable``.
   * For Cython subclasses, register state attributes using ``self.addkey("key")`` where ``key`` is the attribute name.
   * For pure Python subclasses, no additional registration is required.
