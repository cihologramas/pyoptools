Installing pyOpTools
====================

Requirements
------------

* Python 3.10 or newer (tested on Python 3.10 - 3.13)
* Linux, macOS, or Windows
* A C/C++ compiler and Eigen headers (if installing from source or building extensions)

Installing with uv (Recommended)
--------------------------------

If you use `uv <https://docs.astral.sh/uv/>`_:

.. code-block:: bash

    # In a uv project:
    uv add pyoptools

    # Or in a virtual environment:
    uv pip install pyoptools

Installing via pip
------------------

You can install the latest release of pyOpTools directly from PyPI using ``pip``:

.. code-block:: bash

    pip install pyoptools

For interactive notebooks with 3D visualization, you can also install JupyterLab and Plotly:

.. code-block:: bash

    pip install pyoptools jupyterlab plotly

Building from Source
--------------------

To install pyOpTools from source for development:

.. code-block:: bash

    # 1. Clone repository
    git clone https://github.com/cihologramas/pyoptools.git
    cd pyoptools

    # 2. Install dependencies and compile extensions
    uv sync
    uv run python setup.py build_ext --inplace

    # Or with pip:
    pip install -e .

Installing in Debian / Ubuntu as a System Package
-------------------------------------------------

To generate and install a Debian ``.deb`` package:

.. code-block:: bash

    make deb
    sudo dpkg -i ../python3-pyoptools_<version>_<platform>.deb
