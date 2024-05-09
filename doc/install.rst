Installing pyOpTools
====================

Getting it
----------

pyOpTools can be downloaded from the project GitHub repository at:

    https://github.com/cihologramas/pyoptools

Installing via pip
------------------

No matter if you're using Windows, Mac OS, or Linux, you can easily install 
pyOpTools. Jupyter Lab is a dependency of the pyoptools package, so it will 
also be installed automatically when you install pyoptools using `pip`. Using 
the installed Jupyter Lab you should be able to  run all the examples in this
documentation without any trouble.

Just run this command:

.. code-block:: bash

    pip install pyoptools

This will install the latest version of pyoptools along with Jupyter Lab, so 
you're all set to run the examples hassle-free.


Installing it in Debian 12 (should work on any Debian derivative)
-----------------------------------------------------------------

The preferred method for installing pyOpTools, particularly on "Debian derivative"
Linux distributions, is to install it as a system package. To generate the 
pyoptools .deb package, simply navigate to the project root and execute the
following command:

.. code-block:: bash

    make deb

Executing this command will generate a .deb package outside the pyOpTool's 
source root tree. To install the generated package, use the following command:

.. code-block:: bash

    dpkg -i python3-pyoptools_<version>_<platform>.deb


Make sure to adjust the <version> and <platform> parameters accordingly. This 
method eliminates the need for virtual environments and ensures a smooth 
installation process.

