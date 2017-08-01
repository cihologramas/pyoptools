# pyOpTools

pyOpTools is a set of packages that allow the simulation of optical systems by raytracing as well as some calculations involving wavefronts, currently under development. It is written in Python and Cython, and is being developed by the technological development group of Combustión Ingenieros S.A.S, and the applied optics group of the Universidad Nacional de Colombia.

The pyOpTools is divided in several packages that contain the different library's functionalities:

    pyoptools.raytrace

This package contains the classes and functions used to perform the simulation of optical systems using 3D non-sequential raytracing algorithms. 

    pyoptools.misc

This package contains miscellaneous classes and functions used by the other packages, but that can not be classified in any of them. 

## Requirements

The following steps work to install the packages required to run pyoptools in an ipython notebook under debian 8. These
instructions are far from complete, but they will give an idea.

1. Install stdeb. in 
https://pypi.python.org/pypi/stdeb/0.8.5#install-or-using-stdeb-to-create-an-stdeb-installer you can find instructions 
on how to do it.
2. Install ipython: sudo pypi-install ipython --release 3.2.1 (had problems with the latest one)
3. sudo apt-get install python-setuptools
4. sudo pypi-install backports.ssl_match_hostname
5. sudo pypi-install certifi
6. sudo pypi-install tornado
7. sudo apt-get install python-jsonschema
8. sudo apt-get install libosmesa6
9. sudo apt-get install python-wxgtk3.0
10. sudo pypi-install PyOpenGL --release 3.1.0

**NOTE:** These instructions will only work for debian using stdeb, adapt accordingly for your platform/OS.

##How to install

Pyoptools can be installing using one of the following ways

### Debian/Ubuntu (recommended)
We develop pyoptools using almost debian/*buntu Linux distributions exclusively so this is the way that's most tested 
and have better chances of working. To compile a deb package you just need to run:

    make deb

This will create a .deb installer file inside the `deb_dist` directory in your source root tree of pyoptools. Which you
can install using the command (with admin privileges), changing the version/release/platform accordingly:

    dpkg -i python-pyoptools_<version>-<release>_<platform>.deb

Using this method you will have control of versions installed via dpkg/APT package manager.

### PyPI install (using pip)

Pyoptools is hosted in the Python Package Index since Jul-17, therefore you can just use `pip` to install it via command
line with the following command

    pip install pyoptools

### Virtual machine (VBox)

We have created a virtualbox image of a debian 8 installation with the pyoptools running. This image can be downloaded from https://drive.google.com/open?id=0B6vN2VIpMQ48THUwX1NJSjM0cEE

* user name: usuario
* password: usuario

This is to help people interested in testing this tool, while the documentation and installation instructions get updated.


For Debian/Ubuntu based distro we provide debs packages in the debs folder. You
can also create the deb file by running 'make deb' in the root directory of pyOpTools.

For development `sudo python setup.py develop` this will install systemwide from current directory, making changes to the *py files directly availables.