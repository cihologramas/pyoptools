# pyOpTools

pyOpTools is a set of packages that allow the simulation of optical systems by raytracing as well as some calculations involving wavefronts, currently under development. It is written in Python and Cython, and is being developed by the technological development group of Combustión Ingenieros S.A.S, and the applied optics group of the Universidad Nacional de Colombia.

The pyOpTools is divided in several packages that contain the different library's functionalities:

    pyoptools.raytrace

This package contains the classes and functions used to perform the simulation of optical systems using 3D non-sequential raytracing algorithms. 

    pyoptools.misc

This package contains miscellaneous classes and functions used by the other packages, but that can not be classified in any of them. 


## How to install

For Debian/Ubuntu based distro we provide debs packages in the debs folder. You
can also create the deb file by running 'make deb' in the root directory of pyOpTools.
