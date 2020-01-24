General Description
===================

pyOpTools is a free and open source `python <http://python.org>`_ library to the easily simulate optical systems. It is composed by the following packages:

:raytrace:
    Package that contain classes and functions used to simulate optical systems by 3D non sequential raytracing. The classes in the library, allow you to define optical components as a set of bounding surfaces that enclose some material (typically some kind of optical glass). With such components (using its positions and orientations) it is then possible to define an optical system. When the system is ready, it is possible to include light rays, and then to propagate them to obtain the result of the simulation of the system.  

    This is the most developed package of the library.

:wavefront:
   Package containing classes and routines used to simulate the propagation of wavefronts through optical systems. This package is very incomplete.

:misc:
   Package that contains auxiliary functions (not related to optics), that are used in the raytrace and wavefront package.

:gui:
   Package that contains classes and functions that allow the pyoptools simulations to be displayed in the `jupyter <https://jupyter.org/>`_ notebooks. 
