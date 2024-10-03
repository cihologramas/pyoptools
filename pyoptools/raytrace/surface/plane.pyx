# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Plane surface definition module
# Symbols Defined: Plane
# ------------------------------------------------------------------------------

'''Module that defines a plane surface class
'''


import cython

from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_nan_to_vector3d

cdef class Plane(Surface):
    '''Class to define a plane surface.

    Description:

    Plane is a class to define rectangular plane surfaces.
    The surface shape is given by the shape attribute

    Example:

        >>> ps=Plane(shape=Rectangular(size=(25,15)))
    '''

    def __init__(self, *args, **kwargs):
        Surface.__init__(self, *args, **kwargs)

    cdef inline double topo_cy(self, double x, double y) noexcept nogil:
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:
        """Returns the intersection point between a ray and an the XY plane

        """
        cdef double u

        # If ray is parallel to the surface, there is no intersection point.
        if incident_ray._direction(2) == 0:
            assign_nan_to_vector3d(intersection_point)
        else:
            # u = dot(N_,-_origin)/dot(N_,_direction)
            u = - incident_ray._origin(2) / incident_ray._direction(2)

            # retval = _origin+u*_direction
            # for the moment u * incident_ray._direction (mult to the left,
            # does not work) it has to be written incident_ray._direction * u
            intersection_point = incident_ray._direction * u+ incident_ray._origin

    cdef void _calculate_normal(self,
                                Vector3d& intersection_point,
                                Vector3d& normal) noexcept nogil:
        """Method that returns the normal to the surface
        """
        (<double*>(&normal(0)))[0] = 0.
        (<double*>(&normal(1)))[0] = 0.
        (<double*>(&normal(2)))[0] = 1.

    def _repr_(self):
        '''Return an string with the representation of the optical plane
        '''

        return "Plane(shape="+str(self.shape)+",reflectivity=" + \
            str(self.reflectivity)+")"
