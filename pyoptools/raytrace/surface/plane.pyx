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
from numpy import array, float64, zeros, asarray

# from ray_trace.surface import Surface
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.definitions import *

from pyoptools.misc.cmisc.linalg cimport Vector3, Matrix3x3, \
     vector3_from_python_object, vector3_to_tuple, normalize_vector3, \
     matrix3x3_vector3_dot, compute_rotation_matrix, compute_rotation_matrix_i, \
     add_vector3, substract_vector3, vector3_magnitude, vector3_equals, vector3_dot_product, \
     set_nan_vector3, vector3_times_scalar, negate_vector3_inplace

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

    cpdef topo(self, x, y):
        return zeros(asarray(x).shape)

    # ~ def __reduce__(self):
        # ~
        # ~ args=(self.reflectivity, self.shape)
        # ~ return(type(self),args,self.__getstate__())

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void _calculate_intersection(self, Ray incident_ray, Vector3 *intersection_point_ptr) noexcept nogil:
        """Returns the intersection point between a ray and an the XY plane

        """
        cdef double u
        cdef Vector3 tmp_vector

        # If ray is parallel to the surface, there is no intersection point.
        if incident_ray._direction.data[2] == 0:
            set_nan_vector3(intersection_point_ptr)
        else:
            # u=dot(N_,-_origin)/dot(N_,_direction) =-_origin[2]/_direction[2]
            u = - incident_ray._origin.data[2] / incident_ray._direction.data[2]

            #  retval = _origin+u*_direction
            
            vector3_times_scalar(&(incident_ray._direction), u, &tmp_vector)
            add_vector3(&(incident_ray._origin),&tmp_vector,intersection_point_ptr)

    cdef void _calculate_normal(self, Vector3 *intersection_point_ptr, Vector3 *normal_ptr) noexcept nogil:
        """Method that returns the normal to the surface
        """
        normal_ptr[0].data[0]= 0.
        normal_ptr[0].data[1]= 0.
        normal_ptr[0].data[2]= 1.

    def _repr_(self):
        '''Return an string with the representation of the optical plane
        '''

        return "Plane(shape="+str(self.shape)+",reflectivity=" + \
            str(self.reflectivity)+")"
