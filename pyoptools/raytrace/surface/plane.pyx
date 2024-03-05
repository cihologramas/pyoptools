#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# cython: profile=True

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
from numpy import array, dot, inf, float64, zeros, asarray
#from enthought.traits.api import Tuple,Float
#from enthought.traits.ui.view import Group,Item

#from ray_trace.surface import Surface
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.definitions import *

cimport numpy as np
np.import_array()


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
 # ~

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _intersection(self, Ray A):
        """Returns the intersection point between a ray and an the XY plane

        """
        # N_=array([0.,0.,1.])

        # Punto que pertenece al rayo "Origen" del rayo
        cdef np.ndarray[np.float64_t, ndim=1] P1 = A.pos
        # Vector paralelo a la linea
        cdef np.ndarray[np.float64_t, ndim=1] L1 = A.dir

        # if dot(N_,L1) ==0 : return inf_vect
        if L1[2] == 0:
            return inf_vect

        # print N_,P1,L1
        # print dot(N_,-P1),dot(N_,L1)
        # u=dot(N_,-P1)/dot(N_,L1)
        cdef double u = -P1[2]/L1[2]
        # Si u es muy grande, no hay intersección

        retval = P1+u*L1
        #from sys import exit
        # if isnan(retval[0]):
        #    print P1,u,L1
        #    print A.dir
        #    print type(A.orig_surf)
        #    exit(0)

        return retval

    cpdef np.ndarray normal(self, ri):
        """Method that returns the normal to the surface
        """
        N_ = array((0., 0., 1.)).astype(float64)
        return (N_)

    def _repr_(self):
        '''Return an string with the representation of the optical plane
        '''

        return "Plane(shape="+str(self.shape)+",reflectivity=" + \
            str(self.reflectivity)+")"
