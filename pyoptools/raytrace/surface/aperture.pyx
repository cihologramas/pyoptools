#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# cython: profile=True

#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Stops definitión module
# Symbols Defined: Stop, StopCircAp, StopRecAp
#------------------------------------------------------------------------------

'''Module that defines stops classes
'''



from numpy import power, arange, pi, zeros, cos, sin, tan
#from enthought.traits.api import Float, Tuple,  Instance

#from ray_trace.surface import Plane
from pyoptools.raytrace.surface.plane cimport Plane

from pyoptools.raytrace.shape.shape  cimport Shape
from pyoptools.raytrace.ray.ray cimport Ray


cdef class Aperture(Plane):
    '''
    Class to define a surface with an aperture
    
    This class is used to define stops in the optical system. It receives two 
    parameters.
    
    **ARGUMENTS:**
    
    ======== ================================================================
    shape    It is a subclass of Shape, and defines the external shape of the 
             stop
    ap_shape It is a subclass of Shape, and defines the internal shape of the 
             stop (the aperture shape).
    ======== ================================================================
    
    **EXAMPLE**
    
    Creation of an aperture surface::
    
        ap=Aperture(shape=Rectangular(size=(60,60)),ap_shape=Circular(radius=2.5))
        
    Note: To create an stop component, use the Stop class from the comp_lib 
    module. It creates the aperture surface and encapsulate it in an component 
    that can be used in a System.
    '''
    
    # Attribute that defines the aperture shape (the shape of the hole)
    #ap_shape=Instance(Shape)
    cdef public Shape ap_shape
    
    def __init__(self,ap_shape=None,*args, **kwargs):
        Plane.__init__(self,*args, **kwargs)
        self.ap_shape=ap_shape
        
        #Add items to the state list
        self.addkey("ap_shape")

        
    #~ def __reduce__(self):
        #~ 
        #~ args=(self.ap_shape, self.reflectivity, self.shape)
        #~ return(type(self),args,self.__getstate__())
    #~ #def __setstate__(self, state):        
    #~ #    for child in state:
    #~ #        self.add_child(child)
    #~ 

    cpdef propagate(self,Ray ri,double ni,double nr):
        """
        The OptSurf.propagate is overloaded so it can be decided if the rays
        continue propagating or not.
        
        Warning: This surface only checks if the ray continues or not. It does 
        not calculate refraction or reflection. It must not be used to create
        lenses or mirrors.
        """
        
        # Calculate the intersection point and the surface normal
        PI,P=self.int_nor(ri)
        X,Y,Z=PI
        
        if self.ap_shape.hit(PI)==True:
            i=ri.intensity
        else:
            i=0.
        ret_ray=Ray(pos=PI,dir=ri.dir,intensity=i,wavelength=ri.wavelength,
                n=ri.n, label= ri.label, orig_surf=self.id)
        return [ret_ray]
