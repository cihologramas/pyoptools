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
# Description:     Spherical surface definitión module
# Symbols Defined: Spherical
#------------------------------------------------------------------------------

cdef extern from "math.h":
    double sqrt(double)
    double abs(double)
    

    
'''Module that defines an optical spherical surface
'''

from numpy import array, inf, sqrt as npsqrt, absolute, float64, dot, arcsin, pi,  zeros
cimport numpy as np

from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray


cdef class Spherical(Surface):
    '''**Class to define spherical surfaces.**
    
    Spherical is a class to define spherical optical surfaces.

    To define the spherical surface you should pass the shape
    of the aperture, and the radius of curvature of the sphere.

    The vertex of the surface is located at the origin of coordinates (0, 0, 0)
    and the aperture is centered also about the origin

    Example:

        >>> cs=Spherical(shape=Circular(radius=60),curvature=0.15)

    See Surface documentation for other options
    '''

    # Curvature of the spherical surface
#    curvature=Float(0.)
    
    # params is included in the view of Surface to include the specific
    # attributes of this class in the edit_traits window

#    params=Group(Item(name="radius",label="Radio de la apertura"),
#                 Item(name="curvature",label="Curvatura de la Superficie"))
#   
    cdef public double curvature
    def __init__(self,curvature=0., *args, **kwargs):
        Surface.__init__(self,*args, **kwargs)
        self.curvature=curvature
        self.addkey("curvature")
    #~ def __reduce__(self):
        #~ 
        #~ args=(self.curvature, self.reflectivity, self.shape)
        #~ return(type(self),args,self.__getstate__())
    #~
        
    cpdef _intersection(self,Ray A):
        '''**Point of intersection between a ray and the sphere**

        This method returns the point of intersection  between the surface
        and the ray. This intersection point is calculated in the coordinate
        system of the surface.
        
           iray -- incident ray

        iray must be in the coordinate system of the surface
        '''
            
        #P1=A.pos     # Punto que pertenece al rayo "Origen" del rayo
        #L1= A.dir    #Vector paralelo a la linea
        cdef double x1,y1,z1,x2,y2,z2,x21,y21,z21,z3,a,b,c,b2ac,u1,u2,X1,Y1,Z1,X2,Y2,Z2
        
        cdef np.float64_t* pos= <np.float64_t*>(np.PyArray_DATA(A.pos))
        cdef np.float64_t* dir= <np.float64_t*>(np.PyArray_DATA(A.dir))
        
        #x1,y1,z1=A.pos
        #x2,y2,z2=A.pos+A.dir #P1+L1
        
        x1=pos[0];y1=pos[1];z1=pos[2];
        x2=x1+dir[0];y2=y1+dir[1];z2=z1+dir[2];
        
        
        #x21=x2-x1; y21=y2-y1; z21=z2-z1
        
        
        #x21, y21, z21= A.dir #L1
        x21=dir[0]; y21=dir[1]; z21=dir[2];
        
        #x3,y3=(0.,0.)
        z3=1./self.curvature
        
        a=x21**2+y21**2+z21**2
        b=2*((x21)*(x1) + (y21)*(y1) + (z21)*(z1 - z3))
        c=x1**2+y1**2+z1**2-2*(z3*z1)
        b2ac=b**2-4*a*c
        
        #print a, b, c, b2ac
        
        if b2ac<=0. : return array([inf,inf,inf])

        
        u1=(-b+sqrt(b2ac))/(2*a)
        u2=(-b-sqrt(b2ac))/(2*a)

        X1 = x1 + u1 *(x21)
        Y1 = y1 + u1 *(y21)
        Z1 = z1 + u1 *(z21)

        X2 = x1 + u2 *(x21)
        Y2 = y1 + u2 *(y21)
        Z2 = z1 + u2 *(z21)

        ##No hay interseccion
        ##o el rayo es tangente y no hay interseccion 
            
        # TODO: This can have problems if the ray is propagates in the X or Y direcction 
        # Need to find a beter solution 

        if abs(Z2)<abs(Z1):
            X,Y,Z=X2,Y2,Z2
        else:
            X,Y,Z=X1,Y1,Z1

        # This spherical surface is only half surface, so:

        if abs(Z)>abs(z3):
            return array([inf,inf,inf])


        return array([X,Y,Z]).astype(float64)

    cpdef np.ndarray normal(self,ip):
        """**Return the vector normal to the surface**
        
        This method returns the vector normal to the sphere at a point 
        ``int_p=(x,y,z)``.
        """

        N1=ip-array((0.,0.,1./self.curvature)).astype(float64)
        N_=N1/sqrt(dot(N1,N1))
        return N_

    cpdef topo(self, x, y):
        """**Returns the Z value for a given X and Y**
        
        This method returns the topography of the spherical surface to be 
        used to plot the surface.
        """
        z=npsqrt((1./self.curvature)**2 -x**2 -y**2)-absolute(1./self.curvature)
        if self.curvature>0: z=-z
        return z
        
    def _repr_(self):
        '''Return an string with the representation of the optical spherical surface
        '''

        return "Spherical(shape="+str(self.shape)+",reflectivity="+\
                          str(self.reflectivity)+",curvture="+\
                          str(self.curvature)+")"

