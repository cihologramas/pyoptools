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
# Description:     Cylindrical surface definitión module
# Symbols Defined: CylindricalSurface
#------------------------------------------------------------------------------
#
'''Module that defines a class to describe cylindrical surfaces
'''

from numpy import power, array, inf, absolute, float64, dot,sqrt as npsqrt
cimport numpy as np

#from enthought.traits.api import Float, Tuple
#from enthought.traits.ui.view import Group, Item

cdef extern from "math.h":
    double sqrt(double)
    double abs(double)

 
#from ray_trace.surface.surface import Surface
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray

cdef class Cylindrical(Surface):
    '''**Class to define cylindrical surfaces.**
    
    Cylindrical is a class to define cylindrical optical surfaces,
    with an aperture given by the shape attribute.

    To define the cylindrical surface you should pass the size of the
    aperture, and the radius of curvature of the cylinder.

    The vertex of the surface is located at the origin of coordinates (0, 0, 0).

    Example:

        >>> cs=Cylindrical(shape=Rectangular(size=(10,20)),curvature=0.15)

    See Surface documentation for other options

    '''

    # Curvature of the cylindrical surface
    #curvature=Float(1.)

    # params is included in the view of Surface to include the specific
    # attributes of this class in the edit_traits window
    
    #params=Group(Item(name="curvature",label="Curvatura de la superficie"))
    cdef public double curvature
    def __init__(self,curvature=0.,*args, **kwargs):
        Surface.__init__(self,*args, **kwargs)
        self.curvature=curvature
        
        self.addkey("curvature")
        #self.shape.topo=self.topo
        
        #~ #Add attributes to the state list        
        #~ self.state.append(self.curvature)
        
    #~ def __reduce__(self):
       #~ 
        #~ args=(self.curvature, self.reflectivity, self.shape)
        #~ return(type(self),args,self.__getstate__())
    
    
    cpdef _intersection(self,Ray A):
        """**Point of intersection between a ray and a surface.**

        This method returns the point of intersection  between the surface
        and the ray. This intersection point is calculated in surface coordinate
        system.
        
        If there is not a point of intersection ie: ray is outside of the
        element aperture, it returns (numpy.inf,numpy.inf,numpy.inf)
        
        Arguments:
        ==========
        
            iray -- incident ray

        iray must be in the coordinate system of the surface

        This function is overloaded from the Surface superclass. It does not 
        check for the aperture.
        
        This function must not be called directly. You should call 
        Surface.intersection instead.
        """
        cdef double x1,y1,z1,x2,y2,z2,z3,x21,y21,z21
        cdef double X1,Y1,Z1,X2,Y2,Z2,X,Y,Z
        cdef double a,b,c,b2ac,u1,u2
        #P1=A.pos     # Punto que pertenece al rayo "Origen" del rayo
        #L1= A.dir    #Vector paralelo a la linea

        #x1,y1,z1=P1
        #x2,y2,z2=P1+L1
        
        

        cdef np.float64_t* pos= <np.float64_t*>(np.PyArray_DATA(A.pos))
        cdef np.float64_t* dir= <np.float64_t*>(np.PyArray_DATA(A.dir))
        
        x1=pos[0];y1=pos[1];z1=pos[2];
        x2=x1+dir[0];y2=y1+dir[1];z2=z1+dir[2];
        x21=dir[0];y21=dir[1];z21=dir[2];
      
        #x3,y3=(0.,0.)
        z3=1./self.curvature

        a= (x2-x1)*(x21)+(z21)*(z21) #power(x2-x1,2)+power(z2-z1,2)
        b=2*((x21)*(x1) + (z21)*(z1 - z3))
        c=x1*x1+z1*z1-2*(z3*z1)#power(x1,2)+power(z1,2)-2*(z3*z1)
        b2ac=b*b-4*a*c
        
        #No hay interseccion
        #o el rayo es tangente y no hay interseccion
        #TODO: Check the first condition
        if a==0:return array([inf,inf,inf])
        if b2ac<=0. : return array([inf,inf,inf])

        u1=(-b+sqrt(b2ac))/(2*a)
        u2=(-b-sqrt(b2ac))/(2*a)

        X1 = x1 + u1 *(x2 - x1)
        Y1 = y1 + u1 *(y2 - y1)
        Z1 = z1 + u1 *(z2 - z1)

        X2 = x1 + u2 *(x2 - x1)
        Y2 = y1 + u2 *(y2 - y1)
        Z2 = z1 + u2 *(z2 - z1)

        if abs(Z2)<abs(Z1):
            X,Y,Z=X2,Y2,Z2
        else:
            X,Y,Z=X1,Y1,Z1
        return array([X,Y,Z],dtype=float64)

    cpdef np.ndarray normal(self,IP):
        """**Normal vector at the point ip.**
        
        This method returns the normal vector at a specific intersection point.
        
        This method is overloaded from the superclass Surface.
        """
        
        N1=IP-array((0.,0.,1./self.curvature)).astype(float64)
        N1[1]=0.
        N_=N1/sqrt(dot(N1,N1))
        return N_

    cpdef topo(self, x, y):
        '''**Method that returns the topography of the surface**
        
        The matrix returned is z=f(x,y). This method is overloaded from the 
        superclass Surface.
        '''
        z=npsqrt((1./self.curvature)**2 -x**2 )-absolute(1./self.curvature)
        if self.curvature>0: z=-z
        return z

    def _repr_(self):
        '''Return an string with the representation of the optical cylindrical surface
        '''

        return "Cylindrical(shape="+str(self.shape)+",reflectivity="+\
                          str(self.reflectivity)+",curvture="+\
                          str(self.curvature)+")"

