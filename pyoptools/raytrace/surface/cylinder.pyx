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
# Description:     Tube or hollow cylinder surface definition module
# Symbols Defined: Cylinder
#------------------------------------------------------------------------------
#
'''Module that defines a class to describe cylindrical surfaces
'''



from numpy import pi, power, array, inf, sqrt, absolute, float64, dot, sin, cos
cimport numpy as np
#from enthought.traits.api import Float, Tuple
#from enthought.traits.ui.view import Group, Item
#from enthought.tvtk.api import tvtk

#from ray_trace.surface.surface import Surface
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray

cdef class Cylinder(Surface):
    '''**Class to define cylinder shaped surfaces.**

    Description:

    Cylinder is a class to define a tube or  a hollow cylinder surface.

    To define the cylinder surface you should pass the radius of the cylinder, 
    and its length

    The center of the cylinder is located at the origin of coordinates (0, 0, 0)
    and its length is parallel to the Z axis.

    Example:

        >>> cs=Cylinder(radius=5.,length=10.)

    See Surface documentation for other options

    '''

    # Radius of the cylinder
    #radius=Float(1.)
    
    # Length of the cylinder
    #length=Float(1.)

    # params is included in the view of Surface to include the specific
    # attributes of this class in the edit_traits window
    
    #params=Group(Item(name="radius",label="Radio del cllindro"),
    #             Item(name="length",label="Longitud del cilindro"))
    cdef public double radius, length
    def __init__(self,radius=10., length=10.,*args, **kwargs):
        Surface.__init__(self,*args, **kwargs)
        self.radius=radius
        self.length=length   

        #Add attributes to the key list        
        self.addkey("radius")
        self.addkey("length")
        
    #~ def __reduce__(self):
     #~ 
        #~ args=(self.radius, self.length, self.reflectivity, self.shape)
        #~ return(type(self),args,self.__getstate__())

        

    cpdef np.ndarray intersection(self, Ray A):
        '''**Point of intersection between a ray and the cylinder**

        This method returns the point of intersection  between the surface
        and the ray. This intersection point is calculated in the coordinate
        system of the surface.
        
        
        
           iray -- incident ray

        iray must be in the coordinate system of the surface
        
        Note: Because of the way the cylinder is defined, it does not use
        shapes to define its boundary, for that reason, the ``intersection``
        method and not the ``_intersection method`` was overloaded.
        '''
        
        
        P1=A.pos     # Punto que pertenece al rayo "Origen" del rayo
        L1= A.dir    #Vector paralelo a la linea

        x1,y1,z1=P1
        dx,dy,dz=L1

        # This is done by obtaining L from 2 sets of equation
        # x^2+y^2-R^2=0
        # and
        # x=x1+L*dx
        # y=y1+L*dx
        a=power(dx,2)+power(dy,2)
        b=2*(dx*x1 + dy*y1)
        c=power(x1,2)+power(y1,2)-power(self.radius,2)
        b2ac=power(b,2)-4*a*c
        #No hay interseccion
        #o el rayo es tangente y no hay interseccion
        if b2ac<=0. : return array([inf,inf,inf])

        u1=(-b+sqrt(b2ac))/(2*a)
        u2=(-b-sqrt(b2ac))/(2*a)
        
        # if u1 or u2 = 0, the ray already intersected the cylinder, so u1 and
        # u2 must be >0 
        if u1<1e-10:u1=inf
        if u2<1e-10:u2=inf
        
        X1 = x1 + u1 * dx
        Y1 = y1 + u1 * dy
        Z1 = z1 + u1 * dz

        X2 = x1 + u2 * dx
        Y2 = y1 + u2 * dy
        Z2 = z1 + u2 * dz

        # Check for the nearest intersection
        d1=power(X1-x1,2)+power(Y1-y1,2)+power(Z1-z1,2)
        d2=power(X2-x1,2)+power(Y2-y1,2)+power(Z2-z1,2)

        if d2<d1 :
            X,Y,Z=X2,Y2,Z2
        else:
            X,Y,Z=X1,Y1,Z1
        
        if Z<-self.length/2. or Z>self.length/2.:
            return array([inf,inf,inf])
        else: return array([X,Y,Z]).astype(float64)

    cpdef polylist(self):
        """
        Because this is a closed surface, the method had to be overloaded
        """    
        Z1=-self.length/2.
        Z2=self.length/2.
        points=[]
        polys=[]
        for theta in range(40):
            th=theta*2.*pi/40.
            X=self.radius*cos(th)
            Y=self.radius*sin(th)
            points.append((X,Y,Z1))
            points.append((X,Y,Z2))
            
        for i in range(80-2):
            polys.append((i,i+1,i+2))
        
        polys.append((78,79,0))
        polys.append((79,0,1))
        return points,polys
        
    cpdef np.ndarray normal(self,  int_p):
        """**Normal vector at the point int_p.**
        
        This method returns the normal vector at a specific intersection point, 
        given by int_p.
        """
         
        N1=int_p.copy() # Need to do a copy so the intersection point does not get modified
        N1[2]=0 # the normal should not have Z component
        N_=N1/sqrt(dot(N1,N1))
        return N_
 

    def _repr_(self):
        '''Return an string with the representation of the cylinder
        '''
        return "Cylindrical(shape="+str(self.shape)+",reflectivity="+\
                          str(self.reflectivity)+",radius="+\
                          str(self.radius)+",length="+str(self.length)+")"

