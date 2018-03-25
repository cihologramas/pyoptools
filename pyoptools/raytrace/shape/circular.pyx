#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# cython: profile=True

#------------------------------------------------------------------------------
# Copyright (c) 2007,  2008, 2009 Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amezquita Orozco
# Description:     Rectangle definition module
# Symbols Defined: Circular
#------------------------------------------------------------------------------
#

""" Module that defines the Circular superclass """


#from enthought.traits.api import HasTraits, Float, Tuple, Int
#from enthought.traits.ui.view import View, Item, Include,Group
from numpy import power, sin, cos, pi, where 
from pyoptools.raytrace.shape.shape cimport Shape


from numpy import arange, meshgrid

cdef class Circular(Shape):
    ''' Class defining an circular shape. '''
    
    # Radius of the circular shape
    #radius=Float(1.)
    
    # Tuple that holds the number of samples to be used to build the mesh. 
    # The first number gives radial samples, the second angular samples
    #samples=Tuple(Int(10), Int(36))
    
    def __init__(self,radius=1.,samples=(10,36),*args, **kwargs):
        Shape.__init__(self,*args, **kwargs)
        self.radius=radius
        self.samples=samples    
        self.addkey("radius")
        self.addkey("samples")
        
    def __reduce__(self):
        
        state=None
        args=(self.radius, self.samples)
        return(type(self),args)
    
    cpdef hit(self, p):
        """Method that returns True if a p=(X,Y,Z) point is inside the aperture,
        if not it returns False.
        """
        x, y, z=p
        return where (power(x,2)+power(y,2)> power(self.radius,2), False, True)
        
    cpdef bint fhit(self,double px,double py,double pz):
        """This method returns TRUE if an p=(x,y,z)point is inside the surface 
        aperture if not it must return FALSE.
        This is implemented for a point, in cython, to make it fast
        """
             
        if px**2+py**2 > self.radius**2 :
            return False
        else: return True
         
    #~ cpdef polylist(self, topo):
        #~ """Method that returns a tuple (point_list, poly_list) for a circular mesh. 
        #~ 
        #~ Attributes:
        #~ ===========
        #~ 
        #~ topo    Z=topo(x,y) is the function that gives the surface topography
        #~ 
        #~ The point list is a list of tuples (X,Y,Z) containing the coordinates of
        #~ the points used to build the surface mesh.
        #~ The poly_list is a list of tuples (n1,n2,n3,n3) containing the indices 
        #~ of the points in the polylist used to build each polygon that will be 
        #~ used to visualize the mesh.
        #~ """
        #~ 
        #~ nr, na=self.samples
        #~ 
        #~ R=arange(self.radius*(1./nr),self.radius*(1.+1./(2*nr)),self.radius/nr)
        #~ A=arange(0,2*pi,2*pi/na)
#~ 
        #~ # Create the point list, and locate the vertex as the first point
        #~ 
        #~ RM, AM=meshgrid(R,A)
        #~ #print RM
        #~ #print AM
#~ 
        #~ X=RM*cos(AM)
        #~ Y=RM*sin(AM)
        #~ 
        #~ #print X
        #~ #print Y
        #~ 
        #~ Z=topo(X, Y)
        #~ # Create the list leaving the surface open at the vertex, and at a=0
        #~ xs,ys= Z.shape
#~ 
        #~ points=[(0., 0., topo(0., 0.)), ]
        #~ polylist=[]
        #~ i=1
        #~ for x in range(0,xs):
            #~ for y in range(0,ys):
                #~ points.append([X[x,y],Y[x,y],Z[x,y]]);
                #~ if y>=1 and x>=1 and Z[x, y]:
                    #~ polylist.append([i-1,i,i-ys,i-ys-1])
                #~ i=i+1
        #~ # Close the vertex
        #~ 
        #~ for x in range(0,xs-1):
            #~ polylist.append([0,x*ys+1,(x+1)*ys+1])
        #~ polylist.append([0, 1, ys*(xs-1)+1])
        #~ 
        #~ # Close the surface
        #~ for y in range(1, ys):
            #~ polylist.append([y,y+1, ys*(xs-1)+y+1, ys*(xs-1)+y ])
        #~ return points, polylist
        
    
    cpdef pointlist(self):
        
        nr, na=self.samples
        
        R=arange(self.radius*(1./nr),self.radius*(1.+1./(2*nr)),self.radius/nr)
        A=arange(0,2*pi,2*pi/na)
        
        RM, AM=meshgrid(R,A)
        #print RM
        #print AM

        X=RM*cos(AM)
        Y=RM*sin(AM)
        
        return X.ravel(),Y.ravel()
        

    cpdef limits(self):
        """
        Returns the minimum limits for the aperture
        """
        return -self.radius,self.radius,-self.radius,self.radius

