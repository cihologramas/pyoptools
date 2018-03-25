#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# cython: profile=True

#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amezquita Orozco
# Description:     Rectangle definition module
# Symbols Defined: Rectangle
#------------------------------------------------------------------------------
#

"""Module that defines the Rectangular class
"""

#from enthought.traits.api import HasTraits, Float, Tuple, Int
#from enthought.traits.ui.view import View, Item, Include,Group

from pyoptools.raytrace.shape.shape cimport Shape


from numpy import arange, meshgrid, where,linspace

cdef class Rectangular(Shape):
    ''' Class defining an rectangular shape. '''
    
    # Tuple that holds the size of the rectangle
    #size=Tuple(Float(1.),Float(1.))
    
    # Tuple that holds the number of samples to be used to build the mesh
    #samples=Tuple(Int(30), Int(30))
    
    def __init__(self,size=(1.,1.),samples=(30,30), offset=(0,0),*args, **kwargs):
        Shape.__init__(self,*args, **kwargs)
        self.size=(float(size[0]),float(size[1]))
        self.samples=(float(samples[0]),float(samples[1]))
        self.offset=(float(offset[0]),float(offset[1]))
        self.addkey("size")
        self.addkey("samples")
        self.addkey("offset")
    
    cpdef hit(self, p):
        """Method  that returns True if a p=(x,y,z) point is inside the rectangle,
        if not it returns False.
        """
        x, y, z=p
        dx,dy=self.size
        ox,oy=self.offset
        x=x-ox
        y=y-oy
        retval= where((x<-dx/2.) | (x>dx/2.) | (y<-dy/2.) | (y>dy/2.), False, True)
        return retval
        
    cpdef bint fhit(self,double px,double py,double pz):
        """This method returns TRUE if an p=(x,y,z)point is inside the surface 
        aperture if not it must return FALSE.
        This is implemented for a point, in cython, to make it fast
        """
        cdef double dx,dy,ox,oy,opx,opy
        dx,dy=self.size
        ox,oy=self.offset
        opx=px-ox
        opy=py-oy
        if (opx<-dx/2.) | (opx>dx/2.) | (opy<-dy/2.) | (opy>dy/2.):
            return False
        else: return True
        
        
    #~ cpdef polylist(self, topo):
        #~ """Method that returns a tuple (point_list, poly_list) for a rectangular mesh. 
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
        #~ dx,dy=self.size
        #~ nx, ny=self.samples
        #~ Xl=linspace(-dx/2.,dx/2.,nx)
        #~ Yl=linspace(-dy/2.,dy/2.,ny)
        #~ X, Y=meshgrid(Xl,Yl)
        #~ Z=topo(X, Y)
        #~ 
        #~ xs,ys= Z.shape
#~ 
        #~ points=[]
        #~ polylist=[]
        #~ i=0
        #~ for x in range(0,xs):
            #~ for y in range(0,ys):
                #~ points.append([X[x,y],Y[x,y],Z[x,y]]);
                #~ if y>0 and x>0:
                    #~ polylist.append([i-1,i,i-ys,i-ys-1])
                #~ i=i+1
        #~ 
        #~ return points, polylist
    
    cpdef pointlist(self):
        dx,dy=self.size
        nx, ny=self.samples
        ox, oy=self.offset
        Xl=linspace(-dx/2.+ox,dx/2.+ox,nx)
        Yl=linspace(-dy/2.+oy,dy/2.+oy,ny)
        X, Y=meshgrid(Xl,Yl)
        return X.ravel(),Y.ravel()
        

    cpdef limits(self):
        """
        Returns the minimum limits for the aperture
        """
        dx, dy=self.size
        ox, oy=self.offset
        return -dx/2+ox, dx/2+ox, -dy/2+oy, dy/2+oy
