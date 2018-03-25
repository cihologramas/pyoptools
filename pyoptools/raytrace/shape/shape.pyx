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
# Description:     Shape definition module
# Symbols Defined: Shape
#------------------------------------------------------------------------------
#

"""Module that defines the Shape superclass
"""

from warnings import warn

#from enthought.traits.api import HasTraits, Float, Trait, List, Function,  Method
#from enthought.traits.ui.view import View, Item, Include,Group

from numpy import arange,  meshgrid, linspace
from pyoptools.misc.picklable.picklable cimport Picklable

cdef class Shape(Picklable):
    ''' Abstract superclass for all optical surfaces shapes


    Description:

    Shape is an abstract superclass to define the different surfaces shapes 
    (circular, rectangulas, etc.).
    This class defines an API that all subclasses must support.

    The methods that must be overloaded in the subclasses are:

        Shape.hit
            This method must return TRUE if an X,Y point is inside the surface 
            aperture if not it must return FALSE
        
        Shape.polylist
            This method should return a point_list, and a poly_list. The point 
            list is a list of tuples (X,Y,Z) containing the coordinates of the 
            points used to build the surface mesh.
            The poly_list is a list of tuples (n1,n2,n3,n3) containing the 
            indices of the points in the polylist used to build each polygon 
            that will be used to visualize the mesh.
            
    This methods will be called by the Shape.hit method, and to create the 
    Shape.polylist tuple, that will a point_list and a poly_list.
    The point_list is a list of tuples (X,Y,Z) containing the coordinates of
    the points used to build the surface mesh.
    The poly_list is a list of tuples (n1,n2,n3,n3) containing the indices 
    of the points in the polylist used to build each polygon that will be 
    used to visualize the mesh. 
    
    The attribute topo, will be initialized for each instance of the class with 
    a function Z(x,y).
    '''

    def __init__(self):
        #self.topo=None
        Picklable.__init__(self)
                
    
    cpdef hit(self, p):
        """This method must return TRUE if an p=(x,y,z)point is inside the surface 
        aperture if not it must return FALSE. This must work for a list of points.
        
        It must be overloaded.
        """
        warn("Method hit, from class Shape, should be overloaded"+
             " in class "+self.__class__.__name__)
    
    cpdef bint fhit(self,double px,double py,double pz):
        """This method returns TRUE if an p=(x,y,z)point is inside the surface 
        aperture if not it must return FALSE.
        This is implemented for a point, in cython, to make it fast.
        
        It must be overloaded.
        """
        warn("Method hitf, from class Shape, should be overloaded"+
             " in class "+self.__class__.__name__)
             
    
    #~ cpdef polylist(self, topo):
        #~ """This method should return a tuple (point_list, poly_list). 
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
        #~ 
        #~ It must be overloaded, for each shape
        #~ """
        #~ 
        #~ warn("Method polylist, from class Surface, should be overloaded"+
             #~ " in class "+self.__class__.__name__)
             
    cpdef pointlist(self):
        """This method should return a point_list that samples adecuately the shape. 
        
        The return value must be a tuple (X,Y) where X contains the X coordinates 
        of the points, and Y the y coordinates. 
        It must be overloaded, for each shape
        """
        
        warn("Method pointlist, from class Surface, should be overloaded"+
             " in class "+self.__class__.__name__)
        

    def mesh(self, size=None, ndat=(100, 100), topo=None):
        """
        The X,Y are a mesh limited by size=(xi,xf,yi,yf)  
        if size is not given (size=None), the size is given by the aperture.
        topo is the topography function of the surface under study.
        if not given, a 0 will be returned for points outside the aperture, and
        1 inside the aperture
        
        The ndat tuple, gives the number of points in the mesh.
        """
        if size==None:
            xi, xf, yi, yf=self.limits()
        else:
            xi,xf, yi, yf =size
            
        nx, ny = ndat
        X_M=linspace(xi, xf, nx)
        Y_M=linspace(yi, yf, ny)
        #X_M=arange(xi,xf*(1.+1/nx),1./nx)
        #Y_M=arange(yi,yf*(1.+1/ny),1./ny)

        XM, YM=meshgrid(X_M,Y_M)
        if topo==None:
            Z=self.hit((XM, YM, 0))
        else:
            Z=topo(XM, YM)
        return XM, YM, Z
        
    cpdef limits(self):
        """
        Returns the minimum limits for the aperture
        It must be overloaded, for each shape
        """


