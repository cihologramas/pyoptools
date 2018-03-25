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
# Symbols Defined: Triangle
#------------------------------------------------------------------------------
#

"""Module that defines the Triangular class
"""


from pyoptools.raytrace.shape.shape cimport Shape


from numpy import arange, meshgrid, where, dot, array

cdef class Polygon(Shape):
    ''' Class defining a triangular shape. 

    coord -> tuple containing the n (x,y) coordinates of the corners of
             polygone
    samples -> number of divitions per side used to sample the triangle
    '''


    
    def __init__(self,coord=((0,0),(0,100),(100,0)),samples=10,*args, **kwargs):
        Shape.__init__(self,*args, **kwargs)
        self.coord=coord
        self.samples=samples
        
        #Register picklable attributes
        self.addkey("coord")
        self.addkey("samples")
        
    def __reduce__(self):
       
        state=None
        args=(self.coord, self.samples)
        return(type(self),args)
    
        
    cpdef hit(self, p):
        """Method  that returns True if a p=(x,y,z) point is inside the triangle,
        if not it returns False.
        taken from http://www.blackpawn.com/texts/pointinpoly/default.html
        """
        x, y, z=p
        P=array((x,y))
        A=array(self.coord[0])
        B=array(self.coord[1])
        C=array(self.coord[2])
        
        v0=C-A
        v1=B-A
        v2=P-A
        
        dot00=dot(v0,v0)
        dot01=dot(v0,v1)
        dot02=dot(v0,v2)
        dot11=dot(v1,v1)
        dot12=dot(v1,v2)
        
        invDenom=1./(dot00 * dot11 - dot01 * dot01)

        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom

        #Check if point is in triangle
        return (u > 0) and (v > 0) and (u + v < 1)

        
    cpdef bint fhit(self,double px,double py,double pz):
        """This method returns TRUE if an p=(x,y,z)point is inside the surface 
        aperture if not it must return FALSE.
        This is implemented for a point, in cython, to make it fast
        """
        cdef double dot00,dot01,dot02,dot11,dot12,invDenom,u,v
        # This one needs to be optimized
        P=array((px,py))
        A=array(self.coord[0])
        B=array(self.coord[1])
        C=array(self.coord[2])
        
        v0=C-A
        v1=B-A
        v2=P-A
        
        
        dot00=dot(v0,v0)
        dot01=dot(v0,v1)
        dot02=dot(v0,v2)
        dot11=dot(v1,v1)
        dot12=dot(v1,v2)
        
        invDenom=1./(dot00 * dot11 - dot01 * dot01)

        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom

        # Check if point is in triangle
        return (u > 0) and (v > 0) and (u + v < 1)

        
    #~ cpdef polylist(self, topo): #Falta organizar el polilist
        #~ """Method that returns a tuple (point_list, poly_list) for a triangular mesh. 
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
        #~ cdef int i,j
        #~ 
        #~ A=array(self.coord[0])
        #~ B=array(self.coord[1])
        #~ C=array(self.coord[2])
        #~ 
        #~ #Get the mesh points
        #~ points=[]
        #~ for i in range(self.samples+1):
            #~ P0= A+i*(B-A)/self.samples
            #~ P1= A+i*(C-A)/self.samples
            #~ for j in range(i+1):
                #~ if i!=0:
                    #~ P=P0+(P1-P0)*j/i
                #~ else:
                    #~ P=P0
                #~ Z=topo(P[0],P[1])
                #~ points.append((P[0],P[1],Z))
                #~ 
        #~ from matplotlib.delaunay import delaunay
        #~ 
        #~ #Need to find a beter way to do this not using delaunay# or maybe to generate all using triangulations????
        #~ 
        #~ x=[p[0] for p in points] 
        #~ y=[p[1] for p in points]
        #~ cs,e,trip,trin=delaunay(x,y)
        #~ return points, trip
        
    cpdef pointlist(self):
        
        cdef int i,j
        
        A=array(self.coord[0])
        B=array(self.coord[1])
        C=array(self.coord[2])
        
        #Get the mesh points
        X=[]
        Y=[]
        for i in range(self.samples+1):
            P0= A+i*(B-A)/self.samples
            P1= A+i*(C-A)/self.samples
            for j in range(i+1):
                if i!=0:
                    P=P0+(P1-P0)*j/i
                else:
                    P=P0
                X.append(P[0])
                Y.append(P[1])
        return X,Y

    cpdef limits(self):
        """
        Returns the minimum limits for the aperture
        """
        dx, dy=self.size
        return -dx/2, dx/2, -dy/2, dy/2
