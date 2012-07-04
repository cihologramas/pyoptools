#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Prism definitión module
# Symbols Defined: RightAnglePrism
#------------------------------------------------------------------------------
#
'''
Definition of a prism object and helper functions
'''

#from enthought.traits.api import Float, Instance, HasTraits, Range
from numpy import sqrt, pi, absolute

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Plane
from pyoptools.raytrace.shape import Rectangular,Triangular

class RightAnglePrism(Component):

    ''' **Class to define a Right Angle Prism Lens.**
      
    *Attributes:*
    
    *width*
            width of the prism face
    *height*
            height of the prism face
    *material*
            to calculate the refraction index of the prism (inerited from component)
    *reflectivity*
            reflectivity of the hipotenuse

    The origin of the cordinate systema is located at the center of hipotenuse
    face of the prism
    '''

    # Width of the prism
    #width = Float(50.)

    # Height of the prism
    #height = Float(10)

    #Reflectivity of the hipotenuse
    #reflectivity= Range(0.,1.)

    # Private attributes

    # Entrance and exit surfaces

    #__a_face = Instance(Plane)
    #__b_face = Instance(Plane)

    # Hipotenuse surface

    #__h_face = Instance(Plane)


    #Note the prism is not closed in the top or bottom needs to be fixed


    def __init__(self, width=50, height=10., reflectivity=0, *args, **kwargs):
        
        Component.__init__(self,*args, **kwargs)

        self.width=width
        self.height=height
        self.reflectivity=reflectivity
        
        
        __a_face= Plane (shape=Rectangular(size=(self.width,self.height)))
        __b_face= Plane (shape=Rectangular(size=(self.width,self.height)))

        h=sqrt(2.)*self.width

        __h_face= Plane (shape=Rectangular(size=(h,self.height)),reflectivity=self.reflectivity)
        
        w2=self.width/2.
        __e1=Plane (shape=Triangular(((-w2,w2),(-w2,-w2),(w2,-w2))))
        __e2=Plane (shape=Triangular(((-w2,w2),(-w2,-w2),(w2,-w2))))

        
        
        self.surflist["S1"]=(__a_face,(0,0,-self.width/2),(0,0,0))
        self.surflist["S2"]=(__b_face,(self.width/2,0,0),(0,pi/2,0))
        self.surflist["S3"]=(__h_face,(0,0,0),(0,-pi/4,0))
        self.surflist["S4"]=(__e1,(0,self.height/2,0),(pi/2,-pi/2,0) )
        self.surflist["S5"]=(__e2,(0,-self.height/2,0),(pi/2,-pi/2,0) )
        
    #~ def __reduce__(self):
        #~ args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
        #~ return(type(self),args,self.__getstate__())
    #~ 
    #~ 
    #~ #TODO: Check if there is a better way to do this, because we are 
    #~ #rewriting the constructor values here
    #~ 
    #~ def __getstate__(self):
        #~ return self.width, self.height, self.reflectivity, self.__a_face, \
               #~ self.__b_face, self.__h_face, self.__e1, self.__e2,\
               #~ self.surflist
        #~ 
    #~ def __setstate__(self,state):
        #~ self.width, self.height, self.reflectivity, self.__a_face, \
        #~ self.__b_face, self.__h_face, self.__e1, self.__e2, self.surflist=state
