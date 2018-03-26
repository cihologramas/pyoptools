#!/usr/bin/env python
# -*- coding: UTF-8 -*-
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
# Description:     Beam Splitting Cube definitión module
# Symbols Defined: BeamSplitingCube
#------------------------------------------------------------------------------
#
'''
Definition of beam splitting cube object and helper functions
'''

#from enthought.traits.api import Tuple, Float, Instance, HasTraits, Trait
from numpy import pi

from pyoptools.raytrace.system import System
from pyoptools.raytrace.component import Component 
from pyoptools.raytrace.comp_lib import RightAnglePrism
from pyoptools.raytrace.mat_lib import Material
from pyoptools.raytrace.surface import Plane
from pyoptools.raytrace.shape import Rectangular

class Block(Component):
    '''**Class to define a Glass Block**
    
    This class defines a component containing a glass block
    
    *Attributes:*
    
    *size*
        Tuple containing the width, height and length of the glass
        block
    
    '''
    
    #Tuple containing the width, height and length of the glass block
    #size=Tuple(Float(10), Float(10), Float(10))
    
    #__a_surf = Trait(None,Instance(Plane))
    #__p_surf = Trait(None,Instance(Plane))
    #__u_surf = Trait(None,Instance(Plane))
    #__l_surf = Trait(None,Instance(Plane))
    #__lf_surf = Trait(None,Instance(Plane))
    #__rg_surf = Trait(None,Instance(Plane))
    
    def __init__(self, size=(10,10,10), **traits):
        Component.__init__(self,**traits)
        self.size=size
        w, h, l=self.size
        __a_surf = Plane(shape=Rectangular(size=(w, h)))
        __p_surf = Plane(shape=Rectangular(size=(w, h)))
        
        __u_surf = Plane(shape=Rectangular(size=(w, l)))
        __l_surf = Plane(shape=Rectangular(size=(w, l)))
        
        __lf_surf = Plane(shape=Rectangular(size=(l, h)))
        __rg_surf = Plane(shape=Rectangular(size=(l, h)))
        
        self.surflist["S1"]=(__a_surf,(0,0,-l/2),(0,0,0))
        self.surflist["S2"]=(__p_surf,(0,0,l/2),(0,0,0))
        
        self.surflist["S3"]=(__u_surf,(0,h/2,0),(pi/2,0,0))
        self.surflist["S4"]=(__l_surf,(0,-h/2,0),(pi/2,0,0))
        
        self.surflist["S5"]=(__lf_surf,(-w/2,0,0),(0,pi/2,0))
        self.surflist["S6"]=(__rg_surf,(w/2,0,0),(0,pi/2,0))

    #~ def __reduce__(self):
        #~ args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
        #~ return(type(self),args,self.__getstate__())
    
    
    #TODO: Check if there is a better way to do this, because we are 
    #rewriting the constructor values here
    
    #~ def __getstate__(self):
        #~ return self.size, self.__a_surf, self.__p_surf, self.__u_surf,\
               #~ self.__l_surf, self.__lf_surf, self.__rg_surf, \
               #~ self.surflist
                       #~ 
    #~ def __setstate__(self,state):
        #~ self.size, self.__a_surf, self.__p_surf, self.__u_surf,\
        #~ self.__l_surf, self.__lf_surf, self.__rg_surf, \
        #~ self.surflist=state        
        
        
class BeamSplitingCube(System):
    '''**Class to define a BeamSplitingCube.**
    
    This class defines an System object containing the components to define an
    BeamSplitingCube

    *Attributes:*

    *size*
                Side of the cube
    *reflectivity*
                Reflectivity of the hypothenuse
    *material*
                to calculate the refraction index of the cube

    The origin of the cordinate systema is located at the center of the cube
    in the optical axis (center of the hypotenuse).
    '''


    def __init__(self, size=50.,reflectivity=0.5,material=1., **traits):
        System.__init__(self,**traits)
        self.size=size
        self.reflectivity=reflectivity
        self.material=material
            
    
        
        __prism1= RightAnglePrism(width=self.size,height=self.size,
                                        material=self.material,
                                        reflectivity=self.reflectivity)

        __prism2= RightAnglePrism(width=self.size,height=self.size,
                                        material=self.material,reflectivity=0)

        self.complist["C1"]=(__prism1,(0,0,0),(0,-pi/2,0))
        self.complist["C2"]=(__prism2,(0,0,0),(0,-3*pi/2,0))
        
    #~ def __reduce__(self):
        #~ args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
        #~ return(type(self),args,self.__getstate__())
    #~ 
    #~ 
    #~ #TODO: Check if there is a better way to do this, because we are 
    #~ #rewriting the constructor values here
    #~ 
    #~ def __getstate__(self):
        #~ return self.size, self.reflectivity, self.material, self.__prism1,\
               #~ self.__prism2, self.complist
        #~ 
    #~ def __setstate__(self,state):
        #~ self.size, self.reflectivity, self.material, self.__prism1,\
        #~ self.__prism2, self.complist = state
