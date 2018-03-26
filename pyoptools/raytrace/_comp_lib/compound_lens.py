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
# Description:     Compount Lens definitión module
# Symbols Defined: Doublet
#------------------------------------------------------------------------------
#
'''
Definition of compound lens objects and helper functions
'''

#from enthought.traits.api import Float, Instance, HasTraits, Trait
from numpy import pi

from pyoptools.raytrace.system import System
from pyoptools.raytrace.comp_lib import SphericalLens
from pyoptools.raytrace.mat_lib import Material

class Doublet(System):
    '''**Class to define a Doublet Lens**
    
    This class is used to define a System with containing the components needed
    to define a Doublet lens.
    
    *Attributes:*

    *radius*
        Radius of the doublet
    *curvature_s1*    
        Curvature of the anterior surface
    *curvature_s2*    
        Curvature of the middle surface
    *curvature_s3*    
        Curvature of the last surface
    *thickness_l1*    
        thickness of the anterior lens at the optical axis
    *thickness_l2*    
        thickness of the posterior lens at the optical axis
    *material_l1*     
        Material of the anterior lens
    *material_l2*     
        Material of the posterior lens

    The origin of the cordinate system is located at the center of the doublet
    in the optical axis.
    '''

    # radius of the doublet (aperture radius)
    #~ radius = Float(25.)
#~ 
    #~ # Curvature of the anterior surface
    #~ curvature_as= Float(0.01)
    #~ 
    #~ # Curvature of the midle surface
    #~ curvature_ms= Float(0.01)
    #~ 
    #~ # Curvature of the porterior surface
    #~ curvature_ps= Float(0.01)
#~ 
    #~ # Thickness of the anterior lens, measured at the optical axis
    #~ thickness_al= Float(5)
    #~ 
    #~ # Thickness of the posterior lens, measured at the optical axis
    #~ thickness_pl= Float(5)
#~ 
    #~ # Material of the anterior lens
    #~ material_al=Trait(1.,Float(),Material())
    #~ 
    #~ # Material of the posterior lens
    #~ material_pl=Trait(1.,Float(),Material())
#~ 
    #~ # Private attributes
#~ 
    #~ # Anterior lens, posterior lens
#~ 
    #~ __a_lens = Instance(SphericalLens)
    #~ __p_lens = Instance(SphericalLens)


    def __init__(self, radius = 25.,curvature_s1= 0.01,curvature_s2= 0.01,
                        curvature_s3= 0.01, thickness_l1= 5, thickness_l2= 5,
                        material_l1=1., material_l2=1.,*args,**kwarks):
        System.__init__(self,*args,**kwarks)
        self.radius=radius
        self.curvature_as= curvature_s1
        self.curvature_ms=curvature_s2
        self.curvature_ps=curvature_s3
        self.thickness_al=thickness_l1
        self.thickness_pl=thickness_l2
        self.material_al=material_l1
        self.material_pl=material_l2
        
        __a_lens= SphericalLens(curvature_s1=self.curvature_as,
                                     curvature_s2=self.curvature_ms,
                                     thickness=self.thickness_al,
                                     radius =self.radius,
                                     material=self.material_al)

        __p_lens= SphericalLens(curvature_s1=self.curvature_ms,
                                     curvature_s2=self.curvature_ps,
                                     thickness=self.thickness_pl,
                                     radius =self.radius,
                                     material=self.material_pl)

        self.complist["C1"]=(__a_lens,(0,0,-self.thickness_pl/2.),(0,0,0))
        self.complist["C2"]=(__p_lens,(0,0, self.thickness_al/2.),(0,0,0))

#    def __reduce__(self):
#        args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
#        return(type(self),args,self.__getstate__())
    
    
    #TODO: Check if there is a better way to do this, because we are 
    #rewriting the constructor values here
    
#    def __getstate__(self):
#        return self.radius, self.curvature_as, self.curvature_ms,\
#               self.curvature_ps, self.thickness_al, self.thickness_pl,\
#               self.material_al, self.material_pl, self.__a_lens, \
#               self.__p_lens, self.complist
        
#    def __setstate__(self,state):
#        self.radius, self.curvature_as, self.curvature_ms,\
#        self.curvature_ps, self.thickness_al, self.thickness_pl,\
#        self.material_al, self.material_pl, self.__a_lens, \
#        self.__p_lens, self.complist=state

class AirSpacedDoublet(System):
    '''**Class to define a an Air Spaced Doublet Lens**

    This class is used to define a System with containing the components needed
    to define a Doublet lens.

    *Attributes:*

    *radius*
        Radius of the doublet
    *curvature_s1*
        Curvature of the anterior surface of the first lens
    *curvature_s2*
        Curvature of the posterior surface of the first lens
    *curvature_s3*
        Curvature of the anterior surface of the first lens
    *curvature_s4*
        Curvature of the posterior surface of the first lens

    *thickness_l1*
        thickness of the anterior lens at the optical axis
    *thickness_l2*
        thickness of the posterior lens at the optical axis
    *air_gap*
        Distance between the 2 lenses
    *material_l1*
        Material of the anterior lens
    *material_l2*
        Material of the posterior lens

    The origin of the cordinate system is located at the center of the doublet
    in the optical axis.
    '''


    def __init__(self, radius = 25.,curvature_s1= 0.01,curvature_s2= 0.01,
                        curvature_s3= 0.01,curvature_s4= 0.01, thickness_l1= 5,air_gap=5 , thickness_l2= 5,
                        material_l1=1., material_l2=1.,*args,**kwarks):
        System.__init__(self,*args,**kwarks)
        self.radius=radius
        self.curvature_as1= curvature_s1
        self.curvature_ps1=curvature_s2

        self.curvature_as2=curvature_s3
        self.curvature_ps2=curvature_s4

        self.thickness_l1=thickness_l1
        self.thickness_l2=thickness_l2
        self.air_gap=air_gap

        self.material_l1=material_l1
        self.material_l2=material_l2

        __a_lens= SphericalLens(curvature_s1=self.curvature_as1,
                                     curvature_s2=self.curvature_ps1,
                                     thickness=self.thickness_l1,
                                     radius =self.radius,
                                     material=self.material_l1)

        __p_lens= SphericalLens(curvature_s1=self.curvature_as2,
                                     curvature_s2=self.curvature_ps2,
                                     thickness=self.thickness_l2,
                                     radius =self.radius,
                                     material=self.material_l2)

        self.complist["C1"]=(__a_lens,(0,0,-(self.thickness_l2+self.air_gap)/2.),(0,0,0))
        self.complist["C2"]=(__p_lens,(0,0, (self.thickness_l1+self.air_gap)/2.),(0,0,0))