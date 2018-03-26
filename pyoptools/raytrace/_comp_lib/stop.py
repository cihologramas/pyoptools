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
# Description:     Prism definitión module
# Symbols Defined: StopC
#------------------------------------------------------------------------------
#
'''
Definition of stop components
'''

#from enthought.traits.api import Float, Instance, HasTraits, Range
from numpy import sqrt, pi, absolute

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Aperture
from pyoptools.raytrace.shape import Shape

class Stop(Component):

    ''' **Class to define an stop component.**
      
    *Attributes:*
    
    *shape*
            External shape of the diaphragm
    *ap_shape*
            Aperture shape
    
    Note: The aperture shape must be contained by the external shape, but this 
    is not checked.
    '''

    # External shape of the diaphragm
    #shape=Instance(Shape)
    
    # Aperture shape
    #ap_shape=Instance(Shape)

    # Private attributes

    # Surfaces

    #__face = Instance(Aperture)

    def __init__(self,shape=None,ap_shape=None,**traits):
        Component.__init__(self,**traits)
        #self.shape=shape
        #self.ap_shape=ap_shape
        face= Aperture (shape=shape,  ap_shape=ap_shape)
        self.surflist["S1"]=(face,(0,0,0),(0,0,0))
        
    #~ def __reduce__(self):
        #~ args=(None,None) #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
        #~ return(type(self),args,self.__getstate__())
    #~ 
    #~ 
    #~ #TODO: Check if there is a better way to do this, because we are 
    #~ #rewriting the constructor values here
    #~ 
    #~ def __getstate__(self):
        #~ return self.shape, self.ap_shape, self.__face,self.surflist 
        #~ 
    #~ def __setstate__(self,state):
        #~ self.shape, self.ap_shape, self.__face,self.surflist=state
