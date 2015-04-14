#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Modulo con clases y funciones auxiliares.
"""



from pyoptools.raytrace.surface.idealsurface import IdealSurface
from pyoptools.raytrace.surface.idealpplanes import IdealPPlanes

from pyoptools.raytrace.shape import Rectangular
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface.aperture import Aperture
from pyoptools.raytrace.system.system import System

def IdealLens(shape=Rectangular(size=(50,50)), f=100):
    """Funcion envoltorio que representa una lente ideal
    """
    S1=IdealSurface(shape=shape, f=f)
    L1=Component(surflist=[(S1,(0,0,0),(0,0,0))])
    return L1


def IdealTLens(shape=Rectangular(size=(50,50)), ap_shape=Rectangular(size=(40,40)), f=100,d=20):
    """Funcion envoltorio que representa una lente ideal gruesa
    """
    S1=IdealPPlanes(shape=shape, f=f,d=d)
    S2=Aperture(shape=Rectangular(size=(50,50)),ap_shape=ap_shape)
    S3=Aperture(shape=Rectangular(size=(50,50)),ap_shape=ap_shape)
    A1=Component(surflist=[(S2,(0,0,0),(0,0,0)),])
    A2=Component(surflist=[(S3,(0,0,0),(0,0,0)),])
    L1=Component(surflist=[(S1,(0,0,0),(0,0,0)),])
    Si=System(complist=[(L1,(0,0,0),(0,0,0)),
                        (A1,(0,0,-1.001*d/2),(0,0,0)),
                        (A2,(0,0, 1.001*d/2),(0,0,0)),],n=1)

    return Si

