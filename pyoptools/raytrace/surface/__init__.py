#!/usr/bin/env python
# -*- coding: UTF-8 -*-

'''Module that defines all the clases that describe the optical surfaces
'''

from .surface import Surface
from .cylindrical import Cylindrical
from .plane import Plane
from .spherical import Spherical
from .detector import ArrayDetector
from .aperture import Aperture
from .taylor_poly import TaylorPoly
from .cylinder import Cylinder
from .aspherical import Aspherical
from .powell import Powell
from .plane_mask import RPPMask
from .idealsurface import IdealSurface
from .idealpplanes import IdealPPlanes

__all__ = ["Surface",
           "Cylindrical",
           "Plane",
           "Spherical",
           "ArrayDetector",
           "Aperture",
           "TaylorPoly",
           "Cylinder",
           "Aspherical",
           "Powell",
           "RPPMask",
           "IdealSurface",
           "IdealPPlanes"]
