"""Module that defines all the classes that describe the optical surfaces"""

from .aspherical import Aspherical
from .cylinder import Cylinder
from .cylindrical import Cylindrical
from .detector import ArrayDetector
from .idealpplanes import IdealPPlanes
from .idealsurface import IdealSurface
from .opticalstop import Aperture, OpticalStop
from .plane import Plane
from .plane_mask import RPPMask
from .powell import Powell
from .simpledmd import SimpleDMD
from .spherical import Spherical
from .surface import Surface
from .taylor_poly import TaylorPoly

__all__ = [
    "Surface",
    "Cylindrical",
    "Plane",
    "Spherical",
    "ArrayDetector",
    "OpticalStop",
    "Aperture",
    "TaylorPoly",
    "Cylinder",
    "Aspherical",
    "Powell",
    "RPPMask",
    "IdealSurface",
    "IdealPPlanes",
    "SimpleDMD",
]
