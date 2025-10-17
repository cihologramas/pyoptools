# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Spherical lens definition module
# Symbols Defined: SphericalLens
# ------------------------------------------------------------------------------
#
"""
Definition of a Diffraction grating object and helper functions
"""

from math import pi, sin, cos

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Plane, RPPMask
from pyoptools.raytrace.shape import Rectangular
from pyoptools.misc.function_2d.poly_2d.poly_2d import Poly2D


class RectGratting(Component):
    def __init__(
        self,
        size=(50.0, 50.0, 10.0),
        reflectivity=0.0,
        lpmm=600,
        angle=0,
        M=[1],
        *args,
        **kwargs
    ):
        Component.__init__(self, *args, **kwargs)
        self.size = size
        width, height, length = self.size
        lpmmx = lpmm * cos(angle)
        lpmmy = lpmm * sin(angle)
        phf = Poly2D([0, 2 * pi * lpmmx, 2 * pi * lpmmy])
        __a_surf = RPPMask(
            shape=Rectangular(size=(width, height)), phm=phf, reflectivity=reflectivity, M=M
        )
        __p_surf = Plane(shape=Rectangular(size=(width, height)))

        __u_surf = Plane(shape=Rectangular(size=(width, length)))
        __l_surf = Plane(shape=Rectangular(size=(width, length)))

        __lf_surf = Plane(shape=Rectangular(size=(length, height)))
        __rg_surf = Plane(shape=Rectangular(size=(length, height)))

        self.surflist["S1"] = (__a_surf, (0, 0, 0), (0, 0, 0))
        self.surflist["S2"] = (__p_surf, (0, 0, length), (0, 0, 0))

        self.surflist["S3"] = (__u_surf, (0, height / 2, length / 2.0), (pi / 2, 0, 0))
        self.surflist["S4"] = (__l_surf, (0, -height / 2, length / 2.0), (pi / 2, 0, 0))

        self.surflist["S5"] = (__lf_surf, (-width / 2, 0, length / 2.0), (0, pi / 2, 0))
        self.surflist["S6"] = (__rg_surf, (width / 2, 0, length / 2.0), (0, pi / 2, 0))
