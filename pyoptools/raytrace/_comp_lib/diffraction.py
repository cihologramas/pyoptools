#!/usr/bin/env python
# -*- coding: UTF-8 -*-
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
from pyoptools.misc.Poly2D.Poly2D import poly2d


class RectGratting(Component):
    def __init__(self, size=(50., 50., 10.), reflectivity=0., lpmm=600, angle=0, M=[1], *args, **kwargs):
        Component.__init__(self, *args, **kwargs)
        self.size = size
        w, h, l = self.size
        lpmmx = lpmm * cos(angle)
        lpmmy = lpmm * sin(angle)
        phf = poly2d([0, 2 * pi * lpmmx, 2 * pi * lpmmy])
        __a_surf = RPPMask(shape=Rectangular(size=(w, h)), phm=phf,
                           reflectivity=reflectivity, M=M)
        __p_surf = Plane(shape=Rectangular(size=(w, h)))

        __u_surf = Plane(shape=Rectangular(size=(w, l)))
        __l_surf = Plane(shape=Rectangular(size=(w, l)))

        __lf_surf = Plane(shape=Rectangular(size=(l, h)))
        __rg_surf = Plane(shape=Rectangular(size=(l, h)))

        self.surflist["S1"] = (__a_surf, (0, 0, 0), (0, 0, 0))
        self.surflist["S2"] = (__p_surf, (0, 0, l), (0, 0, 0))

        self.surflist["S3"] = (__u_surf, (0, h / 2, l / 2.), (pi / 2, 0, 0))
        self.surflist["S4"] = (__l_surf, (0, -h / 2, l / 2.), (pi / 2, 0, 0))

        self.surflist["S5"] = (__lf_surf, (-w / 2, 0, l / 2.), (0, pi / 2, 0))
        self.surflist["S6"] = (__rg_surf, (w / 2, 0, l / 2.), (0, pi / 2, 0))
