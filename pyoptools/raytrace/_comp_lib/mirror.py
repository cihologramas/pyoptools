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
# Description:     Spherical lens definitión module
# Symbols Defined: SphericalLens
# ------------------------------------------------------------------------------
#
"""
Definition of a mirror object and helper functions
"""

from math import pi

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Cylindrical, Plane
from pyoptools.raytrace.shape import Circular, Rectangular


class RoundMirror(Component):
    def __init__(self, radius=50., thickness=10, reflectivity=0.5, *args, **kwargs):
        Component.__init__(self, *args, **kwargs)
        __a_surf = Plane(shape=Circular(radius=radius), reflectivity=reflectivity)
        __p_surf = Plane(shape=Circular(radius=radius))

        self.surflist["S1"] = (__a_surf, (0, 0, 0), (0, 0, 0))
        self.surflist["S2"] = (__p_surf, (0, 0, thickness), (0, 0, 0))

        __c_surf_1 = Cylindrical(shape=Rectangular(size=(2. * radius, thickness)),
                                 curvature=1. / radius)
        __c_surf_2 = Cylindrical(shape=Rectangular(size=(2 * radius, thickness)),
                                 curvature=1. / radius)

        self.surflist["B1"] = (__c_surf_1, (-radius, 0, thickness / 2.),
                               (pi / 2., 0, pi / 2))
        self.surflist["B2"] = (__c_surf_2, (radius, 0, thickness / 2.),
                               (-pi / 2., 0, pi / 2))


class RectMirror(Component):
    def __init__(self, size=(50., 50., 10.), reflectivity=0.5, *args, **kwargs):
        Component.__init__(self, *args, **kwargs)
        self.size = size
        w, h, l = self.size
        __a_surf = Plane(shape=Rectangular(size=(w, h)), reflectivity=reflectivity)
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
