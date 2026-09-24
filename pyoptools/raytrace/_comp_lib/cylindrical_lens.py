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
# Description:     Cylindrical lens definition module
# Symbols Defined: CylindricalLens
# ------------------------------------------------------------------------------
#
"""
Definition of a cylindrical lens object and helper functions
"""

from math import pi

import numpy as np

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.shape import Rectangular
from pyoptools.raytrace.surface import Cylindrical, Plane


class CylindricalLens(Component):
    """Class to define a rectangular shaped cylindrical Lens.

    :param size: Size (sx,sy) of the lens in mm
    :type size: tuple(float,float)
    :param thickness: Thickness of the lens at the center
    :type thickness: float
    :param curvature_s1: Curvature of the anterior surface of the lens, in mm^-1.
    :type curvature_s1: float
    :param curvature_s2: Curvature of the posterior surface of the lens, in mm^-1.
    :type curvature_s2: float
    :param material: Material of the lens
    :type material: float or
        :class:`~pyoptools.raytrace.mat_lib.material.Material`'s subclass
        instance
    """

    def __init__(
        self,
        size=(20, 20),
        thickness=10,
        curvature_s1=1.0 / 200,
        curvature_s2=1.0 / 200,
        *args,
        **kwargs,
    ):
        Component.__init__(self, *args, **kwargs)
        self.size = size
        w, h = self.size
        self.thickness = thickness
        self.curvature_s1 = curvature_s1
        self.curvature_s2 = curvature_s2

        if self.curvature_s1 != 0.0:
            __a_surf = Cylindrical(
                shape=Rectangular(size=(w, h)), curvature=self.curvature_s1
            )
        else:
            __a_surf = Plane(shape=Rectangular(size=(w, h)))

        if self.curvature_s2 != 0:
            __p_surf = Cylindrical(
                shape=Rectangular(size=(w, h)), curvature=self.curvature_s2
            )
        else:
            __p_surf = Plane(shape=Rectangular(size=(w, h)))

        self.surflist["S1"] = (__a_surf, (0, 0, -self.thickness / 2.0), (0, 0, 0))
        self.surflist["S2"] = (__p_surf, (0, 0, self.thickness / 2.0), (0, 0, 0))

        def _sag(c, x):
            if c == 0:
                return 0.0
            r = 1.0 / c
            if abs(x) > abs(r):
                raise ValueError(
                    f"Aperture half-width {abs(x)} exceeds curvature radius {abs(r)}"
                )
            return r - np.sign(r) * np.sqrt(r * r - x * x)

        s_a = _sag(self.curvature_s1, w / 2.0)
        s_p = _sag(self.curvature_s2, w / 2.0)

        z_s1_min = min(-self.thickness / 2.0, -self.thickness / 2.0 + s_a)
        z_s1_max = max(-self.thickness / 2.0, -self.thickness / 2.0 + s_a)
        z_s2_min = min(self.thickness / 2.0, self.thickness / 2.0 + s_p)
        z_s2_max = max(self.thickness / 2.0, self.thickness / 2.0 + s_p)

        edge_thickness = (self.thickness / 2.0 + s_p) - (-self.thickness / 2.0 + s_a)
        if edge_thickness <= 0:
            raise ValueError(
                "Lens parameters result in non-physical zero or negative edge thickness"
            )

        z_edge_center = (
            (-self.thickness / 2.0 + s_a) + (self.thickness / 2.0 + s_p)
        ) / 2.0

        z_min = min(z_s1_min, z_s2_min)
        z_max = max(z_s1_max, z_s2_max)
        z_total_thickness = z_max - z_min
        z_center = (z_max + z_min) / 2.0

        # Surfaces closing the lens edges
        # Upper and lower surfaces (y = +h/2, y = -h/2)
        __u_surf = Plane(shape=Rectangular(size=(w, z_total_thickness)))
        __l_surf = Plane(shape=Rectangular(size=(w, z_total_thickness)))
        self.surflist["S3"] = (__u_surf, (0, h / 2.0, z_center), (pi / 2.0, 0, 0))
        self.surflist["S4"] = (__l_surf, (0, -h / 2.0, z_center), (pi / 2.0, 0, 0))

        # Left and right surfaces (x = -w/2, x = +w/2)
        __lf_surf = Plane(shape=Rectangular(size=(edge_thickness, h)))
        __rg_surf = Plane(shape=Rectangular(size=(edge_thickness, h)))
        self.surflist["S5"] = (
            __lf_surf,
            (-w / 2.0, 0, z_edge_center),
            (0, pi / 2.0, 0),
        )
        self.surflist["S6"] = (
            __rg_surf,
            (w / 2.0, 0, z_edge_center),
            (0, pi / 2.0, 0),
        )
