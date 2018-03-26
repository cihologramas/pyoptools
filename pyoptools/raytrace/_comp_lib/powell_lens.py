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
Definition of a spherical lens object and helper functions
"""

# from enthought.traits.api import Float, Instance, HasTraits,  Trait

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Aspherical, Plane
from pyoptools.raytrace.shape import Circular
from pyoptools.misc.Poly2D import *


class PowellLens(Component):
    """ **Class to define a powell Lens**.

    *Attributes:*

    *radius*
        diameter/2. of the lens in the part of the cylinder
    *thickness*
        Thicknes of the lens measured in the center
    *Conicity K*
        Conicity of the aspherical surface
    *curvature R*
        curvature of the aspherical surface
    *material*
        to calculate the refraction index of the lens (inherited from component)


    The origin of the cordinate systema is located at the center of the lens
    in the optical axis (center between vertex).
    """

    def __init__(self, radius=4.445, thickness=7.62, K=-4.302, R=3.00, *args, **kwargs):
        Component.__init__(self, *args, **kwargs)
        self.radius = radius
        self.thickness = thickness
        self.K = K
        self.R = R

        __a_surf = Aspherical(shape=Circular(radius=self.radius), Ax=0, Ay=self.R, Kx=-1, Ky=self.K,
                              poly=poly2d((0, 0)))

        self.surflist["S1"] = (__a_surf, (0, 0, 0), (0, 0, 0))

        __p_surf = Plane(shape=Circular(radius=(self.radius)))

        self.surflist["B1"] = (__p_surf, (0, 0, self.thickness), (0, 0, 0))
