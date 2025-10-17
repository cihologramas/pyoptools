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
Definition of a spherical lens object and helper functions
"""


from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Cylindrical, Plane
from pyoptools.raytrace.shape import Rectangular


class CylindricalLens(Component):
    """Class to define a rectangular shaped cylindrical Lens.

    :param size: Size (sx,sy)  of the lens in mm
    :type size: tuple(float,float)
    :param thickness: Thickness of the lens at the center
    :type thickness: float
    :param curvature_s1: Curvature of the anterior surface of the lens, in mm.
    :type curvature_s1: float
    :param curvature_s2: Curvature of the posterior surface of the lens, in mm.
    :type curvature_s2: float
    :param material: Material of the lens
    :type material: float or
        :class:`~pyoptools.raytrace.mat_lib.material.Material`'s subclass
        instance

    .. warning::
       This component is composed only by the front and back surfaces. The
       surfaces closing the edges are not implemented yet. User must take
       this in to account and avoid rays that can enter or leave the lense by
       the edges.

    """

    def __init__(
        self,
        size=(20, 20),
        thickness=10,
        curvature_s1=1.0 / 200,
        curvature_s2=1.0 / 200,
        *args,
        **kwargs
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

        self.surflist["S1"] = (__a_surf, (0, 0, -self.thickness / 2), (0, 0, 0))
        self.surflist["S2"] = (__p_surf, (0, 0, self.thickness / 2), (0, 0, 0))

        # TODO: Falta cerrrar la lente por los lados
