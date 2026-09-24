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
# Description:     Prism definition module
# Symbols Defined: RightAnglePrism
# ------------------------------------------------------------------------------
#
"""
Definition of a prism object and helper functions
"""

from math import cos, radians, sin

from numpy import pi, sqrt

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.shape import Polygon, Rectangular, Triangular
from pyoptools.raytrace.surface import Plane


class RightAnglePrism(Component):
    """Class to define a Right Angle Prism.

    **ARGUMENTS:**

        ============ ===========================================================
        width        Width of the prism face
        height       Height of the prism face
        material     To calculate the refraction index of the prism (inherited
                     from component)
        reflectivity Reflectivity of the coating of the hypotenuse. For a normal
                     prism it is 0. Note: Total internal reflection works in the
                     prism.
        reflega      Reflectivity of the Leg A of the prism. For a normal prism
                     it is 0.
        reflegb      Reflectivity of the Leg B of the prism. For a normal prism
                     it is 0.
        ============ ===========================================================

    The origin of the coordinate system is located at the center of hypotenuse
    face of the prism
    """

    def __init__(
        self,
        width=50,
        height=10.0,
        reflectivity=0,
        reflega=0,
        reflegb=0,
        *args,
        **kwargs,
    ):
        Component.__init__(self, *args, **kwargs)

        self.width = width
        self.height = height
        self.reflectivity = reflectivity
        self.reflega = reflega
        self.reflegb = reflegb

        __a_face = Plane(
            shape=Rectangular(size=(self.width, self.height)), reflectivity=reflega
        )
        __b_face = Plane(
            shape=Rectangular(size=(self.width, self.height)), reflectivity=reflegb
        )

        h = sqrt(2.0) * self.width

        __h_face = Plane(
            shape=Rectangular(size=(h, self.height)), reflectivity=self.reflectivity
        )

        w2 = self.width / 2.0
        __e1 = Plane(shape=Triangular(((-w2, w2), (-w2, -w2), (w2, -w2))))
        __e2 = Plane(shape=Triangular(((-w2, w2), (-w2, -w2), (w2, -w2))))

        self.surflist["S1"] = (__a_face, (0, 0, -self.width / 2), (0, 0, 0))
        self.surflist["S2"] = (__b_face, (self.width / 2, 0, 0), (0, pi / 2, 0))
        self.surflist["S3"] = (__h_face, (0, 0, 0), (0, -pi / 4, 0))
        self.surflist["S4"] = (__e1, (0, self.height / 2, 0), (pi / 2, -pi / 2, 0))
        self.surflist["S5"] = (__e2, (0, -self.height / 2, 0), (pi / 2, -pi / 2, 0))


class PentaPrism(Component):
    """Class to define a pentaprism

    :param s: Size of the entrance aperture of the pentaprism (the aperture is
        square)
    :type s: float
    :param material: Material of the pentaprism
    :type material: float or
        :class:`~pyoptools.raytrace.mat_lib.material.Material`'s subclass
        instance
    """

    def __init__(self, s, *args, **kwargs):

        Component.__init__(self, *args, **kwargs)

        s1 = Plane(shape=Rectangular(size=(s, s)))
        s2 = Plane(shape=Rectangular(size=(s, s)))
        d = s / cos(radians(22.5))
        s3 = Plane(shape=Rectangular(size=(d, s)), reflectivity=1)
        s4 = Plane(shape=Rectangular(size=(d, s)), reflectivity=1)
        d1 = d * sin(radians(22.5) / 2.0)
        s5 = Plane(shape=Rectangular(size=(2 * sqrt(2) * d1, s)))

        half_s = s / 2.0
        # 5 vertices of the pentaprism cross-section in XZ
        p1 = (half_s, -half_s)
        p2 = (-half_s, -half_s)
        p3 = (-half_s - d1 - (d / 2.0) * sin(radians(22.5)), half_s)
        p4 = (-half_s, half_s + d1 + (d / 2.0) * sin(radians(22.5)))
        p5 = (half_s, half_s)
        poly_coords = (p1, p2, p3, p4, p5)

        s6 = Plane(shape=Polygon(coord=poly_coords))
        s7 = Plane(shape=Polygon(coord=poly_coords))

        self.surflist["S1"] = (s1, (0, 0, -s / 2.0), (0, 0, 0))
        self.surflist["S2"] = (s2, (s / 2.0, 0, 0), (0, pi / 2, 0))
        self.surflist["S3"] = (s3, (0, 0, s / 2.0 + d1), (0, pi / 8, 0))
        self.surflist["S4"] = (s4, (-s / 2.0 - d1, 0, 0), (0, 3 * pi / 8, 0))
        self.surflist["S5"] = (s5, (-s / 2.0 - d1, 0, s / 2.0 + d1), (0, -pi / 4, 0))
        self.surflist["S6"] = (s6, (0, half_s, 0), (pi / 2, 0, 0))
        self.surflist["S7"] = (s7, (0, -half_s, 0), (pi / 2, 0, 0))


class DovePrism(Component):
    """Class to define a dove prism

    :param s: Height and depth of the dove prism
    :type s: float
    :param length: Width of the dove prism (length of the longest side)
    :type length: float
    :param material: Material of the prism
    :type material: float or
        :class:`~pyoptools.raytrace.mat_lib.material.Material`'s subclass
        instance
    """

    def __init__(self, s, length, *args, **kwargs):
        # s alto o profundidad del prisma
        # length Longitud del lado mas largo del prisma

        # La referencia del prisma de dove esta en el centro del prisma

        Component.__init__(self, *args, **kwargs)

        d = 1.4142135623730951 * s

        # Diagonales del prisma

        s1 = Plane(shape=Rectangular(size=(d, s)))
        s2 = Plane(shape=Rectangular(size=(d, s)))

        # Lado largo del prisma
        s3 = Plane(shape=Rectangular(size=(length, s)))
        # lado corto del prisma
        s4 = Plane(shape=Rectangular(size=(length - 2 * s, s)))

        half_l = length / 2.0
        top_half_l = half_l - s
        half_s = s / 2.0
        # Trapezoid vertices in XZ plane
        poly_coords = (
            (-half_l, -half_s),
            (half_l, -half_s),
            (top_half_l, half_s),
            (-top_half_l, half_s),
        )
        s5 = Plane(shape=Polygon(coord=poly_coords))
        s6 = Plane(shape=Polygon(coord=poly_coords))

        sp1 = (length - s) / 2.0
        self.surflist["S1"] = (s1, (-sp1, 0, 0), (0, -pi / 4, 0))
        self.surflist["S2"] = (s2, (sp1, 0, 0), (0, pi / 4, 0))
        self.surflist["S3"] = (s3, (0, 0, -s / 2.0), (0, 0, 0))
        self.surflist["S4"] = (s4, (0, 0, s / 2.0), (0, 0, 0))
        self.surflist["S5"] = (s5, (0, half_s, 0), (pi / 2, 0, 0))
        self.surflist["S6"] = (s6, (0, -half_s, 0), (pi / 2, 0, 0))
