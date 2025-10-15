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

from numpy import sqrt, pi, absolute

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Plane
from pyoptools.raytrace.shape import Rectangular, Triangular
from math import cos, radians, sin


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
        **kwargs
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

    # ~ def __reduce__(self):
    # ~ args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
    # ~ return(type(self),args,self.__getstate__())
    # ~
    # ~
    # ~ #TODO: Check if there is a better way to do this, because we are
    # ~ #rewriting the constructor values here
    # ~
    # ~ def __getstate__(self):
    # ~ return self.width, self.height, self.reflectivity, self.__a_face, \
    # ~ self.__b_face, self.__h_face, self.__e1, self.__e2,\
    # ~ self.surflist
    # ~
    # ~ def __setstate__(self,state):
    # ~ self.width, self.height, self.reflectivity, self.__a_face, \
    # ~ self.__b_face, self.__h_face, self.__e1, self.__e2, self.surflist=state


class PentaPrism(Component):
    """Class to define a pentaprism

    :param s: Size of the entrance aperture of the pentaprism (the aperture is
        square)
    :type s: float
    :param material: Material of the pentaprism
    :type material: float or
        :class:`~pyoptools.raytrace.mat_lib.material.Material`'s subclass
        instance

    .. warning::

        The pentaprism has no upper or lower surface. Care must be taken to
        avoid rays entering or exiting by such apertures.

        .. todo:: Fix this

    """

    # TODO: El pentaprisma está abierto por arriba y por abajo. Hay que definir las superficies para cerrarlo

    def __init__(self, s, *args, **kwargs):

        Component.__init__(self, *args, **kwargs)

        s1 = Plane(shape=Rectangular(size=(s, s)))
        s2 = Plane(shape=Rectangular(size=(s, s)))
        d = s / cos(radians(22.5))
        s3 = Plane(shape=Rectangular(size=(d, s)), reflectivity=1)
        s4 = Plane(shape=Rectangular(size=(d, s)), reflectivity=1)
        d1 = d * sin(radians(22.5) / 2.0)
        s5 = Plane(shape=Rectangular(size=(2 * sqrt(2) * d1, s)))

        self.surflist["S1"] = (s1, (0, 0, -s / 2.0), (0, 0, 0))
        self.surflist["S2"] = (s2, (s / 2.0, 0, 0), (0, pi / 2, 0))
        self.surflist["S3"] = (s3, (0, 0, s / 2.0 + d1), (0, pi / 8, 0))
        self.surflist["S4"] = (s4, (-s / 2.0 - d1, 0, 0), (0, 3 * pi / 8, 0))
        self.surflist["S5"] = (s5, (-s / 2.0 - d1, 0, s / 2.0 + d1), (0, -pi / 4, 0))
        # ,material=get_material("N-BK7"))


class DovePrism(Component):
    """Class to define a dove prism

    :param s: Height and depth of the dove prism
    :type s: float
    :param l: Width of the dove prism (length of the longest side)
    :type l: float
    :param material: Material of the prism
    :type material: float or
        :class:`~pyoptools.raytrace.mat_lib.material.Material`'s subclass
        instance

    .. warning::

        The pentaprism has no upper or lower surface. Care must be taken to
        avoid rays entering or exiting by such apertures.

        .. todo:: Fix this
    """

    # TODO: El prismadove está abierto por arriba y por abajo. Hay que definir las superficies para cerrarlo
    def __init__(self, s, length, *args, **kwargs):
        # s alto o profundidad del prisma
        # l Longitud del lado mas largo del prisma

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

        sp1 = (length - s) / 2.0
        self.surflist["S1"] = (s1, (-sp1, 0, 0), (0, -pi / 4, 0))
        self.surflist["S2"] = (s2, (sp1, 0, 0), (0, pi / 4, 0))
        self.surflist["S3"] = (s3, (0, 0, -s / 2.0), (0, 0, 0))
        self.surflist["S4"] = (s4, (0, 0, s / 2.0), (0, 0, 0))
