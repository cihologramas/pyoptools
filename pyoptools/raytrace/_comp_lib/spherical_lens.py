#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

# from enthought.traits.api import Float, Instance, HasTraits,  Trait
from numpy import sqrt, pi, absolute, inf

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Spherical, Cylindrical, Plane
from pyoptools.raytrace.shape import Circular, Rectangular


class SphericalLens(Component):
    """Helper class to define spherical lenses.

    **ARGUMENTS:**
        ============ ==========================================================
        radius       Aperture radious of the lens given in mm
        thickness    Thickness of the lens measured in the center given in mm
        curvature_s1 curvature of the anterior surface given in 1/mm
        curvature_s2 curvature of the posterior surface given in 1/mm
        material     Used to calculate the refraction index of the lens
                     (inherited from component). Can be a material instance or
                     floating point number if the refraction index is constant.
        ============ ==========================================================

    **RETURN VALUE:**

    Returns a Component subclass that behaves as a spherical lens. The origin of
    the Component's coordinate system is located at the center of the lens
    in the optical axis (mid-point between vertex of the Spherical surfaces).
    """

    # Aperture of the lens
    # radius = Float(50.)

    # Thickness of the lens at the optical axis
    # thickness= Float(10)

    # Curvature of the anterior surface
    # curvature_s1 = Float( 1./200)

    # Curvature of the posterior surface
    # curvature_s2 = Float( 1./200)

    # Private attributes

    # anterior surface
    # __a_surf = Trait(None, Instance(Spherical), Instance(Plane))

    # posterior surface
    # __p_surf = Trait(None, Instance(Spherical), Instance(Plane))

    # cylindrical surfaces to close the lens
    # __c_surf_1= Instance(Cylindrical)
    # __c_surf_2= Instance(Cylindrical)

    def __init__(
        self,
        radius=50.0,
        thickness=10,
        curvature_s1=1.0 / 200,
        curvature_s2=1.0 / 200,
        *args,
        **kwargs
    ):
        Component.__init__(self, *args, **kwargs)
        self.radius = radius
        self.thickness = thickness
        self.curvature_s1 = curvature_s1
        self.curvature_s2 = curvature_s2

        if self.radius * abs(self.curvature_s1) > 1.0:
            raise ValueError("Aperture radius can not be larger than s1 curvature radius")

        if self.radius * abs(self.curvature_s2) > 1.0:
            raise ValueError("Aperture radius can not be larger than s2 curvature radius")

        if self.curvature_s1 != 0.0:
            __a_surf = Spherical(
                shape=Circular(radius=self.radius), curvature=self.curvature_s1
            )
        else:
            __a_surf = Plane(shape=Circular(radius=self.radius))

        if self.curvature_s2 != 0:
            __p_surf = Spherical(
                shape=Circular(radius=self.radius), curvature=self.curvature_s2
            )
        else:
            __p_surf = Plane(shape=Circular(radius=self.radius))

        self.surflist["S1"] = (__a_surf, (0, 0, -self.thickness / 2), (0, 0, 0))
        self.surflist["S2"] = (__p_surf, (0, 0, self.thickness / 2), (0, 0, 0))

        if self.curvature_s1 != 0:
            r_a = 1.0 / self.curvature_s1
            s_a = absolute(r_a) - sqrt(r_a * r_a - self.radius * self.radius)
            if (r_a) < 0:
                s_a = -s_a
        else:
            s_a = 0.0

        if self.curvature_s2 != 0:
            r_p = 1.0 / self.curvature_s2
            s_p = absolute(r_p) - sqrt(r_p * r_p - self.radius * self.radius)
            if (r_p) > 0:
                s_p = -s_p
        else:
            s_p = 0.0

        # Ojo, falta verificar si la lente es fisicamente posible es decir th1>0
        th1 = self.thickness - s_a - s_p

        zp = float(-self.thickness / 2 + s_a + th1 / 2.0)
        __c_surf_1 = Cylindrical(
            shape=Rectangular(size=(2.0 * self.radius, th1)),
            curvature=1.0 / self.radius,
        )

        __c_surf_2 = Cylindrical(
            shape=Rectangular(size=(2 * self.radius, th1)), curvature=1.0 / self.radius
        )

        self.surflist["B1"] = (__c_surf_1, (-self.radius, 0, zp), (pi / 2.0, 0, pi / 2))

        self.surflist["B2"] = (__c_surf_2, (self.radius, 0, zp), (-pi / 2.0, 0, pi / 2))

    # ~ def __reduce__(self):
    # ~ args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
    # ~ return(type(self),args,self.__getstate__())
    # ~
    # ~
    # ~ #TODO: Check if there is a better way to do this, because we are
    # ~ #rewriting the constructor values here
    # ~
    # ~ def __getstate__(self):
    # ~
    # ~ return self.radius, self.thickness, self.curvature_s1, self.curvature_s2, \
    # ~ self.__a_surf, self.__p_surf, self.surflist
    # ~
    # ~
    # ~ def __setstate__(self,state):
    # ~ self.radius, self.thickness, self.curvature_s1, self.curvature_s2, \
    # ~ self.__a_surf, self.__p_surf, self.surflist = state

    def paraxial_constants(self, wavelength=0.58929, n=1.0):
        """Method to calculate the paraxial constants of a spherical lens

        **ARGUMENT:**

            ========== ============================================
            wavelength Wavelength used for the calculations
            n          Refraction index of the surrounding media
            ========== ============================================

        **RETURN VALUE:**

            3 element tuple (f, afl, pfl) containing

            === =====================================================
            f   Effective Focal length
            afl Anterior focal length  (negative for positive lenses)
            pfl Posterior focal length (positive for positive lenses)
            === =====================================================
        """

        nl = self.material.n(wavelength)

        # Anterior surface focal length (measured inside the lens)

        asf = inf if self.curvature_s1 == 0 else nl / ((nl - n) * self.curvature_s1)

        # Posterior surface focal length (measured inside the lens)

        psfp = inf if self.curvature_s2 == 0 else nl / ((n - nl) * self.curvature_s2)

        # Posterior surface focal length (measured outside the lens)

        psf = inf if self.curvature_s2 == 0 else n / ((n - nl) * self.curvature_s2)

        # Focal Length (it is the same at both sides because the media
        # on both sides is the same)

        f = n / (nl / asf + n / psf - self.thickness * n / (asf * psf))

        # Anterior Focal length

        afl = -f * (1.0 - self.thickness / (psfp))

        # Posterior Focal length

        pfl = f * (1.0 - self.thickness / asf)

        return (f, afl, pfl)
