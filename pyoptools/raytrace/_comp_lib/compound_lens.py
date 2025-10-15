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
# Description:     Compound Lens definition module
# Symbols Defined: Doublet
# ------------------------------------------------------------------------------
#
"""
Definition of compound lens objects and helper functions
"""

from math import isnan

from pyoptools.raytrace.system import System
from pyoptools.raytrace.comp_lib import SphericalLens
from pyoptools.raytrace.mat_lib import Material, material


class Doublet(System):
    """Class to define a Doublet Lens

    This class is used to define a System with containing the components needed
    to define a Doublet lens.

    :param radius: Radius of the doublet. Given in mm
    :type radius: float
    :param curvature_s1: Curvature of the anterior surface. Given in 1/mm
    :type curvature_s1: float
    :param curvature_s2: Curvature of the middle surface. Given in 1/mm
    :type curvature_s2: float
    :param curvature_s3: Curvature of the last surface
    :type curvature_s3: float
    :param thickness_l1: Thickness of the anterior lens at the optical axis
    :type thickness_l1: float
    :param thickness_l2: Thickness of the posterior lens at the optical axis
    :param material_l1: Material of the anterior lens
    :type material_l1: float or
                        :class:`~pyoptools.raytrace.mat_lib.material.Material`
                        subclass instance

    :param material_l2: Material of the posterior lens
    :type material_l1: float or
                        :class:`~pyoptools.raytrace.mat_lib.material.Material`
                        subclass instance

    The origin of the coordinate system is located at the center of the doublet
    in the optical axis.
    """

    def __init__(
        self,
        radius=25.0,
        curvature_s1=0.01,
        curvature_s2=0.01,
        curvature_s3=0.01,
        thickness_l1=5,
        thickness_l2=5,
        material_l1=1.0,
        material_l2=1.0,
        *args,
        **kwarks
    ):
        System.__init__(self, *args, **kwarks)
        self.radius = radius
        self.curvature_as = curvature_s1
        self.curvature_ms = curvature_s2
        self.curvature_ps = curvature_s3
        self.thickness_al = thickness_l1
        self.thickness_pl = thickness_l2
        self.material_al = material_l1
        self.material_pl = material_l2

        __a_lens = SphericalLens(
            curvature_s1=self.curvature_as,
            curvature_s2=self.curvature_ms,
            thickness=self.thickness_al,
            radius=self.radius,
            material=self.material_al,
        )

        __p_lens = SphericalLens(
            curvature_s1=self.curvature_ms,
            curvature_s2=self.curvature_ps,
            thickness=self.thickness_pl,
            radius=self.radius,
            material=self.material_pl,
        )

        self.complist["C1"] = (__a_lens, (0, 0, -self.thickness_pl / 2.0), (0, 0, 0))
        self.complist["C2"] = (__p_lens, (0, 0, self.thickness_al / 2.0), (0, 0, 0))


#    def __reduce__(self):
#        args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
#        return(type(self),args,self.__getstate__())


# TODO: Check if there is a better way to do this, because we are
# rewriting the constructor values here

#    def __getstate__(self):
#        return self.radius, self.curvature_as, self.curvature_ms,\
#               self.curvature_ps, self.thickness_al, self.thickness_pl,\
#               self.material_al, self.material_pl, self.__a_lens, \
#               self.__p_lens, self.complist

#    def __setstate__(self,state):
#        self.radius, self.curvature_as, self.curvature_ms,\
#        self.curvature_ps, self.thickness_al, self.thickness_pl,\
#        self.material_al, self.material_pl, self.__a_lens, \
#        self.__p_lens, self.complist=state


class AirSpacedDoublet(System):
    """Class to define a an Air Spaced Doublet Lens

    This class is used to define a System with containing the components needed
    to define a Doublet lens.

    :param radius: Radius of the doublet
    :type radius: float
    :param curvature_s1: Curvature of the anterior surface of the first lens
    :type curvature_s1: float
    :param curvature_s2: Curvature of the posterior surface of the first lens
    :type curvature_s2: float
    :param curvature_s3: Curvature of the anterior surface of the first lens
    :type curvature_s3: float
    :param curvature_s4: Curvature of the posterior surface of the first lens
    :type curvature_s4: float
    :param thickness_l1: Thickness of the anterior lens at the optical axis
    :type thickness_l1: float
    :param thickness_l2: Thickness of the posterior lens at the optical axis
    :type thickness_l2: float
    :param air_gap: Distance between the 2 lenses
    :type air_gap: float
    :param material_l1: Material of the anterior lens
    :type material_l1: float or
                        :class:`~pyoptools.raytrace.mat_lib.material.Material`
                        subclass instance
    :param material_l2: Material of the posterior lens
    :type material_l2: float or
                        :class:`~pyoptools.raytrace.mat_lib.material.Material`
                        subclass instance

    The origin of the coordinate system is located at the center of the doublet
    in the optical axis.
    """

    def __init__(
        self,
        radius=25.0,
        curvature_s1=0.01,
        curvature_s2=0.01,
        curvature_s3=0.01,
        curvature_s4=0.01,
        thickness_l1=5,
        air_gap=5,
        thickness_l2=5,
        material_l1=1.0,
        material_l2=1.0,
        *args,
        **kwarks
    ):
        System.__init__(self, *args, **kwarks)
        self.radius = radius
        self.curvature_as1 = curvature_s1
        self.curvature_ps1 = curvature_s2

        self.curvature_as2 = curvature_s3
        self.curvature_ps2 = curvature_s4

        self.thickness_l1 = thickness_l1
        self.thickness_l2 = thickness_l2
        self.air_gap = air_gap

        self.material_l1 = material_l1
        self.material_l2 = material_l2

        __a_lens = SphericalLens(
            curvature_s1=self.curvature_as1,
            curvature_s2=self.curvature_ps1,
            thickness=self.thickness_l1,
            radius=self.radius,
            material=self.material_l1,
        )

        __p_lens = SphericalLens(
            curvature_s1=self.curvature_as2,
            curvature_s2=self.curvature_ps2,
            thickness=self.thickness_l2,
            radius=self.radius,
            material=self.material_l2,
        )

        self.complist["C1"] = (
            __a_lens,
            (0, 0, -(self.thickness_l2 + self.air_gap) / 2.0),
            (0, 0, 0),
        )
        self.complist["C2"] = (
            __p_lens,
            (0, 0, (self.thickness_l1 + self.air_gap) / 2.0),
            (0, 0, 0),
        )


class MultiLens(System):
    """Class to define a multilens system from a table of parameters as given
    in standard raytracing programs.

    Attributes
    ----------
    sd : list
        List of tuples (ty, rad, thick, semid, matcat, marref)
        that contain each of the parameters needed to define
        each spherical surface using the same format used by the standard
        raytracing software.

        * ty: str
              String representing the surface type. For the moment only
              "spherical" surfaces are valid.
        * rad: float
              Radius of the surface in mm
        * thick: float
              Thickness of the material in mm. Distance from this surface to
              the next. For the last surface it has no real meaning.
        * semid: float
              Semi-diameter (radius) in mm of the round aperture that limits
              the surface.
        * matcat: str
               String with the material catalog where the material is defined.
               If matcat = "", the first material with name matref will be
               used.
        * matref: str of float
              String with the material name, or or number representing the
              constant refraction index of the material. matref = "" means
              there is no material between the 2 surfaces.


    The origin of the optical system is located in the center of the Multilens
    in the optical axis. The center is the mid-point between the 2 most
    external vertices
    """

    def __init__(self, sd, *args, **kwarks):
        System.__init__(self, *args, **kwarks)
        # Calculate the lens total thickness
        TT = 0
        for s in sd[:-1]:
            TT = TT + s[2]

        p = -TT / 2
        nn = 1
        for n in range(1, len(sd)):
            t0, r0, th0, s0, mc0, mt0 = sd[n - 1]
            t1, r1, th1, s1, mc1, mt1 = sd[n]

            if isnan(r0) or r0 == 0:
                c0 = 0
            else:
                c0 = 1 / r0

            if isnan(r1) or r1 == 0:
                c1 = 0
            else:
                c1 = 1 / r1

            if mt0 != "":
                if mc0 == "Value":
                    mat = float(mt0.replace(",", "."))
                else:
                    mat = getattr(material, mc0)[mt0]

                lens = SphericalLens(
                    curvature_s1=c0,
                    curvature_s2=c1,
                    thickness=th0,
                    radius=s0,
                    material=mat,
                )
                0, 0, p + th0 / 2

                self.complist["C{}".format(nn)] = (lens, (0, 0, p + th0 / 2), (0, 0, 0))
                nn = nn + 1
            p = p + th0
