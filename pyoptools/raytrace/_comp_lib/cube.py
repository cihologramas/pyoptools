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
# Description:     Beam Splitting Cube definition module
# Symbols Defined: BeamSplitingCube
# ------------------------------------------------------------------------------
#
"""
Definition of beam splitting cube object and helper functions
"""

import warnings
from numpy import pi

from pyoptools.raytrace.system import System
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.comp_lib import RightAnglePrism
from pyoptools.raytrace.mat_lib import Material
from pyoptools.raytrace.surface import Plane
from pyoptools.raytrace.shape import Rectangular


class Block(Component):
    """Class to define a Glass Block

    This class defines a component containing a glass block

    :param size: Tuple (W,H,L) containing the width, height and length of the
        glass block. Dimensions given in mm.
    :type size: tuple(float, float, float)

    """

    def __init__(self, size=(10, 10, 10), **traits):
        Component.__init__(self, **traits)
        self.size = size
        width, height, length = self.size
        __a_surf = Plane(shape=Rectangular(size=(width, height)))
        __p_surf = Plane(shape=Rectangular(size=(width, height)))

        __u_surf = Plane(shape=Rectangular(size=(width, length)))
        __l_surf = Plane(shape=Rectangular(size=(width, length)))

        __lf_surf = Plane(shape=Rectangular(size=(length, height)))
        __rg_surf = Plane(shape=Rectangular(size=(length, height)))

        self.surflist["S1"] = (__a_surf, (0, 0, -length / 2), (0, 0, 0))
        self.surflist["S2"] = (__p_surf, (0, 0, length / 2), (0, 0, 0))

        self.surflist["S3"] = (__u_surf, (0, height / 2, 0), (pi / 2, 0, 0))
        self.surflist["S4"] = (__l_surf, (0, -height / 2, 0), (pi / 2, 0, 0))

        self.surflist["S5"] = (__lf_surf, (-width / 2, 0, 0), (0, pi / 2, 0))
        self.surflist["S6"] = (__rg_surf, (width / 2, 0, 0), (0, pi / 2, 0))

    # ~ def __reduce__(self):
    # ~ args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
    # ~ return(type(self),args,self.__getstate__())

    # TODO: Check if there is a better way to do this, because we are
    # rewriting the constructor values here

    # ~ def __getstate__(self):
    # ~ return self.size, self.__a_surf, self.__p_surf, self.__u_surf,\
    # ~ self.__l_surf, self.__lf_surf, self.__rg_surf, \
    # ~ self.surflist
    # ~
    # ~ def __setstate__(self,state):
    # ~ self.size, self.__a_surf, self.__p_surf, self.__u_surf,\
    # ~ self.__l_surf, self.__lf_surf, self.__rg_surf, \
    # ~ self.surflist=state


class BeamSplittingCube(System):
    """Class to define a BeamSplittingCube.

    This class defines an System object containing the components to define an
    BeamSplitingCube.

    :param size: Side dimension of the cube
    :type size: float
    :param reflectivity: Reflectivity of the hypotenuse (between 0 and 1).
    :type reflectivity: float
    :param material: Material used to make the cube. Used to calculate the
       refraction index of the cube
    :type material: float or
       :class:`~pyoptools.raytrace.mat_lib.material.Material` subclass instance

    The origin of the coordinate system is located at the center of the cube
    in the optical axis (center of the hypotenuse).
    """

    def __init__(self, size=50.0, reflectivity=0.5, material=1.0, **traits):
        System.__init__(self, **traits)
        self.size = size
        self.reflectivity = reflectivity
        self.material = material

        __prism1 = RightAnglePrism(
            width=self.size,
            height=self.size,
            material=self.material,
            reflectivity=self.reflectivity,
        )

        __prism2 = RightAnglePrism(
            width=self.size, height=self.size, material=self.material, reflectivity=0
        )

        self.complist["C1"] = (__prism1, (0, 0, 0), (0, -pi / 2, 0))
        self.complist["C2"] = (__prism2, (0, 0, 0), (0, -3 * pi / 2, 0))

    # ~ def __reduce__(self):
    # ~ args=() #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
    # ~ return(type(self),args,self.__getstate__())
    # ~
    # ~
    # ~ #TODO: Check if there is a better way to do this, because we are
    # ~ #rewriting the constructor values here
    # ~
    # ~ def __getstate__(self):
    # ~ return self.size, self.reflectivity, self.material, self.__prism1,\
    # ~ self.__prism2, self.complist
    # ~
    # ~ def __setstate__(self,state):
    # ~ self.size, self.reflectivity, self.material, self.__prism1,\
    # ~ self.__prism2, self.complist = state


class BeamSplitingCube(BeamSplittingCube):
    """Deprecated class, please use the one with the correct spelling
    :class:`~pyoptools.raytrace.comp_lib.BeamSplittingCube`

    .. warning:: Will be removed in the future"""

    def __init__(self, *argv, **kwargs):

        warnings.warn(
            "This class will be deprecated (the spelling is not "
            "correct) . Please fix your code by using "
            "BeamSplittingCube instead",
            DeprecationWarning,
            stacklevel=2
        )
        super().__init__(*argv, **kwargs)
