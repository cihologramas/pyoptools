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
# Symbols Defined: StopC
# ------------------------------------------------------------------------------
#
"""
Definition of stop components
"""

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import OpticalStop


class Stop(Component):
    """Class to define an stop component.

    This component is used to simulate apertures or diafragms in an optical
    system

    **ARGUMENTS**

        ======== =====================
        shape    Stop external shape
        ap_shape Aperture (hole) shape
        ======== =====================

    **shape** and **ap_shape** are instances of any sub-class of
    :class:`~pyoptools.raytrace.shape.Shape`.

    .. warning::
        The aperture shape must be contained by the external shape, but this
        is not checked.
    """

    # External shape of the diaphragm
    # shape=Instance(Shape)

    # Aperture shape
    # ap_shape=Instance(Shape)

    # Private attributes

    # Surfaces

    # __face = Instance(Aperture)

    def __init__(self, shape=None, ap_shape=None, **traits):
        Component.__init__(self, **traits)
        # self.shape=shape
        # self.ap_shape=ap_shape
        face = OpticalStop(shape=shape, ap_shape=ap_shape)
        self.surflist["S1"] = (face, (0, 0, 0), (0, 0, 0))
