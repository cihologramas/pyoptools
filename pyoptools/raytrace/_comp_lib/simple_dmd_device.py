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
# Description:     SimpleDMD Device component module
# Symbols Defined: SimpleDMDDevice
# ------------------------------------------------------------------------------

"""
Definition of SimpleDMD device component.

This module provides a complete DMD (Digital Micromirror Device) component
with a parallelepiped structure consisting of a SimpleDMD surface as the front
face and OpticalStop surfaces on the other five faces.
"""

from math import pi

from pyoptools.raytrace.component import Component
from pyoptools.raytrace.shape import Rectangular
from pyoptools.raytrace.surface import OpticalStop, SimpleDMD


class SimpleDMDDevice(Component):
    """Component representing a complete DMD device with parallelepiped structure.

    SimpleDMDDevice provides a convenient wrapper for integrating DMD (Digital
    Micromirror Device) elements into optical systems. The component consists of
    6 surfaces arranged as a parallelepiped:
    - Front face: SimpleDMD surface (active DMD with state-dependent reflection)
    - Other 5 faces: OpticalStop surfaces (full stops with no apertures)

    **Origin Pattern**: This component follows the "principal surface" pattern
    where the origin is placed at the center of the optically active surface
    (the DMD front face at z=0). This is consistent with mirrors and diffraction
    gratings, but differs from lenses which use the geometric center pattern.

    Parameters
    ----------
    tilt_angle : float
        Tilt angle in radians for the SimpleDMD surface. This is the angle between
        the surface normal and the perpendicular direction (Z-axis). When tilt_angle=0,
        the normal is perpendicular to the surface. Use math.radians() to convert
        from degrees.
    on_direction_angle : float
        Direction angle in radians defining where the surface normal is rotated
        to for the ON state. Measured counter-clockwise from the +X axis in the
        XY plane when viewed from +Z. For example, 0° points toward +X, 90° points
        toward +Y, 180° points toward -X, and 270° points toward -Y.
    off_direction_angle : float
        Direction angle in radians defining where the surface normal is rotated
        to for the OFF state. Measured counter-clockwise from the +X axis in the
        XY plane when viewed from +Z. For example, 0° points toward +X, 90° points
        toward +Y, 180° points toward -X, and 270° points toward -Y.
    state : str, optional
        Initial DMD state: "on", "off", or "flat". Default is "flat".
    width : float, optional
        Width of the DMD device in mm (X dimension). Default is 10.8 mm,
        matching typical Texas Instruments DLP chip active areas.
    height : float, optional
        Height of the DMD device in mm (Y dimension). Default is 10.8 mm.
    thickness : float, optional
        Thickness of the DMD device in mm (Z dimension). Default is 2.0 mm.
    **kwargs
        Additional arguments passed to Component base class.

    Attributes
    ----------
    surflist : dict
        Dictionary mapping surface names to (surface, position, rotation) tuples.
        Keys are: "front", "back", "left", "right", "top", "bottom".

    Examples
    --------
    Create a DMD device with 12° tilt magnitude::

        from math import radians
        from pyoptools.all import SimpleDMDDevice

        dmd = SimpleDMDDevice(
            tilt_angle=radians(12),           # 12° tilt magnitude
            on_direction_angle=radians(45),   # Normal tilts toward northeast
            off_direction_angle=radians(225), # Normal tilts toward southwest
            state="on"
        )

    Create a DMD with custom dimensions::

        dmd = SimpleDMDDevice(
            tilt_angle=radians(12),
            on_direction_angle=radians(45),
            off_direction_angle=radians(225),
            width=15.0,   # Custom width
            height=12.0,  # Custom height
            thickness=3.0 # Custom thickness
        )

    Change DMD state at runtime::

        # Access the SimpleDMD surface and change state
        dmd.surflist["front"][0].state = "off"

    Notes
    -----
    - **Default dimensions**: 10.368×5.832×1.0 mm (DLP4710 active area dimensions)
    - **Origin at DMD surface**: z=0 is at the center of the front (DMD) face,
      following the "principal surface" pattern used by mirrors/gratings
    - **SimpleDMD model**: Uses the simplified DMD surface model suitable for
      system-level design. Does not model individual pixels or diffraction effects.
    - **Side faces**: All five non-front faces are full stops (OpticalStop with
      ap_shape=None) that block rays completely.
    - **Direction angles**: These specify where the surface normal is rotated to,
      which may differ from manufacturer datasheets. Adjust empirically to match
      your desired optical behavior.

    See Also
    --------
    SimpleDMD : The underlying DMD surface class
    Component : Base class for optical components
    OpticalStop : Stop surface used for the side faces
    """

    def __init__(
        self,
        tilt_angle,
        on_direction_angle,
        off_direction_angle,
        state="flat",
        width=10.368,
        height=5.832,
        thickness=2.0,
        **kwargs,
    ):
        """Initialize SimpleDMDDevice component.

        Parameters described in class docstring.
        """
        Component.__init__(self, **kwargs)

        # Validate dimensions
        if width <= 0:
            raise ValueError(
                f"Width must be positive, got {width}. Specify width in millimeters."
            )
        if height <= 0:
            raise ValueError(
                f"Height must be positive, got {height}. Specify height in millimeters."
            )
        if thickness <= 0:
            raise ValueError(
                f"Thickness must be positive, got {thickness}. "
                "Specify thickness in millimeters."
            )

        # Store dimensions for __repr__
        self._width = width
        self._height = height
        self._thickness = thickness
        self._tilt_angle = tilt_angle
        self._on_direction_angle = on_direction_angle
        self._off_direction_angle = off_direction_angle
        self._state = state

        # Create front face - SimpleDMD surface at z=0 (origin at DMD surface)
        front_surface = SimpleDMD(
            tilt_angle=tilt_angle,
            on_direction_angle=on_direction_angle,
            off_direction_angle=off_direction_angle,
            state=state,
            shape=Rectangular(size=(width, height)),
        )
        self.surflist["front"] = (front_surface, (0, 0, 0), (0, pi, 0))

        # Create back face - OpticalStop at z=-thickness
        back_surface = OpticalStop(
            shape=Rectangular(size=(width, height)),
            ap_shape=None,  # Full stop, no aperture
        )
        self.surflist["back"] = (back_surface, (0, 0, thickness), (0, 0, 0))

        # Create left face - OpticalStop at x=-width/2
        left_surface = OpticalStop(
            shape=Rectangular(size=(thickness, height)), ap_shape=None
        )
        self.surflist["left"] = (
            left_surface,
            (-width / 2, 0, -thickness / 2),
            (0, pi / 2, 0),
        )

        # Create right face - OpticalStop at x=width/2
        right_surface = OpticalStop(
            shape=Rectangular(size=(thickness, height)), ap_shape=None
        )
        self.surflist["right"] = (
            right_surface,
            (width / 2, 0, -thickness / 2),
            (0, -pi / 2, 0),
        )

        # Create top face - OpticalStop at y=height/2
        top_surface = OpticalStop(
            shape=Rectangular(size=(width, thickness)), ap_shape=None
        )
        self.surflist["top"] = (
            top_surface,
            (0, height / 2, -thickness / 2),
            (-pi / 2, 0, 0),
        )

        # Create bottom face - OpticalStop at y=-height/2
        bottom_surface = OpticalStop(
            shape=Rectangular(size=(width, thickness)), ap_shape=None
        )
        self.surflist["bottom"] = (
            bottom_surface,
            (0, -height / 2, -thickness / 2),
            (pi / 2, 0, 0),
        )

    def __repr__(self):
        """Return string representation of SimpleDMDDevice.

        Returns
        -------
        str
            String showing component class, dimensions, and DMD state.
        """
        return (
            f"SimpleDMDDevice("
            f"width={self._width}, "
            f"height={self._height}, "
            f"thickness={self._thickness}, "
            f"state='{self._state}', "
            f"tilt_angle={self._tilt_angle:.5f}, "
            f"on_direction_angle={self._on_direction_angle:.5f}, "
            f"off_direction_angle={self._off_direction_angle:.5f})"
        )
