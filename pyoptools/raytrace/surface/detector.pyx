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
# Description:     Detector surface definition module
# Symbols Defined: CCD
# ------------------------------------------------------------------------------

from numpy import zeros

from pyoptools.raytrace.surface.plane cimport Plane

from pyoptools.misc.pmisc import wavelength2RGB
from pyoptools.raytrace.shape.rectangular cimport Rectangular


cdef class ArrayDetector(Plane):
    """
    Class to simulate a CCD-like detector surface.

    The `ArrayDetector` class acts similarly to a real CCD (Charge-Coupled Device)
    in an optical system. It is designed to capture and store the coordinates
    of rays that hit the detector surface, effectively simulating the behavior
    of a CCD detector.

    Parameters
    ----------
    size : tuple of float, optional
        A tuple `(width, height)` representing the physical size of the detector
        surface in mm. The default size is `(10, 10)`.
    *args : tuple, optional
        Additional positional arguments passed to the `Plane` superclass.
    **kwargs : dict, optional
        Additional keyword arguments passed to the `Plane` superclass.

    Attributes
    ----------
    size : tuple of float
        The physical size of the detector surface, as specified during initialization.

    Examples
    --------
    Creating an ArrayDetector with a specified size::

        >>> detector = ArrayDetector(size=(10, 10))

    Notes
    -----
    - The `ArrayDetector` captures the coordinates of rays that hit the detector
      surface and can be used to simulate the detection process in optical simulations.
    - The `size` parameter defines the physical dimensions of the detector, which
      determines the area over which it can detect incoming rays.
    - Refer to the `Plane` class documentation for additional options and functionality
      inherited by the `ArrayDetector`.
    """

    # CCD physical size
    cdef public tuple size

    def __init__(self, size=(10, 10), *args, **kwargs):

        # Create a detector with a dummy size|
        Plane.__init__(self, shape=Rectangular(size=size), *args, **kwargs)
        self.size = size
        self.addkey("size")

    def get_histogram(self, resolution=(256, 256)):
        """
        Generate a histogram of ray impacts per unit area on the detector surface.

        This method returns a 2D array (histogram) representing the number of ray
        impacts per unit area of the detector. The histogram provides a pixelated
        view of how rays are distributed across the detector surface.

        Parameters
        ----------
        resolution : tuple of int, optional
            A tuple `(px, py)` representing the resolution of the histogram in pixels.
            `px` is the number of pixels along the X-axis, and `py` is the number
            of pixels along the Y-axis. The default value is `(256, 256)`. The physical
            size of the detector is determined during the surface's creation.

        Returns
        -------
        ndarray
            A 2D numpy array of shape `(px, py)` where each element represents the
            number of ray impacts in that pixel area.

        Notes
        -----
        - The `resolution` parameter controls the resolution of the histogram.
        Increasing the resolution results in a finer grid with more pixels,
        while decreasing the resolution results in a coarser grid with fewer
        pixels.
        - The physical size of the detector, defined during its creation,
        determines the actual area over which the impacts are distributed.
        - The method assumes that the Z-coordinate of all impacts is zero,
        consistent with a planar detector surface.
        - Rays that hit the detector are recorded in the `_hit_list`, which stores
        their impact coordinates.
        """
        cdef double x, y, _z, sx, sy
        cdef int px, py, nx, ny
        px, py = resolution
        sx, sy = self.size

        retval = zeros((px, py))

        # Iterate over the hit list to accumulate ray impacts
        for i in self._hit_list:
            x, y, _z = i[0]  # z should always be 0
            nx = <int>(px*(x+sx/2.)/sx)
            ny = <int>(py*(y+sy/2.)/sy)
            retval[nx, ny] += 1
        return retval

    def get_color_histogram(self, resolution=(256, 256)):
        """
        Generate a color histogram of ray impacts per unit area on the detector surface.

        This method returns a 3D array representing the number of ray impacts per
        unit area of the detector, simulating a color image. Each impact contributes
        to the RGB color channels based on the ray's wavelength.

        Parameters
        ----------
        resolution : tuple of int, optional
            A tuple `(px, py)` representing the resolution of the histogram in pixels.
            `px` is the number of pixels along the X-axis, and `py` is the number
            of pixels along the Y-axis. The default value is `(256, 256)`. The physical
            size of the detector is determined during the surface's creation.

        Returns
        -------
        ndarray
            A 3D numpy array of shape `(px, py, 3)` where each element contains the
            accumulated RGB values corresponding to ray impacts at that pixel position.

        Notes
        -----
        - The `resolution` parameter controls the resolution of the histogram.
        Increasing the resolution results in a finer grid with more pixels,
        while decreasing the resolution results in a coarser grid with fewer
        pixels.
        - The physical size of the detector, defined during its creation, determines
        the actual area over which the impacts are distributed.
        - The method does not currently account for the intensity of the rays, which
        should be considered in future improvements.
        - The Z-coordinate of the impact is assumed to be zero, consistent with a
        planar detector surface.
        - Rays that hit the detector are recorded in the `_hit_list`, which stores
        their impact coordinates and wavelengths.
        """
        cdef int px, py, nx, ny
        cdef double sx, sy
        px, py = resolution
        sx, sy = self.size

        retval = zeros((px, py, 3))
        for i in self._hit_list:
            x, y, _z = i[0]
            # z must be always 0
            nx = <int>(px*(x+sx/2.)/sx)
            ny = <int>(py*(y+sy/2.)/sy)

            # TODO: The intensities are not taken into account. This should be fixed
            r, g, b = wavelength2RGB(i[1].wavelength)
            retval[nx, ny, 0] += r
            retval[nx, ny, 1] += g
            retval[nx, ny, 2] += b
        return retval
