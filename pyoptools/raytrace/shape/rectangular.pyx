# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amezquita Orozco
# Description:     Rectangle definition module
# Symbols Defined: Rectangle
# ------------------------------------------------------------------------------
#

"""Module that defines the Rectangular class
"""


from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.cmisc.eigen cimport Vector3d

cdef class Rectangular(Shape):
    """
    Class defining a rectangular shape.

    The `Rectangular` class represents a rectangular aperture or surface shape,
    defined by its size. It inherits from the `Shape` class and implements
    methods specific to rectangular shapes.

    Parameters
    ----------
    size : tuple of float, optional
        A tuple `(width, height)` that defines the size of the rectangle.
        Defaults to `(1.0, 1.0)`.
    samples : tuple of int, optional
        A tuple `(nx, ny)` that determines the number of samples used to
        discretize the rectangle along its width and height, respectively.
        Defaults to `(30, 30)`.
    offset : tuple of float, optional
        A tuple `(x_offset, y_offset)` that specifies the offset of the
        rectangle from the origin. Defaults to `(0.0, 0.0)`.

    Attributes
    ----------
    size : tuple of float
        A tuple `(width, height)` representing the dimensions of the rectangle.
    samples : tuple of int
        A tuple `(nx, ny)` representing the number of samples along the width
        and height of the rectangle.
    offset : tuple of float
        A tuple `(x_offset, y_offset)` representing the offset of the rectangle
        from the origin.
    """
    def __init__(self, size=(1., 1.), samples=(30, 30), offset=(0, 0)):
        Shape.__init__(self)
        self.size = (float(size[0]), float(size[1]))
        self.samples = (int(samples[0]), int(samples[1]))
        self.offset = (float(offset[0]), float(offset[1]))
        # self.addkey("size")
        # self.addkey("samples")
        # self.addkey("offset")

    cdef bint hit_cy(self, Vector3d& point) noexcept nogil:
        """
        Determine if a point is inside the rectangular surface aperture.

        This method checks whether a given point `p = (x, y, z)` lies within the
        boundaries of the rectangular aperture. It returns `True` if the point is
        inside the aperture and `False` otherwise. This method is implemented in
        Cython to ensure fast execution.

        Parameters
        ----------
        point : Vector3d&
            A reference to a `Vector3d` object representing the coordinates `(x, y, z)`
            of the point to be checked.

        Returns
        -------
        bint
            `True` if the point is within the rectangular surface aperture,
            `False` otherwise.

        """
        cdef double dx, dy, ox, oy, opx, opy
        cdef double px, py

        px = point(0)
        py = point(1)
        # pz = point(2)

        dx, dy = self.size
        ox, oy = self.offset
        opx = px-ox
        opy = py-oy

        return (opx >= -dx/2.) and (opx <= dx/2.) and (opy >= -dy/2.) and (opy <= dy/2.)

    cpdef pointlist(self):
        """
        Generate a list of points that adequately sample the rectangular shape.

        This method returns two lists, `X` and `Y`, representing the X and Y
        coordinates of points that sample the rectangular shape. The sampling
        resolution is determined by the number of divisions along the width
        and height of the rectangle, specified by the `samples` attribute.

        The method explicitly calculates the grid points, taking into account
        the size of the rectangle and its offset from the origin.

        Returns
        -------
        tuple of lists
            A tuple `(X, Y)` where `X` is a list of X coordinates and `Y` is a
            list of Y coordinates for the sampled points on the rectangular shape.

        Notes
        -----
        - The `samples` attribute defines the resolution of the sampling:
            - `nx` (number of samples along the width) determines the number of
            points along the X-axis.
            - `ny` (number of samples along the height) determines the number of
            points along the Y-axis.
        - The points are calculated such that the grid spans from `-dx/2 + ox`
        to `dx/2 + ox` in the X direction and from `-dy/2 + oy` to `dy/2 + oy`
        in the Y direction, where `dx` and `dy` are the width and height of
        the rectangle, and `ox` and `oy` are the offsets.
        """

        cdef double dx, dy, ox, oy
        cdef int nx, ny
        cdef list X, Y
        cdef int i, j
        cdef double step_x, step_y, x_val, y_val

        dx, dy = self.size
        nx, ny = self.samples
        ox, oy = self.offset

        # Calculate step sizes
        step_x = dx / (nx - 1) if nx > 1 else 0
        step_y = dy / (ny - 1) if ny > 1 else 0

        # Initialize empty lists for X and Y
        X = []
        Y = []

        # Create the grid of points
        for j in range(ny):
            y_val = (-dy / 2.0) + oy + j * step_y
            for i in range(nx):
                x_val = (-dx / 2.0) + ox + i * step_x
                X.append(x_val)
                Y.append(y_val)

        return X, Y

    cpdef limits(self):
        """
        Return the minimum and maximum limits of the rectangular aperture.

        This method returns the minimum and maximum X and Y coordinates that
        define the bounding box of the rectangular aperture. These limits are
        calculated based on the size of the rectangle and its offset from the
        origin.

        Returns
        -------
        tuple of float
            A tuple `(xmin, xmax, ymin, ymax)` where:
            - `xmin` is the minimum X-coordinate, calculated as `-dx/2 + ox`.
            - `xmax` is the maximum X-coordinate, calculated as `dx/2 + ox`.
            - `ymin` is the minimum Y-coordinate, calculated as `-dy/2 + oy`.
            - `ymax` is the maximum Y-coordinate, calculated as `dy/2 + oy`.

        Notes
        -----
        - `dx` and `dy` represent the width and height of the rectangle,
        respectively.
        - `ox` and `oy` represent the X and Y offsets of the rectangle from
        the origin.
        - These limits define a rectangular bounding box that fully contains
        the aperture.
        """
        cdef float dx, dy, ox, oy
        dx, dy = self.size
        ox, oy = self.offset
        return -dx/2+ox, dx/2+ox, -dy/2+oy, dy/2+oy
