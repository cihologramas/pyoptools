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
# Description:     Polygon shape definition module
# Symbols Defined: Polygon
# ------------------------------------------------------------------------------

from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.cmisc.eigen cimport Vector3d, Vector2d


cdef class Polygon(Shape):
    """Class defining an arbitrary polygonal aperture shape.

    Parameters
    ----------
    coord : sequence of (float, float), optional
        List or tuple of 2D vertex coordinates defining the polygon perimeter
        in counter-clockwise or clockwise order. Minimum 3 vertices.
        Defaults to `((0, 0), (0, 100), (100, 0))`.
    samples : int, optional
        Number of subdivisions for internal mesh sampling. Defaults to 10.
    """

    def __init__(self, coord=((0, 0), (0, 100), (100, 0)), samples=10, *args, **kwargs):
        Shape.__init__(self, *args, **kwargs)

        if len(coord) < 3:
            raise ValueError("A polygon must have at least 3 vertices.")

        self._coord = tuple((float(p[0]), float(p[1])) for p in coord)
        self.samples = samples

        self.poly_points.clear()
        cdef Vector2d pt
        for p in self._coord:
            pt = Vector2d(p[0], p[1])
            self.poly_points.push_back(pt)

    @property
    def coord(self):
        """Return the coordinates of the polygon vertices."""
        return self._coord

    def __reduce__(self):
        args = (self._coord, self.samples)
        return (type(self), args)

    cdef bint hit_cy(self, Vector3d &point) noexcept nogil:
        """Return True if point (x, y, z) lies within the polygon aperture."""
        cdef double px = point(0)
        cdef double py = point(1)
        cdef int n = self.poly_points.size()
        if n < 3:
            return False

        cdef bint inside = False
        cdef int i, j = n - 1
        cdef double xi, yi, xj, yj

        for i in range(n):
            xi = self.poly_points[i](0)
            yi = self.poly_points[i](1)
            xj = self.poly_points[j](0)
            yj = self.poly_points[j](1)

            if ((yi > py) != (yj > py)) and (
                px < (xj - xi) * (py - yi) / (yj - yi) + xi
            ):
                inside = not inside
            j = i

        return inside

    cpdef limits(self):
        """Return the bounding box (xmin, xmax, ymin, ymax) of the polygon."""
        cdef int n = self.poly_points.size()
        if n == 0:
            return 0.0, 0.0, 0.0, 0.0

        cdef double xmin = self.poly_points[0](0)
        cdef double xmax = xmin
        cdef double ymin = self.poly_points[0](1)
        cdef double ymax = ymin
        cdef double x, y
        cdef int i

        for i in range(1, n):
            x = self.poly_points[i](0)
            y = self.poly_points[i](1)
            if x < xmin:
                xmin = x
            if x > xmax:
                xmax = x
            if y < ymin:
                ymin = y
            if y > ymax:
                ymax = y

        return xmin, xmax, ymin, ymax

    cpdef pointlist(self):
        """Return lists (X, Y) of mesh points for surface representation."""
        cdef list X = []
        cdef list Y = []
        cdef int n = self.poly_points.size()
        cdef int i, j

        for i in range(n):
            X.append(self.poly_points[i](0))
            Y.append(self.poly_points[i](1))

        cdef double xmin, xmax, ymin, ymax
        xmin, xmax, ymin, ymax = self.limits()

        cdef int samples = max(self.samples, 2)
        cdef double dx = (xmax - xmin) / <double>samples
        cdef double dy = (ymax - ymin) / <double>samples
        cdef Vector3d test_pt
        cdef double px, py

        if dx > 0 and dy > 0:
            for i in range(1, samples):
                px = xmin + i * dx
                for j in range(1, samples):
                    py = ymin + j * dy
                    test_pt = Vector3d(px, py, 0.0)
                    if self.hit_cy(test_pt):
                        X.append(px)
                        Y.append(py)

        return X, Y
