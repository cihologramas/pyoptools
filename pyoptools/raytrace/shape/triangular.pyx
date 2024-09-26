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
# Symbols Defined: Triangle
# ------------------------------------------------------------------------------
#

"""Module that defines the Triangular class
"""


from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.cmisc.eigen cimport Vector3d, Vector2d, \
    assign_tuple_to_vector2d


cdef class Triangular(Shape):

    """
    Class defining a triangular polygonal shape.

    This class represents a triangular shape defined by the coordinates of
    its three corners. It inherits from the `Shape` class and implements
    methods specific to triangular shapes.

    Parameters
    ----------
    coord : tuple of tuple of float, optional
        A tuple containing the coordinates of the three corners of the
        triangle. Each corner is represented by a `(x, y)` tuple.
        Defaults to `((0, 0), (0, 100), (100, 0))`.
    samples : int, optional
        The number of subdivisions per side used to sample the triangle.
        This determines the resolution of the grid points within the triangle.
        Defaults to `10`.

    Attributes
    ----------
    point_a : Vector2d
        A `Vector2d` object representing the first corner of the triangle.
    point_b : Vector2d
        A `Vector2d` object representing the second corner of the triangle.
    point_c : Vector2d
        A `Vector2d` object representing the third corner of the triangle.
    samples : int
        The number of subdivisions per side used to sample the triangle.
    """

    def __init__(self, coord=((0, 0), (0, 100), (100, 0)), samples=10):
        Shape.__init__(self)

        # self.coord=coord
        assign_tuple_to_vector2d(coord[0], self.point_a)
        assign_tuple_to_vector2d(coord[1], self.point_b)
        assign_tuple_to_vector2d(coord[2], self.point_c)

        self.samples=samples

        # Register picklable attributes
        # self.addkey("point_a")
        # self.addkey("point_b")
        # self.addkey("point_c")

        # self.addkey("samples")

    def __reduce__(self):

        args=(self.coord, self.samples)
        return(type(self), args)

    cdef bint hit_cy(self, Vector3d &point) noexcept nogil:
        """
        Determine if a point is inside the triangular surface aperture.

        This method checks whether a given point `p = (x, y, z)` lies within
        the boundaries of the triangular surface aperture. It returns `True`
        if the point is inside the triangle and `False` otherwise. This method
        is implemented in Cython to ensure fast execution.

        The algorithm used is based on barycentric coordinates, which are computed
        to determine if the point lies inside the triangle formed by the corners
        `point_a`, `point_b`, and `point_c`.

        Parameters
        ----------
        point : Vector3d&
            A reference to a `Vector3d` object representing the coordinates `(x, y, z)`
            of the point to be checked. Only the x and y coordinates are used in
            the calculation.

        Returns
        -------
        bint
            `True` if the point is within the triangular surface aperture,
            `False` otherwise.

        """

        cdef double dot00, dot01, dot02, dot11, dot12, invDenom, u, v

        cdef double px, py, _pz

        px = point(0)
        py = point(1)
        _pz = point(2)

        cdef Vector2d P = Vector2d(px, py)  # = array((px, py))
        # A=array(self.coord[0])
        # B=array(self.coord[1])
        # C=array(self.coord[2])

        cdef Vector2d v0 = self.point_c - self.point_a  # C-A
        cdef Vector2d v1 = self.point_b - self.point_a  # B-A
        cdef Vector2d v2 = P - self.point_a  # P-A

        dot00=v0.dot(v0)  # dot(v0, v0)
        dot01=v0.dot(v1)  # dot(v0, v1)
        dot02=v0.dot(v2)  # dot(v0, v2)
        dot11=v1.dot(v1)  # dot(v1, v1)
        dot12=v1.dot(v2)  # dot(v1, v2)

        invDenom=1./(dot00 * dot11 - dot01 * dot01)

        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom

        # Check if point is in triangle
        return (u > 0) and (v > 0) and (u + v < 1)

    cpdef pointlist(self):
        """
        Generate a list of points that adequately sample the triangular shape.

        This method returns two lists, `X` and `Y`, representing the X and Y
        coordinates of points that sample the triangular shape defined by the
        vertices `point_a`, `point_b`, and `point_c`. The sampling resolution
        is determined by the `samples` attribute, which specifies the number of
        subdivisions per side of the triangle.

        Returns
        -------
        tuple of lists
            A tuple `(X, Y)` where `X` is a list of X coordinates and `Y` is a
            list of Y coordinates for the sampled points within the triangular shape.

        Notes
        -----
        - The method generates a triangular grid of points by linearly interpolating
        between the vertices `point_a`, `point_b`, and `point_c`.
        - The `samples` attribute determines the number of subdivisions along
        each side of the triangle.
        - The outer loop iterates over the number of samples, creating points along
        the edges of the triangle.
        - The inner loop interpolates between the points on the two sides of the
        triangle to fill in the interior points.
        - The method handles the edge case when `i == 0` to ensure that the starting
        vertex is correctly included.
        """

        cdef int i, j

        cdef Vector2d A=self.point_a
        cdef Vector2d B=self.point_b
        cdef Vector2d C=self.point_c
        cdef Vector2d P0, P1, P
        # Get the mesh points
        cdef list X=[]
        cdef list Y=[]

        for i in range(self.samples+1):
            P0= A+((B-A)*<double>i)/<double>self.samples
            P1= A+((C-A)*<double>i)/<double>self.samples
            for j in range(i+1):
                if i!=0:
                    P=P0+(P1-P0)*(<double>j/i)
                else:
                    P=P0
                X.append(P(0))
                Y.append(P(1))
        return X, Y

    cpdef limits(self):
        """
        Return the minimum and maximum limits of the triangular aperture.

        This method returns the minimum and maximum X and Y coordinates that
        define the bounding box of the triangular aperture. These limits are
        calculated based on the coordinates of the three vertices of the triangle.

        Returns
        -------
        tuple of float
            A tuple `(xmin, xmax, ymin, ymax)` where:
            - `xmin` is the minimum X-coordinate among the three vertices.
            - `xmax` is the maximum X-coordinate among the three vertices.
            - `ymin` is the minimum Y-coordinate among the three vertices.
            - `ymax` is the maximum Y-coordinate among the three vertices.

        Notes
        -----
        - The limits are calculated directly from the coordinates of the triangle's
        vertices (`point_a`, `point_b`, and `point_c`).
        - The returned limits define a rectangular bounding box that fully contains
        the triangular shape.
        """
        cdef double xmin, xmax, ymin, ymax

        # Extract coordinates from the vertices
        cdef double ax, ay, bx, by, cx, cy
        ax, ay = self.point_a(0), self.point_a(1)
        bx, by = self.point_b(0), self.point_b(1)
        cx, cy = self.point_c(0), self.point_c(1)

        # Calculate bounding box limits
        xmin = min(ax, bx, cx)
        xmax = max(ax, bx, cx)
        ymin = min(ay, by, cy)
        ymax = max(ay, by, cy)

        return xmin, xmax, ymin, ymax
