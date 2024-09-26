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
# Symbols Defined: Polygon
# ------------------------------------------------------------------------------

from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.cmisc.eigen cimport Vector3d, Vector2d, \
    assign_tuple_to_vector2d

cdef class Polygon(Shape):

    """
    class defining a polygonal shape

     Args:
        coords (tuple): Tuple containing the coordinates of the 3 corners
            of a  triangle. Each coordinate is a(float, float) tuple.
        samples (int): Number of subdivitions per side used to sample the
            triangle.

    Todo:
        * This class is a copy of the Triangular class. Need to be implemented
          correctly.
   """

    def __init__(self, coord=((0, 0), (0, 100), (100, 0)), samples=10, *args, **kwargs):
        Shape.__init__(self, *args, **kwargs)

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
        """This method returns TRUE if an p=(x,y,z)point is inside the surface
        aperture if not it must return FALSE.
        This is implemented for a point, in cython, to make it fast
        """
        cdef double dot00, dot01, dot02, dot11, dot12, invDenom, u, v

        cdef double px, py

        px = point(0)
        py = point(1)

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

        cdef int i, j

        cdef Vector2d A = self.point_a
        cdef Vector2d B = self.point_b
        cdef Vector2d C = self.point_c
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
        Returns the minimum limits for the aperture
        """
        dx, dy=self.size
        return -dx/2, dx/2, -dy/2, dy/2
