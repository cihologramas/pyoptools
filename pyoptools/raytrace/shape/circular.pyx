# ------------------------------------------------------------------------------
# Copyright (c) 2007,  2008, 2009 Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amezquita Orozco
# Description:     Rectangle definition module
# Symbols Defined: Circular
# ------------------------------------------------------------------------------

from libc.math cimport M_PI, sin, cos

from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.cmisc.eigen cimport Vector3d


cdef class Circular(Shape):
    """
    Class defining a circular shape.

    The `Circular` class represents a circular aperture or surface shape,
    defined by a given radius. It inherits from the `Shape` class and
    implements methods specific to circular shapes.

    Parameters
    ----------
    radius : float, optional
        The radius of the circular shape. Defaults to 1.0.
    samples : tuple of int, optional
        A tuple `(radial_samples, angular_samples)` that determines the number
        of samples used to discretize the circle. `radial_samples` defines the
        number of divisions along the radius, and `angular_samples` defines the
        number of divisions around the circle. Defaults to `(10, 36)`.
    *args : tuple, optional
        Additional positional arguments passed to the `Shape` superclass.
    **kwargs : dict, optional
        Additional keyword arguments passed to the `Shape` superclass.

    Attributes
    ----------
    radius : float
        The radius of the circular shape.
    samples : tuple of int
        A tuple `(radial_samples, angular_samples)` representing the sampling
        resolution of the circle.
    """

    def __init__(self, radius=1., samples=(10, 36)):
        Shape.__init__(self)
        self.radius = radius
        self.samples = samples
        # self.addkey("radius")
        # self.addkey("samples")

    def __reduce__(self):

        args = (self.radius, self.samples)
        return(type(self), args)

    cdef bint hit_cy(self, Vector3d& point) noexcept nogil:
        """
        Determine if a point is inside the circular surface aperture.

        This method checks whether a given point `p = (x, y, z)` lies inside the
        circular surface aperture. It returns `True` if the point is within the
        aperture and `False` otherwise. This method is implemented in Cython to
        ensure fast execution.

        Parameters
        ----------
        point : Vector3d&
            A reference to a `Vector3d` object representing the coordinates `(x, y, z)`
            of the point to be checked.

        Returns
        -------
        bint
            `True` if the point is within the surface aperture, `False` otherwise.
        """
        return point.dot(point) < self.radius * self.radius

    cpdef pointlist(self):
        """
        Generate a list of points that adequately sample the circular shape.

        This method returns two lists, `X` and `Y`, representing the X and Y
        coordinates of points that sample the circular shape. The sampling is
        based on the number of radial (`nr`) and angular (`na`) divisions specified
        in the `samples` attribute. The center point `(0, 0)` is always included
        in the list.

        Returns
        -------
        tuple of lists
            A tuple `(X, Y)` where `X` is a list of X coordinates and `Y` is a list
            of Y coordinates for the sampled points on the circular shape.

        Notes
        -----
        - The `samples` attribute defines the resolution of the sampling:
            - `nr` (number of radial samples) determines the number of points
            along the radius of the circle.
            - `na` (number of angular samples) determines the number of points
            around the circle's circumference.
        - The point `(0, 0)` is always included at the start of the lists.
        """
        cdef int nr, na, ir, ia

        nr, na = self.samples

        # Add the center point
        cdef list X = [0.,]
        cdef list Y = [0.,]

        cdef double angle, radius

        for ia in range(0, na):

            angle = ia*2*M_PI/na
            for ir in range(1, nr+1):
                radius = ir * self.radius / nr

                X.append(radius * cos(angle))
                Y.append(radius * sin(angle))

        return X, Y

    cpdef limits(self):
        """
        Return the minimum and maximum limits of the circular aperture.

        This method returns the minimum and maximum X and Y coordinates that
        define the bounding box of the circular aperture. These limits are
        determined by the radius of the circle.

        Returns
        -------
        tuple of float
            A tuple `(xmin, xmax, ymin, ymax)` where:
            - `xmin` is the minimum X-coordinate (`-self.radius`).
            - `xmax` is the maximum X-coordinate (`self.radius`).
            - `ymin` is the minimum Y-coordinate (`-self.radius`).
            - `ymax` is the maximum Y-coordinate (`self.radius`).

        Notes
        -----
        These limits define a square bounding box that fully contains the circular
        aperture.
        """
        return -self.radius, self.radius, -self.radius, self.radius
