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
# Description:     Spherical surface definition module
# Symbols Defined: Spherical
# ------------------------------------------------------------------------------

from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_nan_to_vector3d, \
                                        assign_doubles_to_vector3d

from libc.math cimport sqrt, abs

cdef class Spherical(Surface):
    """
    Class to define spherical optical surfaces.

    The `Spherical` class represents a spherical optical surface. The surface
    is defined by specifying the shape of the aperture and the radius of
    curvature of the sphere.

    The vertex of the spherical surface is located at the origin of coordinates
    (0, 0, 0), and the aperture is centered at the origin as well.

    Parameters
    ----------
    curvature : double, optional
        The curvature of the spherical surface. A positive curvature indicates
        a convex surface, while a negative curvature indicates a concave surface.
        The default is `0.0`.
    *args : tuple, optional
        Additional positional arguments passed to the `Surface` superclass.
    **kwargs : dict, optional
        Additional keyword arguments passed to the `Surface` superclass.

    Attributes
    ----------
    curvature : double
        The curvature of the spherical surface, as specified during initialization.

    Examples
    --------
    Creating a spherical surface with a circular aperture and specific curvature::

        >>> cs = Spherical(shape=Circular(radius=60), curvature=0.15)

    Notes
    -----
    - The `curvature` is defined as the reciprocal of the radius of curvature
      (i.e., `curvature = 1 / radius`). A flat surface would have a curvature of `0`.
    - Refer to the `Surface` documentation for additional options and parameters
      that can be used when defining the surface.
    """
    cdef public double curvature

    def __init__(self, curvature=0., *args, **kwargs):
        Surface.__init__(self, *args, **kwargs)
        self.curvature = curvature
        self.addkey("curvature")

    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:
        """
        Calculate the intersection point between a ray and the spherical surface.

        This method calculates the point of intersection between the surface
        and the given ray. The intersection point is calculated in the coordinate
        system of the surface.

        Parameters
        ----------
        incident_ray : Ray
            The incoming ray for which the intersection with the spherical surface
            is to be calculated.
        intersection_point : Vector3d&
            A reference to a `Vector3d` object where the calculated intersection
            point will be stored. If there is no valid intersection, this vector
            is set to NaN.

        Notes
        -----
        - This method assumes the spherical surface is centered at the origin.
        - The method only handles a spherical surface that represents a half-sphere.
        If the calculated intersection point lies outside this half-sphere,
        the intersection point will be set to NaN.
        - This method uses the quadratic formula to find the intersection points
        and selects the valid one based on the direction of propagation.
        - There may be issues if the ray propagates along the X or Y direction
        (i.e., if the ray is perpendicular to the Z-axis). This can lead
        to numerical instability, and a better solution may be required.

        """

        cdef double x1, y1, z1, x21, y21, z21, z3, a, b, c, b2ac
        cdef double u1, u2, X1, Y1, Z1, X2, Y2, Z2

        x1 = incident_ray._origin(0)
        y1 = incident_ray._origin(1)
        z1 = incident_ray._origin(2)

        x21 = incident_ray._direction(0)
        y21 = incident_ray._direction(1)
        z21 = incident_ray._direction(2)

        z3 = 1./self.curvature

        a = x21**2+y21**2+z21**2
        b = 2*((x21)*(x1) + (y21)*(y1) + (z21)*(z1 - z3))
        c = x1**2+y1**2+z1**2-2*(z3*z1)
        b2ac = b**2-4*a*c

        if b2ac <= 0.:
            assign_nan_to_vector3d(intersection_point)
        else:
            u1 = (-b+sqrt(b2ac))/(2*a)
            u2 = (-b-sqrt(b2ac))/(2*a)

            X1 = x1 + u1 * (x21)
            Y1 = y1 + u1 * (y21)
            Z1 = z1 + u1 * (z21)

            X2 = x1 + u2 * (x21)
            Y2 = y1 + u2 * (y21)
            Z2 = z1 + u2 * (z21)

            if abs(Z2) < abs(Z1):
                X, Y, Z = X2, Y2, Z2
            else:
                X, Y, Z = X1, Y1, Z1

            if abs(Z) > abs(z3):
                assign_nan_to_vector3d(intersection_point)
            else:
                assign_doubles_to_vector3d(X, Y, Z, intersection_point)

    cdef void _calculate_normal(self,
                                Vector3d& intersection_point,
                                Vector3d& normal) noexcept nogil:
        """
        Calculate the normal vector to the spherical surface at a given point.

        This method calculates the normal vector to the spherical surface at the
        specified intersection point. The normal is computed in the coordinate
        system of the surface.

        Parameters
        ----------
        intersection_point : Vector3d&
            A reference to a `Vector3d` object representing the point on the
            spherical surface where the normal vector is to be calculated.
        normal : Vector3d&
            A reference to a `Vector3d` object where the calculated normal vector
            will be stored.

        Notes
        -----
        - The normal vector is calculated as the difference between the intersection
        point and the sphere's center of curvature, which is located at
        `(0, 0, 1/curvature)`.
        - The resulting normal vector is normalized to have a unit length.
        - This method is marked `noexcept` and `nogil`, making it safe for use in
        performance-critical, multi-threaded Cython code.
        """
        normal = intersection_point - Vector3d(0, 0, 1./self.curvature)
        normal.normalize()

    cdef inline double topo_cy(self, double x, double y) noexcept nogil:
        """
        Calculate the z value for given x and y coordinates on the spherical surface.

        This method computes the topography (z value) of the spherical surface
        corresponding to the provided x and y coordinate. It is primarily used
        for plotting the surface.

        """
        cdef double z
        z = sqrt((1./self.curvature)**2 - x**2 - y**2) - abs(1./self.curvature)

        if self.curvature > 0:
            z = -z
        return z

    def __repr__(self):
        """
        Return a string representation of the Spherical optical surface.

        This method provides a detailed string representation of the `Spherical`
        object, including its shape, reflectivity, curvature, and other relevant
        attributes inherited from the `Surface` class.

        Returns
        -------
        str
            A string that represents the current state of the `Spherical` object.
        """
        return (
            f"Spherical(shape={repr(self.shape)}, "
            f"reflectivity={self.reflectivity}, "
            f"curvature={self.curvature}, "
            f"id={self.id})"
        )
