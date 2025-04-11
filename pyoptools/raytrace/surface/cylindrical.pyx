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
# Description:     Cylindrical surface definition module
# Symbols Defined: CylindricalSurface
# ------------------------------------------------------------------------------


from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_nan_to_vector3d, \
                                        assign_doubles_to_vector3d

from libc.math cimport sqrt, abs

cdef class Cylindrical(Surface):
    """
    Class to define cylindrical optical surfaces.

    The `Cylindrical` class represents a cylindrical optical surface. The surface
    is defined by specifying the aperture shape and the radius of curvature
    of the cylinder.

    The vertex of the cylindrical surface is located at the origin of coordinates
    (0, 0, 0).

    Parameters
    ----------
    curvature : double, optional
        The curvature of the cylindrical surface. A positive curvature indicates
        a convex surface, while a negative curvature indicates a concave surface.
        The default is `0.0`.
    *args : tuple, optional
        Additional positional arguments passed to the `Surface` superclass.
    **kwargs : dict, optional
        Additional keyword arguments passed to the `Surface` superclass.

    Attributes
    ----------
    curvature : double
        The curvature of the cylindrical surface, as specified during initialization.

    Examples
    --------
    Creating a cylindrical surface with a rectangular aperture and a specific
    curvature::

        >>> cs = Cylindrical(shape=Rectangular(size=(10, 20)), curvature=0.15)

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
        Calculate the intersection point between a ray and the cylindrical surface.

        This method calculates the point of intersection between the cylindrical
        surface and a given incident ray. The intersection point is calculated
        in the coordinate system of the surface. If no valid intersection exists
        (e.g., the ray is outside the element's aperture or the ray does not
        intersect the surface), the intersection point is set to NaN.

        Parameters
        ----------
        incident_ray : Ray
            The incoming ray for which the intersection with the cylindrical
            surface is to be calculated. The ray must be in the coordinate
            system of the surface.
        intersection_point : Vector3d&
            A reference to a `Vector3d` object where the calculated intersection
            point will be stored. If there is no valid intersection, this vector
            is set to NaN.

        Notes
        -----
        - This method is overloaded from the `Surface` superclass and does not check
        for the aperture. It only calculates the geometric intersection point.
        - The method solves the quadratic equation derived from the intersection
        of the ray with the cylindrical surface. Two potential intersection points
        are computed, and the one closest to the surface is selected.
        - If the calculated intersection point lies outside the bounds of the surface
        (e.g., the Z-coordinate exceeds the cylinder's limits), the intersection
        point is set to NaN.
        - This method should not be called directly. Instead, use `Surface.intersection`
        to perform intersection calculations, which includes aperture checks.

        """
        cdef double x1, y1, z1, x2, y2, z2, z3, x21, z21
        cdef double X1, Y1, Z1, X2, Y2, Z2, X, Y, Z
        cdef double a, b, c, b2ac, u1, u2

        x1 = incident_ray._origin(0)
        y1 = incident_ray._origin(1)
        z1 = incident_ray._origin(2)

        x2 = x1+incident_ray._direction(0)
        y2 = y1+incident_ray._direction(1)
        z2 = z1+incident_ray._direction(2)

        x21 = incident_ray._direction(0)
        z21 = incident_ray._direction(2)

        z3 = 1./self.curvature

        a = (x2-x1)*(x21)+(z21)*(z21)
        b = 2*((x21)*(x1) + (z21)*(z1 - z3))
        c = x1*x1+z1*z1-2*(z3*z1)
        b2ac = b*b-4*a*c

        if a == 0 or b2ac <= 0.:
            assign_nan_to_vector3d(intersection_point)
        else:
            u1 = (-b+sqrt(b2ac))/(2*a)
            u2 = (-b-sqrt(b2ac))/(2*a)

            X1 = x1 + u1 * (x2 - x1)
            Y1 = y1 + u1 * (y2 - y1)
            Z1 = z1 + u1 * (z2 - z1)

            X2 = x1 + u2 * (x2 - x1)
            Y2 = y1 + u2 * (y2 - y1)
            Z2 = z1 + u2 * (z2 - z1)

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
        Calculate the normal vector to the cylindrical surface at a given point.

        This method calculates the normal vector to the cylindrical surface at the
        specified intersection point. The normal is computed in the coordinate
        system of the surface.

        Parameters
        ----------
        intersection_point : Vector3d&
            A reference to a `Vector3d` object representing the point on the
            cylindrical surface where the normal vector is to be calculated.
        normal : Vector3d&
            A reference to a `Vector3d` object where the calculated normal vector
            will be stored.

        Notes
        -----
        - The normal vector is calculated as the difference between the intersection
        point and the center of curvature of the cylinder, which is located at
        `(0, 0, 1/curvature)`.
        - The Y component of the normal vector is set to zero to ensure the normal
        lies in the XZ plane, which is consistent with the geometry of a
        cylindrical surface.
        - The resulting normal vector is then normalized to have a unit length.
        - This method is marked `noexcept` and `nogil`, making it safe for use in
        performance-critical, multi-threaded Cython code.
        """
        normal = intersection_point - Vector3d(0, 0, 1./self.curvature)
        (<double*>(&(normal(1))))[0] = 0
        normal.normalize()

    cdef inline double topo_cy(self, double x, double y) noexcept nogil:
        """
        Calculate the topography (Z values) of the cylindrical surface for
        given X and Y coordinates.

        This method computes the topography of the cylindrical surface, represented
        as a z value, for the provided x and x coordinate. The z values correspond
        to the height of the cylindrical surface at each (x, y) coordinate pair.
        This method is overloaded from the `Surface` superclass.
        """

        cdef double z
        z = sqrt((1./self.curvature)**2 - x**2) - abs(1./self.curvature)
        if self.curvature > 0:
            z = -z
        return z

    def __repr__(self):
        """
        Return a string representation of the Cylindrical optical surface.

        This method provides a detailed string representation of the `Cylindrical`
        object, including its shape, reflectivity, curvature, and other relevant
        attributes inherited from the `Surface` class.

        Returns
        -------
        str
            A string that represents the current state of the `Cylindrical` object.
        """
        return (
            f"Cylindrical(shape={repr(self.shape)}, "
            f"reflectivity={self.reflectivity}, "
            f"curvature={self.curvature}, "
            f"id={self.id})"
        )
