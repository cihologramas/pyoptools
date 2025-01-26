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
# Description:     Tube or hollow cylinder surface definition module
# Symbols Defined: Cylinder
# ------------------------------------------------------------------------------

from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface

from libc.math cimport sqrt, INFINITY, sin, cos, M_PI

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_nan_to_vector3d, \
                                        assign_doubles_to_vector3d


cdef class Cylinder(Surface):
    """
    Class to define cylindrical surfaces.

    The `Cylinder` class represents a tube or hollow cylindrical surface in an
    optical system. The cylinder is defined by its radius and length, with the
    center located at the origin of the coordinate system (0, 0, 0) and its
    length aligned along the Z-axis.

    Parameters
    ----------
    radius : double, optional
        The radius of the cylinder. The default value is `10.0`.
    length : double, optional
        The length of the cylinder. The default value is `10.0`.
    *args : tuple, optional
        Additional positional arguments passed to the `Surface` superclass.
    **kwargs : dict, optional
        Additional keyword arguments passed to the `Surface` superclass.

    Attributes
    ----------
    radius : double
        The radius of the cylinder.
    length : double
        The length of the cylinder.

    Examples
    --------
    Creating a cylinder surface with a specific radius and length::

        >>> cs = Cylinder(radius=5.0, length=10.0)

    Notes
    -----
    - The cylinder is centered at the origin of the coordinate system, with
      its length extending along the Z-axis.
    - This class inherits additional functionality and attributes from the `Surface`
      class, which can be used to further define the optical properties of the cylinder.
    - Refer to the `Surface` class documentation for more details on other options
      and attributes that can be used in conjunction with this class.
    """
    cdef public double radius, length

    def __init__(self, radius=10., length=10., *args, **kwargs):
        Surface.__init__(self, *args, **kwargs)
        self.radius = radius
        self.length = length

        # Add attributes to the key list
        self.addkey("radius")
        self.addkey("length")

    # cpdef np.ndarray intersection(self, Ray A):
    cdef void _calculate_intersection(self, Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:
        """
        Calculate the intersection point between a ray and the cylindrical surface.

        This method calculates the point of intersection between the cylindrical
        surface and an incoming ray. The intersection point is determined in
        the coordinate system of the surface.

        Parameters
        ----------
        incident_ray : Ray
            The incoming ray for which the intersection with the cylindrical surface
            is to be calculated. The ray must be in the coordinate system of the
            surface.
        intersection_point : Vector3d&
            A reference to a `Vector3d` object where the calculated intersection
            point will be stored. If there is no valid intersection, this vector
            should be set to indicate no intersection (typically NaN or an invalid
            coordinate).

        Notes
        -----
        - This method is specifically designed for cylindrical surfaces and does not
        rely on external shapes to define the cylinder's boundary. For this reason,
        the `intersection_cy` method is used instead of the more general
        `_intersection` method that is typically overloaded in other surface classes.
        - The calculation assumes that the cylinder is centered at the origin and
        aligned along the Z-axis.
        - This method should be used in contexts where performance is critical, as
        it is defined as `inline` for efficiency.
        """

        # P1 = A.pos     # Punto que pertenece al rayo "Origen" del rayo
        # L1 = A.dir  # Vector paralelo a la linea
        cdef double x1, y1, z1, dx, dy, dz, a, b, c, b2ac, u1, u2, X, X1, X2, \
            Y, Y1, Y2

        x1 = incident_ray._origin(0)
        y1 = incident_ray._origin(1)
        z1 = incident_ray._origin(2)

        # dx, dy, dz = L1

        dx = incident_ray._direction(0)
        dy = incident_ray._direction(1)
        dz = incident_ray._direction(2)

        # This is done by obtaining L from 2 sets of equation
        # x^2+y^2-R^2=0
        # and
        # x=x1+L*dx
        # y=y1+L*dx

        a = dx*dx + dy*dy  # power(dx, 2)+power(dy, 2)
        b = 2*(dx*x1 + dy*y1)

        # power(x1, 2)+power(y1, 2)-power(self.radius, 2)
        c = x1*x1 + y1*y1 - self.radius *self.radius
        b2ac = b*b-4*a*c
        # No hay interseccion
        # o el rayo es tangente y no hay interseccion
        if b2ac <= 0.:
            assign_nan_to_vector3d(intersection_point)
            return

        u1 = (-b+sqrt(b2ac))/(2*a)
        u2 = (-b-sqrt(b2ac))/(2*a)

        # if u1 or u2 = 0, the ray already intersected the cylinder, so u1 and
        # u2 must be >0
        if u1 < 1e-10:
            u1 = INFINITY
        if u2 < 1e-10:
            u2 = INFINITY

        X1 = x1 + u1 * dx
        Y1 = y1 + u1 * dy
        Z1 = z1 + u1 * dz

        X2 = x1 + u2 * dx
        Y2 = y1 + u2 * dy
        Z2 = z1 + u2 * dz

        # Check for the nearest intersection
        # d1 = power(X1-x1, 2)+power(Y1-y1, 2)+power(Z1-z1, 2)
        d1 = (X1-x1)*(X1-x1) + (Y1-y1)*(Y1-y1)+ (Z1-z1)*(Z1-z1)
        # d2 = power(X2-x1, 2)+power(Y2-y1, 2)+power(Z2-z1, 2)
        d2 = (X2-x1)*(X2-x1) + (Y2-y1)*(Y2-y1)+ (Z2-z1)*(Z2-z1)

        if d2 < d1:
            X, Y, Z = X2, Y2, Z2
        else:
            X, Y, Z = X1, Y1, Z1

        if Z < -self.length/2. or Z > self.length/2.:
            assign_nan_to_vector3d(intersection_point)
            return
        else:
            assign_doubles_to_vector3d(X, Y, Z, intersection_point)

    cpdef polylist(self):
        """
        Generate a list of points and polygons representing the cylindrical surface.

        This method generates a list of points and polygons that can be used to
        visualize the cylindrical surface. It is specifically designed for closed
        surfaces, which is why it is overloaded in this class.

        Returns
        -------
        tuple of lists
            - `points`: A list of tuples, where each tuple represents a point
            `(X, Y, Z)` in 3D space that forms part of the cylindrical surface.
            - `polys`: A list of tuples, where each tuple contains three indices
            corresponding to points in the `points` list. These indices define
            the triangular polygons that make up the surface of the cylinder.

        Notes
        -----
        - The method discretizes the cylindrical surface by sampling it at 40
        evenly spaced angles around the Z-axis. For each angle, two points
        are generated: one at the lower end of the cylinder (`Z1`) and one
        at the upper end (`Z2`).
        - The polygons are then created by connecting adjacent points to form
        triangular surfaces, ensuring that the entire cylindrical surface is
        covered.
        - The method is designed to handle the closed nature of the cylindrical
        surface, which is why it is overloaded in this specific class.

        Examples
        --------
        Generating the points and polygons for a cylinder surface::

            points, polys = cylinder.polylist()
        """
        Z1 = -self.length/2.
        Z2 = self.length/2.
        points = []
        polys = []
        for theta in range(40):
            th = theta*2.*M_PI/40.
            X = self.radius*cos(th)
            Y = self.radius*sin(th)
            points.append((X, Y, Z1))
            points.append((X, Y, Z2))

        for i in range(80-2):
            polys.append((i, i+1, i+2))

        polys.append((78, 79, 0))
        polys.append((79, 0, 1))
        return points, polys

    cdef void _calculate_normal(self,
                                Vector3d& intersection_point,
                                Vector3d& normal) noexcept nogil:
        """
        Calculate the normal vector at a given intersection point on the cylindrical
        surface.

        This method computes the normal vector at a specified intersection point on the
        cylindrical surface. The normal vector is calculated in the coordinate system
        of the surface.

        Parameters
        ----------
        intersection_point : Vector3d&
            A reference to a `Vector3d` object representing the point on the cylindrical
            surface where the normal vector is to be calculated.
        normal : Vector3d&
            A reference to a `Vector3d` object where the calculated normal vector
            will be stored. The Z-component of the normal is set to zero to ensure
            the normal vector lies in the X-Y plane.

        Notes
        -----
        - The normal vector is calculated by copying the `intersection_point` into
        the `normal` vector and then setting the Z-component to zero. This ensures
        that the normal vector is perpendicular to the cylindrical surface and lies
        in the X-Y plane.
        - The resulting normal vector is then normalized to have a unit length.
        - This method is marked `noexcept` and `nogil`, making it suitable for use
        in performance-critical, multi-threaded Cython code.

        """
        normal = intersection_point
        (<double*>(&(normal(2))))[0] = 0
        normal.normalize()

    def __repr__(self):
        """
        Return a string representation of the Cylinder optical surface.

        This method provides a detailed string representation of the `Cylinder`
        object, including its shape, reflectivity, radius, length, and other
        relevant attributes inherited from the `Surface` class.

        Returns
        -------
        str
            A string that represents the current state of the `Cylinder` object.
        """
        return (
            f"Cylinder("
            f"reflectivity={self.reflectivity}, "
            f"radius={self.radius}, "
            f"length={self.length}, "
            f"id={self.id})"
        )
