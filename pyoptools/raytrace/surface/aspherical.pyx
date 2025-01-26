# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:         Ricardo Amézquita Orozco
# Description:    Aspherical surface definition module
# ------------------------------------------------------------------------------
#
'''Module that defines support for Aspherical optical surfaces
'''

from pyoptools.misc.function_2d.function_2d cimport Function2D
from pyoptools.misc.function_2d.poly_2d.poly_2d cimport Poly2D
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface

from pyoptools.misc.cmisc.eigen cimport Vector3d, VectorXd, convert_vectorXd_to_list, \
     assign_nan_to_vector3d

from libc.math cimport sqrt, INFINITY
import warnings

cimport cython

cdef class Aspherical(Surface):
    """
    Class representing a high-order aspherical optical surface.

    An aspherical surface is defined by the equation:

    Z = (Ax*x**2 + Ay*y**2) /
        (1 + sqrt(1 - (1+Kx)*Ax**2*x**2 - (1+Ky)*Ay**2*y**2))
        + Poly2D(x, y)

    The aspherical surface combines a base aspheric shape with an additional
    polynomial deformation term defined by Function2D.

    Parameters
    ----------
    shape : Shape object
        The shape of the aspherical surface.
    Ax : float, optional
        X-axis curvature coefficient. Default is 0.
    Ay : float, optional
        Y-axis curvature coefficient. Default is 0.
    Kx : float, optional
        X-axis conic constant. Default is 0.
    Ky : float, optional
        Y-axis conic constant. Default is 0.
    poly : Poly2D, optional
        Polynomial surface deformation term. If None, a zero-order polynomial
        is used (no deformation).
    *args, **kwargs
        Additional arguments passed to the parent Surface class.

    Attributes
    ----------
    xmin, xmax : float
        Bounds of the surface in the x direction.
    ymin, ymax : float
        Bounds of the surface in the y direction.
    zmin, zmax : float
        Bounds of the surface in the z direction.
    DX, DY : Poly2D
        Partial derivatives of the polynomial deformation term.
    Examples
    --------
    >>> from pyoptools.raytrace.surface import Aspherical
    >>> from pyoptools.raytrace.shape import Rectangle
    >>> from pyoptools.misc.function_2d.poly_2d.poly_2d import Poly2D
    >>> surface = Aspherical(
    ...     shape=Rectangle(size=(5,5)),
    ...     Ax=0.5,
    ...     Ay=0.3,
    ...     Kx=0.1,
    ...     Ky=0.1,
    ...     poly=Poly2D((0,1,1))
    ... )

    Notes
    -----
    The surface shape is determined by the combination of the base aspherical
    equation and the polynomial deformation term. The polynomial term allows
    for modeling complex surface irregularities or corrections.
    """

    cdef public double Ax, Ay, Kx, Ky
    cdef public Function2D poly, DX, DY
    cdef public double xmin, xmax, ymin, ymax, zmax, zmin

    def __init__(self, Ax=0., Ay=0., Kx=0., Ky=0., poly=None, *args, **kwargs):
        Surface.__init__(self, *args, **kwargs)
        self.Ax = Ax
        self.Ay = Ay
        self.Kx = Kx
        self.Ky = Ky

        cdef VectorXd zero_c = VectorXd(1)

        if poly is None:
            self.poly = Poly2D(convert_vectorXd_to_list(zero_c))
        else:
            self.poly = <Function2D>poly

        self.DX, self.DY = self.poly.dxdy()

        # Find the X,Y,Z volume where the Aspherical Surface is defined
        # The X Y limits are provided by the shape. The Z limit will be aproxi
        # by sampling the surface
        cdef double xmax, xmin, ymax, ymin, zmax, zmin

        xmin, xmax, ymin, ymax = self.shape.limits()

        cdef int ndat = 200
        cdef double dx = (xmax -xmin)/ndat
        cdef double dy = (ymax -ymin)/ndat
        cdef int ix, iy
        cdef double x, y, z

        # Sample z in the middle of the interval
        zmax = self.topo_cy((xmax+xmin)/2, (ymax+ymin)/2)
        zmin = zmax
        for ix in range(ndat+1):
            x=xmin + ix* dx
            for iy in range(ndat+1):
                y=ymin + iy* dy
                z = self.topo_cy(x, y)
                zmax = zmax if zmax > z else z
                zmin = zmin if zmin < z else z

        cdef double dz = zmax - zmin

        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        # Increase a little the bound box just for safety
        self.zmax = zmax+0.01*dz
        self.zmin = zmin-0.01*dz

        # Add attributes to the state list
        self.addkey("Ax")
        self.addkey("Ay")
        self.addkey("Kx")
        self.addkey("Ky")
        self.addkey("poly")
        self.addkey("zmax")
        self.addkey("zmin")

    @cython.cdivision(True)
    cdef double topo_cy(self, double x, double y) noexcept nogil:
        """
        Calculate the topography of the aspherical surface.

        This method computes the z-coordinate of a point on the aspherical surface
        given its x and y coordinates.

        Parameters
        ----------
        x : double
            The x-coordinate of the point.
        y : double
            The y-coordinate of the point.

        Returns
        -------
        double
            The z-coordinate of the point on the aspherical surface.

        Notes
        -----
        The aspherical surface is defined by the equation:
        Z = (Ax*x**2 + Ay*y**2) / (1 + sqrt(1 - (1+Kx)*Ax**2*x**2 - (1+Ky)*Ay**2*y**2))
            + Poly2D()

        Where Ax, Ay are the curvature coefficients, Kx, Ky are the conic constants,
        and Poly2D() is an additional polynomial term for surface deformation.
        """
        cdef double Ax = self.Ax
        cdef double Ay = self.Ay
        cdef double Kx = self.Kx
        cdef double Ky = self.Ky

        cdef double Z0 = (Ax*x**2+Ay*y**2) / \
            (1+sqrt(1-(1+Kx)*Ax**2*x**2+(-(1+Ky))*Ay**2*y**2))

        cdef double Z1 = self.poly.eval_cy(x, y)

        return Z0+Z1

    cdef void _calculate_normal(self,
                                Vector3d& intersection_point,
                                Vector3d& normal) noexcept nogil:
        """
        Calculate the normal vector to the aspherical surface at a given point.

        This method computes the vector normal to the aspherical surface at the
        intersection point.

        Parameters
        ----------
        intersection_point : Vector3d
            The point on the surface where the normal is to be calculated.
        normal : Vector3d
            The output vector where the calculated normal will be stored.

        Returns
        -------
        None
            The method modifies the `normal` vector in-place.

        Notes
        -----
        The method uses the x and y coordinates of the intersection point to
        calculate the z value and the normal vector. The normal vector is
        normalized before returning.

        The calculation takes into account both the basic aspherical equation
        and the additional polynomial deformation term.
        """
        cdef double Ax, Ay, Kx, Ky, x, y, _z, dxA, dyA, dxP, dyP

        Ax = self.Ax
        Ay = self.Ay
        Kx = self.Kx
        Ky = self.Ky

        x = intersection_point(0)
        y = intersection_point(1)

        dxA = (2*Ax*x)/(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1) + \
              (Ax**2*(Kx+1)*x*(Ay*y**2+Ax*x**2)) / \
              (sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1) *
               (sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)**2)

        dyA = (2*Ay*y)/(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)- \
              (Ay**2*(-Ky-1)*y*(Ay*y**2+Ax*x**2)) / \
              (sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1) *
               (sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)**2)

        dxP = self.DX.eval_cy(x, y)
        dyP = self.DY.eval_cy(x, y)

        normal = Vector3d(dxA+dxP, dyA+dyP, -1)
        normal.normalize()

    cdef double __f1(self, double t, Ray iray) noexcept nogil:
        """
        Auxiliary function to find the intersection point.

        This function calculates the difference between the z-coordinate of a point
        on the ray and the corresponding z-coordinate on the aspherical surface.
        It is used in the root-finding algorithm to determine the intersection point.

        Parameters
        ----------
        t : double
            The parameter along the ray.
        iray : Ray
            The incident ray.

        Returns
        -------
        double
            The difference between the z-coordinate on the ray and the z-coordinate
            on the aspherical surface at the given (x, y) point.
        """
        cdef double Ox, Oy, Oz, Dx, Dy, Dz
        cdef double X, Y, Z

        # Ox, Oy, Oz = iray.pos
        Ox = iray._origin(0)
        Oy = iray._origin(1)
        Oz = iray._origin(2)

        # Dx, Dy, Dz = iray.dir
        Dx = iray._direction(0)
        Dy = iray._direction(1)
        Dz = iray._direction(2)

        X = Dx * t + Ox
        Y = Dy * t + Oy
        Z = Dz * t + Oz
        return self.topo_cy(X, Y) - Z

    @cython.cdivision(True)
    cdef bint find_t_range(self, Ray incident_ray, double& t_min_out,
                           double& t_max_out) noexcept nogil:
        """
        Find the range of t values for which the ray intersects the bounding box of the
        aspherical surface.

        This method calculates the minimum and maximum t values for which the given ray
        intersects the bounding box of the aspherical surface. It checks intersections
        with the bounding planes in X, Y, and Z dimensions.

        Parameters
        ----------
        incident_ray : Ray
            The incident ray for which to calculate the intersection range.
        t_min_out : double&
            Output parameter for the minimum t value.
        t_max_out : double&
            Output parameter for the maximum t value.

        Returns
        -------
        bint
            True if a valid intersection range is found, False otherwise.

        Notes
        -----
        This method uses the ray equation P(t) = O + tD, where O is the ray origin
        and D is the ray direction, to calculate intersections with the bounding box.
        """
        cdef double t_min = -INFINITY
        cdef double t_max = INFINITY
        cdef double Ox = incident_ray._origin(0)
        cdef double Oy = incident_ray._origin(1)
        cdef double Oz = incident_ray._origin(2)
        cdef double Dx = incident_ray._direction(0)
        cdef double Dy = incident_ray._direction(1)
        cdef double Dz = incident_ray._direction(2)
        cdef double t1, t2

        cdef double epsilon = 1e-10

        # For X dimension
        if abs(Dx) > epsilon:  # Dx != 0
            t1 = (self.xmin - Ox) / Dx
            t2 = (self.xmax - Ox) / Dx
            if t1 > t2:
                t1, t2 = t2, t1
            t_min = max(t_min, t1)
            t_max = min(t_max, t2)
        elif Ox < self.xmin or Ox > self.xmax:
            return False

        # For Y dimension
        if abs(Dy) > epsilon:  # Dy != 0:
            t1 = (self.ymin - Oy) / Dy
            t2 = (self.ymax - Oy) / Dy
            if t1 > t2:
                t1, t2 = t2, t1
            t_min = max(t_min, t1)
            t_max = min(t_max, t2)
        elif Oy < self.ymin or Oy > self.ymax:
            return False

        # For Z dimension
        if abs(Dz) > epsilon:  # Dz != 0:
            t1 = (self.zmin - Oz) / Dz
            t2 = (self.zmax - Oz) / Dz
            if t1 > t2:
                t1, t2 = t2, t1
            t_min = max(t_min, t1)
            t_max = min(t_max, t2)
        elif Oz < self.zmin or Oz > self.zmax:
            return False

        # Check if we have a valid range
        if t_max < 0:
            return False

        t_min_out = t_min
        t_max_out = t_max
        return True

    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:
        """
        Calculate the intersection point between a ray and the aspherical surface.

        This method uses a combination of bracketing and the secant method to find
        the intersection point of an incident ray with the aspherical surface.

        Parameters
        ----------
        incident_ray : Ray
            The incident ray to intersect with the surface.
        intersection_point : Vector3d
            A vector to store the calculated intersection point.

        Returns
        -------
        None
            The method modifies the intersection_point vector in-place.

        Notes
        -----
        The method first finds a range of t values where the intersection might occur,
        then uses a bracketing method to narrow down the range, and finally applies
        the secant method for precise intersection calculation.

         Limitations
        ----------
        - The algorithm assumes a single intersection point. In cases where there are
          multiple intersections (odd or even number), the result may be incorrect or
          the ray may be ignored:
          * Even number of intersections: the ray will be ignored (set to NaN)
          * Odd number of intersections: the algorithm will proceed but may return
            an incorrect intersection point, as it cannot distinguish between multiple
            intersection points
        - The algorithm may fail to converge after the maximum number of iterations
          (default: 100). In such cases, the intersection point will be set to NaN
          and the ray will be ignored.
        - The method relies on a convergence tolerance (epsilon) of 1e-10. Points
          requiring higher precision may not be accurately calculated.
        - The intersection must occur within the bounds determined by the bounding box.
          Intersections outside these bounds will not be detected.

        If any of these limitations are encountered, the method will either return
        potentially incorrect results or set the intersection_point to NaN values
        and issue a warning.
        """
        # Initialize variables for the intersection calculation
        cdef double t
        cdef double f
        cdef int max_iterations = 100  # Maximum number of iterations for convergence
        cdef double epsilon = 1e-10    # Tolerance for convergence
        cdef Vector3d origin = incident_ray._origin        # Ray origin point
        cdef Vector3d direction = incident_ray._direction  # Ray direction vector

        # Variables to store the range of possible intersection points
        cdef double t_min, t_max

        # Check if the ray intersects the bounding box
        if not self.find_t_range(incident_ray, t_min, t_max):
            assign_nan_to_vector3d(intersection_point)
            return

        # Initialize bracketing interval endpoints and their function values
        cdef double ta = t_min  # Start of interval
        cdef double tb = t_max  # End of interval
        cdef double fa = self.__f1(ta, incident_ray)  # Function value at start
        cdef double fb = self.__f1(tb, incident_ray)  # Function value at end

        # Check if there's an odd number of intersections
        # If fa*fb > 0, the ray either doesn't intersect or intersects an even
        # number of times
        if fa*fb > 0:
            with gil:
                warnings.warn(
                    "The ray is intersecting the surface an even number of times. "
                    "This is not supported at the moment, and the ray will be ignored.",
                    RuntimeWarning,
                )

            assign_nan_to_vector3d(intersection_point)
            return

        # Secant method iteration loop
        for i in range(max_iterations):
            # Calculate new approximation using secant method
            t = ta - fa*(tb-ta)/(fb-fa)
            f = self.__f1(t, incident_ray)

            # Check if we've reached desired accuracy
            if abs(f) < epsilon:
                break

            # Update bracketing interval
            if fa*f<0:  # If sign change between fa and f
                tb = t
                fb = f
            else:       # If sign change between f and fb
                ta = t
                fa = f

        # If maximum iterations reached without convergence, set to NaN
        if i == max_iterations - 1:
            with gil:
                warnings.warn(
                    "Maximum iterations reached without convergence in aspherical"
                    "surface intersection calculation. The ray will be ignored.",
                    RuntimeWarning,
                )

            assign_nan_to_vector3d(intersection_point)
        else:
            # Calculate the intersection point using the parameter t
            intersection_point = origin + direction * t

    def _repr_(self):
        '''Return an string with the representation of an aspherical surface.
        '''
        return "Aspherical(shape="+str(self.shape)+",reflectivity=" +\
            str(self.reflectivity)+",Kx="+str(self.Kx) +\
            ",Ky="+str(self.Ky)+",Ax="+str(self.Ax)+",Ay=" +\
            str(self.Ay)+",poly="+str(self.poly)+")"
