# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco <AUTHOR>
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

from pyoptools.misc.poly_2d.poly_2d cimport Poly2D
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface

from pyoptools.misc.cmisc.eigen cimport Vector3d, VectorXd, convert_vectorXd_to_list, \
     assign_nan_to_vector3d

from libc.math cimport sqrt, INFINITY

cimport cython

# from ray_trace.surface.taylor_poly import eval_poly,  Poly_DyDx

cdef class Aspherical(Surface):
    """**Class that defines a high order aspherical surface**

    An aspherical surface is defined as::

      Z=(Ax*x**2+Ay*y**2)/(1+sqrt(1-(1+Kx)*Ax**2*x**2-(1+Ky)*Ay**2*y**2))+ Poly2D()

    The Poly2D is defined by a array in the same way as it is defined in the
    TaylorPoly Class

    Example
        >>> cs=Aspherical(shape=Rectangle(size=(5,5)),Ax=.5,Ay=.3,Kx=.1, Ky=.1\
                            poly = Poly2D((0,1,1)))
    """

    cdef public double Ax, Ay, Kx, Ky
    # TODO check the correct type
    cdef public Poly2D poly, DX, DY
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
            self.poly = <Poly2D>poly

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
        """**Return the vector normal to the surface**

        This method returns the vector normal to the asphere at a point
        ``int_p=(x,y,z)``.

        Note: It uses ``x`` and ``y`` to calculate the ``z`` value and the normal.
        """
        cdef double Ax, Ay, Kx, Ky, x, y, _z, dxA, dyA, dxP, dyP

        Ax = self.Ax
        Ay = self.Ay
        Kx = self.Kx
        Ky = self.Ky

        x = intersection_point(0)
        y = intersection_point(1)
        # x, y, _z = intersection_point(0)

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
    cdef double __df1(self, double t, Ray ray) noexcept nogil:
        """
        Numerically calculate the derivative of __f1 with respect to t
        using the central difference method.
        """
        cdef double h = 1e-6  # Step size for numerical differentiation
        cdef double f_plus = self.__f1(t + h, ray)
        cdef double f_minus = self.__f1(t - h, ray)
        return (f_plus - f_minus) / (2 * h)

    @cython.cdivision(True)
    cdef bint find_t_range(self, Ray incident_ray, double& t_min_out, double& t_max_out) noexcept nogil:
        cdef double t_min = 0
        cdef double t_max = INFINITY
        cdef double Ox = incident_ray._origin(0)
        cdef double Oy = incident_ray._origin(1)
        cdef double Oz = incident_ray._origin(2)
        cdef double Dx = incident_ray._direction(0)
        cdef double Dy = incident_ray._direction(1)
        cdef double Dz = incident_ray._direction(2)
        cdef double t1, t2

        # For X dimension
        if Dx != 0:
            t1 = (self.xmin - Ox) / Dx
            t2 = (self.xmax - Ox) / Dx
            if t1 > t2:
                t1, t2 = t2, t1
            t_min = max(t_min, t1)
            t_max = min(t_max, t2)
        elif Ox < self.xmin or Ox > self.xmax:
            return False

        # For Y dimension
        if Dy != 0:
            t1 = (self.ymin - Oy) / Dy
            t2 = (self.ymax - Oy) / Dy
            if t1 > t2:
                t1, t2 = t2, t1
            t_min = max(t_min, t1)
            t_max = min(t_max, t2)
        elif Oy < self.ymin or Oy > self.ymax:
            return False

        # For Z dimension
        if Dz != 0:
            t1 = (self.zmin - Oz) / Dz
            t2 = (self.zmax - Oz) / Dz
            if t1 > t2:
                t1, t2 = t2, t1
            t_min = max(t_min, t1)
            t_max = min(t_max, t2)
        elif Oz < self.zmin or Oz > self.zmax:
            return False

        # Check if we have a valid range
        if t_max <= t_min:
            return False


        t_min_out = t_min
        t_max_out = t_max
        return True


    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:
        cdef double t = 0.0
        cdef double t_new = 0.0
        cdef double Z, f, df
        cdef int max_iterations = 100
        cdef double epsilon = 1e-6
        cdef Vector3d origin = incident_ray._origin
        cdef Vector3d direction = incident_ray._direction
        
        cdef double fval_0
        cdef double fval_1 = self.__f1(t, incident_ray)

        while t<200:
            t=t+0.5
            fval_0 = self.__f1(t, incident_ray)

            if fval_0*fval_1<0:
                break
            fval_1 = fval_0

        
        cdef double ta = t-0.5
        cdef double tb = t+0.5
        cdef double fa = self.__f1(ta, incident_ray)
        cdef double fb = self.__f1(tb, incident_ray)

        for i in range(max_iterations):
            t = ta - fa*(tb-ta)/(fb-fa)
            f = self.__f1(t, incident_ray)
            if abs(f) < epsilon:
                break
            if fa*f<0:
                tb = t
                fb = f
            else:
                ta = t
                fa = f

        if i == max_iterations - 1:
            assign_nan_to_vector3d(intersection_point)
        else:
            intersection_point = origin + direction * t


    def _repr_(self):
        '''Return an string with the representation of an aspherical surface.
        '''
        return "Aspherical(shape="+str(self.shape)+",reflectivity=" +\
            str(self.reflectivity)+",Kx="+str(self.Kx) +\
            ",Ky="+str(self.Ky)+",Ax="+str(self.Ax)+",Ay=" +\
            str(self.Ay)+",poly="+str(self.poly)+")"
