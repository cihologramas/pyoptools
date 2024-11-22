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


from scipy.optimize import fsolve
from pyoptools.misc.poly_2d.poly_2d cimport Poly2D
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface

from pyoptools.misc.cmisc.eigen cimport Vector3d, VectorXd, convert_vectorXd_to_list, \
     assign_nan_to_vector3d

from libc.math cimport sqrt

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
    cdef public double zmax, zmin

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

    cdef double __df1(self, double t, Ray ray) noexcept nogil:
        """
        Numerically calculate the derivative of __f1 with respect to t
        using the central difference method.
        """
        cdef double h = 1e-6  # Step size for numerical differentiation
        cdef double f_plus = self.__f1(t + h, ray)
        cdef double f_minus = self.__f1(t - h, ray)
        return (f_plus - f_minus) / (2 * h)

    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:
        """
        Calculate the point of intersection between a ray and the asphere.

        This method returns the point of intersection between the surface
        and the ray. The intersection point is calculated in the coordinate
        system of the surface using an iterative process.

        Parameters
        ----------
        incident_ray : Ray
            The incident ray in the coordinate system of the surface.

        Returns
        -------
        Vector3d
            The point of intersection between the ray and the asphere.

        Notes
        -----
        The incident ray must be in the coordinate system of the surface.
        """
        cdef double t, Z
        with gil:
            t0 = 0
            t_sol, _info, ier, _mesg = fsolve(self.__f1, t0, (incident_ray) ,
                                              full_output=True)

            if ier == 1:
                t = t_sol[0]
                Z=incident_ray._origin(2)+t*incident_ray._direction(2)

                if self.zmax > Z > self.zmin:
                    intersection_point = incident_ray._origin+incident_ray._direction*t
                else:
                    assign_nan_to_vector3d(intersection_point)
            else:
                assign_nan_to_vector3d(intersection_point)

    def _repr_(self):
        '''Return an string with the representation of an aspherical surface.
        '''
        return "Aspherical(shape="+str(self.shape)+",reflectivity=" +\
            str(self.reflectivity)+",Kx="+str(self.Kx) +\
            ",Ky="+str(self.Ky)+",Ax="+str(self.Ax)+",Ay=" +\
            str(self.Ay)+",poly="+str(self.poly)+")"
