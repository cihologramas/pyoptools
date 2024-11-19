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


from scipy.optimize import fsolve, brentq
from pyoptools.misc.poly_2d.poly_2d cimport Poly2D
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface
from numpy import inf

from pyoptools.misc.cmisc.eigen cimport Vector3d, VectorXd, convert_vectorXd_to_list

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
    cdef public Poly2D poly
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
        with gil:
            dxA = (2*Ax*x)/(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1) + \
                (Ax**2*(Kx+1)*x*(Ay*y**2+Ax*x**2)) / \
                (sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1) *
                 (sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)**2)

            dyA = (2*Ay*y)/(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)- \
                (Ay**2*(-Ky-1)*y*(Ay*y**2+Ax*x**2)) / \
                (sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1) *
                 (sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)**2)

            Dx, Dy = self.poly.dxdy()
            dxP = Dx.peval(x, y)
            dyP = Dy.peval(x, y)

        normal = Vector3d(dxA+dxP, dyA+dyP, -1)
        normal.normalize()

    cpdef double __f1(self, double t, Ray iray):
        cdef double Ax, Ay, Kx, Ky, Ox, Oy, Oz, Dx, Dy, Dz
        Ax = self.Ax
        Ay = self.Ay
        Kx = self.Kx
        Ky = self.Ky

        # Ox, Oy, Oz = iray.pos

        Ox = iray.cpos[0]
        Oy = iray.cpos[1]
        Oz = iray.cpos[2]

        # Dx, Dy, Dz = iray.dir

        Dx = iray._dir[0]
        Dy = iray._dir[1]
        Dz = iray._dir[2]

        cdef double X, Y, Z
        X = Dx*t+Ox
        Y = Dy*t+Oy
        Z = Dz*t+Oz
        return (Ay*Y**2+Ax*X**2)/(sqrt(Ay**2*(-Ky-1)*Y**2-Ax**2*(Kx+1)*X**2+1)+1) + \
            self.poly.peval(X, Y) - Z

    cpdef double __f2(self, double t, iray):
        cdef double Ax, Ay, Kx, Ky, Ox, Oy, Oz, Dx, Dy, Dz

        Ax = self.Ax
        Ay = self.Ay
        Kx = self.Kx
        Ky = self.Ky

        # Ox, Oy, Oz = iray.pos
        Ox = iray.cpos[0]
        Oy = iray.cpos[1]
        Oz = iray.cpos[2]

        # Dx, Dy, Dz = iray.dir
        Dx = iray._dir[0]
        Dy = iray._dir[1]
        Dz = iray._dir[2]

        return (Ay*(Dy*t+Oy)**2+Ax*(Dx*t+Ox)**2)/(sqrt(Ay**2*(-Ky-1)*(Dy*t+Oy)**2-
                                                  Ax**2*(Kx+1)*(Dx*t+Ox)**2+1)+1) - \
            Dz*t-Oz

    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:

        '''**Point of intersection between a ray and the asphere**

        This method returns the point of intersection  between the surface
        and the ray. This intersection point is calculated in the coordinate
        system of the surface.

        It uses an iterative process to calculate the intersection point.


           incident_ray -- incident ray

        incident_ray must be in the coordinate system of the surface
        '''

        # z=pz+t*dz t=(z-pz)/dz
        # Find the limits for t
        cdef double ta, tb, t, fa, fb, tm, tta, ttb, dt
        with gil:
            ta = (self.zmax-incident_ray.origin[2])/incident_ray.direction[2]
            tb = (self.zmin-incident_ray.origin[2])/incident_ray.direction[2]

            if self.poly is not None:
                fa = self.__f1(ta, incident_ray)
                fb = self.__f1(tb, incident_ray)
                if (fa < 0 and fb > 0) or (fa > 0 and fb < 0):
                    t = brentq(self.__f1, ta, tb, (incident_ray,), maxiter=1000)
                else:  # there are more than 1 intersection points we are assuming 2
                    # tm=fsolve(self.__f1, 0,(incident_ray,),warning=False)
                    # In new scipy version the warning kw is not supported
                    tm = fsolve(self.__f1, 0, (incident_ray,))

                    if (tm < ta and tm < tb) or (tm > ta and tm > tb):
                        t = inf
                    else:
                        dt = tb-ta
                        tta = tm-0.2*dt
                        ttb = tm+0.2*dt
                        t = brentq(self.__f1, tta, ttb, (incident_ray,), maxiter=1000)

            else:
                fa = self.__f2(ta, incident_ray)
                fb = self.__f2(tb, incident_ray)

                if (fa < 0 and fb > 0) or (fa > 0 and fb < 0):
                    t = brentq(self.__f2, ta, tb, (incident_ray,), maxiter=1000)
                else:  # there are more than 1 intersection points we are assuming 2
                    # In new scipy version the warning kw is not supported
                    # tm=fsolve(self.__f2, 0,(incident_ray,),warning=False)
                    tm = fsolve(self.__f2, 0, (incident_ray,))

                    if (tm < ta and tm < tb) or (tm > ta and tm > tb):
                        t = inf
                    else:
                        dt = tb-ta
                        tta = tm-0.1*dt
                        ttb = tm+0.1*dt
                        t = brentq(self.__f2, tta, ttb, (incident_ray,), maxiter=1000)

        intersection_point = incident_ray._origin+incident_ray._direction*t

    def _repr_(self):
        '''Return an string with the representation of an aspherical surface.
        '''
        return "Aspherical(shape="+str(self.shape)+",reflectivity=" +\
            str(self.reflectivity)+",Kx="+str(self.Kx) +\
            ",Ky="+str(self.Ky)+",Ax="+str(self.Ax)+",Ay=" +\
            str(self.Ay)+",poly="+str(self.poly)+")"
