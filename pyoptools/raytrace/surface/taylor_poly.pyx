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
# Description:    Polynomical surface definition module
# ------------------------------------------------------------------------------
#
'''Module that defines a polynomical optical surface, using taylor polynomials
'''

cdef extern from "math.h":
    double sqrt(double)

from numpy import array, dot, inf

from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray
from scipy.optimize import fsolve, brentq
from pyoptools.misc.function_2d.poly_2d.poly_2d cimport Poly2D

# from ray_trace.surface.taylor_poly import eval_poly,  Poly_DyDx
# from pyoptools.misc.poly_2d.poly_2d cimport *

cdef class TaylorPoly(Surface):
    """**Class that defines a high order polynomical surface**

    An aspherical surface is defined as::

      Z= poly2d(x,y)

    Example
        >>> cs=TaylorPoly(shape=Rectangle(size=(5,5)), poly =poly2d((0,1,1)))
    """

    cdef public Poly2D poly
    cdef public double zmax, zmin

    def __init__(self, poly=None, *args, **kwargs):
        Surface.__init__(self, *args, **kwargs)

        if poly is None:
            self.poly = Poly2D([0])
        else:
            self.poly = <Poly2D>poly

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
        self.addkey("poly")
        self.addkey("zmax")
        self.addkey("zmin")

    cdef double topo_cy(self, double x, double y) noexcept nogil:

        """**Returns the Z value for a given X and Y**

        This method returns the topography of the TaylorPoly surface to be
        used to plot the surface.
        """

        return (<Poly2D>self.poly).eval_cy(x, y)

    cpdef normal(self, int_p):
        """**Return the vector normal to the surface**

        This method returns the vector normal to the asphere at a point
        ``int_p=(x,y,z)``.

        Note: It uses ``x`` and ``y`` to calculate the ``z`` value and the normal.
        """
        cdef double x, y, _z, dxP, dyP

        x, y, _z= int_p

        if self.poly is not None:
            Dx, Dy=self.poly.dxdy()
            dxP=Dx.peval(x, y)
            dyP=Dy.peval(x, y)
        else:
            dxP=0.
            dyP=0.

        N_=array((dxP, dyP, -1))

        return N_/sqrt(dot(N_, N_))

    cpdef double __f1(self, double t, Ray iray):
        """Funcion auxiliar que calcula la diferencia en Z entre un punto
        X,Y,Z perteneciente a un rayo, y Z=poly(X,Y). Se usa para encontrar
        los X,Y de interseccion.

        TODO: Hay que buscar una solución analitica
        """
        cdef double Ox, Oy, Oz, Dx, Dy, Dz
        # Ox, Oy, Oz = iray.pos

        Ox=iray.cpos[0]
        Oy=iray.cpos[1]
        Oz=iray.cpos[2]

        # Dx, Dy, Dz = iray.dir

        Dx=iray._dir[0]
        Dy=iray._dir[1]
        Dz=iray._dir[2]

        cdef double X, Y, Z
        X=Dx*t+Ox
        Y=Dy*t+Oy
        Z=Dz*t+Oz
        return self.poly.peval(X, Y)-Z

    cpdef _intersection(self, Ray iray):
        '''**Point of intersection between a ray and the asphere**

        This method returns the point of intersection  between the surface
        and the ray. This intersection point is calculated in the coordinate
        system of the surface.

        It uses an iterative process to calculate the intersection point.


           iray -- incident ray

        iray must be in the coordinate system of the surface
        '''

        # z=pz+t*dz t=(z-pz)/dz
        # Find the limits for t

        cdef double ta, tb, t, fa, fb, tm, tta, ttb, dt

        ta=(self.zmax-iray.origin[2])/iray.direction[2]
        tb=(self.zmin-iray.origin[2])/iray.direction[2]

        if self.poly is not None:
            fa=self.__f1(ta, iray)
            fb=self.__f1(tb, iray)
            if (fa<0 and fb>0) or (fa>0 and fb<0):
                t=brentq(self.__f1, ta, tb, (iray,), maxiter=1000)
            else:  # there are more than 1 intersection points we are assuming 2
                # tm=fsolve(self.__f1, 0,(iray,),warning=False)
                # In new scipy version the warning kw is not supported
                tm=fsolve(self.__f1, 0, (iray,))

                if (tm<ta and tm<tb)or(tm>ta and tm>tb):
                    t=inf
                else:
                    dt=tb-ta
                    tta=tm-0.2*dt
                    ttb=tm+0.2*dt
                    t=brentq(self.__f1, tta, ttb, (iray,), maxiter=1000)

        ret_val= iray.origin+t*iray.direction

        return ret_val

    def _repr_(self):
        '''Return an string with the representation of an aspherical surface.
        '''
        return "TaylorPolt(shape="+str(self.shape)+",reflectivity="+ \
            str(self.reflectivity)+",poly="+str(self.poly)+")"
