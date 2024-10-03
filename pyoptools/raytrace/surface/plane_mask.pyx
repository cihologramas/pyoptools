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
# Description:     Plane surface definition module
# Symbols Defined: Plane
# ------------------------------------------------------------------------------

'''Module that defines a reflective plane phase mask surface class
'''
from numpy import array, inf, float64
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.misc.poly_2d.poly_2d cimport Poly2D
from pyoptools.misc.definitions import inf_vect
import cython
cdef extern from "math.h":
    double sqrt(double)
    double copysign(double, double)
    double M_PI


cdef class RPPMask(Surface):
    '''Class to define a reflective plane phase mask.

    Description:

    RPPMask  is a class to define diffractive plane surfaces.
    If reflectivity is 1 the gratting is reflective, if it is 0 the
    gratting is transmissive. If it is between 0 and 1, both transmitted
    and reflected diffracted rays are shown.

    The surface shape is given by the shape attribute

    The phm attribute is a Poly2D instance, that contains the polinomial
    describing the phase modulation of the gratting. The X and Y input
    values of the polynomial are in microns.

    Example:

        >>> g=RPPMask(shape=Rectangular(size=(10,10)), phm=Poly2D([0,2*pi/15.,0]),M=[1])

    This is a 10 mm x 10 mm linear gratting that has a fringe each 15 microns
    '''

    # cdef np.ndarray k
    cdef public Poly2D phx
    cdef public Poly2D phy
    cdef public Poly2D phm
    cdef public list M
    # def __init__(self,k=(0,2*pi/1e-3),M=[1],**traits):

    def __init__(self, Poly2D phm=None, M=None, *args, **kwargs):
        """
        phm represent the phase variations across the surface
        """
        Surface.__init__(self, *args, **kwargs)
        if phm is None:
            self.phm = Poly2D((0, 0, 0, 0, 0, 0))
        else:
            self.phm = phm
        dxdy = self.phm.dxdy()
        self.phx = <Poly2D > dxdy[0]
        self.phy = <Poly2D > dxdy[1]
        if M is None:
            self.M =[1]
        else:
            self.M = M

        # Add attributes to the state list
        self.addkey("phm")
        self.addkey("phx")
        self.addkey("phy")
        self.addkey("M")

    # ~ def __reduce__(self):
        # ~
        # ~ args=(self.phm, self.M, self.reflectivity, self.shape)
        # ~ return(type(self),args,self.__getstate__())
        # ~

    cdef inline double topo_cy(self, double x, double y) noexcept nogil:
        return 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef _intersection(self, Ray A):
        """Returns the intersection point between a ray and an the XY plane

        """
        # N_=array([0.,0.,1.])

        # Punto que pertenece al rayo "Origen" del rayo
        # cdef np.ndarray[np.float64_t, ndim=1]
        P1 = A.origin
        # Vector paralelo a la linea
        # cdef np.ndarray[np.float64_t, ndim=1]
        L1 = A.direction

        # if dot(N_,L1) ==0 : return inf_vect
        if L1[2] == 0:
            return inf_vect

        # print N_,P1,L1
        # print dot(N_,-P1),dot(N_,L1)
        # u=dot(N_,-P1)/dot(N_,L1)
        cdef double u = -P1[2]/L1[2]
        # Si u es muy grande, no hay intersección

        retval = P1+u*L1

        return retval

    cpdef normal(self, ri):
        """Method that returns the normal to the surface
        """
        N_ = array((0., 0., 1.)).astype(float64)
        return (N_)

    cpdef list propagate(self, Ray ri, double ni, double nr):
        """Method that calculates the propagation of a ray through a diffraction
        gratting.

        This method uses all the units in millimeters
        """
        cdef double l, rx, ry, rz, k, d, kx, ky, ox, oy, oz, oz2, M
        PI, _P = self.int_nor(ri)
        # Express wavelength in millimeters so all the method works in millimeters
        l = ri.wavelength*1e-3
        rx, ry, rz = ri.direction
        K = array((self.phx.peval(PI[0], PI[1]), self.phy.peval(PI[0], PI[1])))
        k = sqrt(K[0]**2+K[1]**2)
        if k != 0:
            d = 2.*M_PI/k
            kx, ky = K/k
        else:
            kx = 0
            ky = 0
            d = inf

        ret = []

        for M in self.M:
            ox = rx+M*l/d * kx
            oy = ry+M*l/d * ky
            oz2 = rz**2 - 2*M*l/d*(rx*kx+ry*ky)-(M*l/d)**2

            if oz2 < 0:
                print("warning: eliminating physically impossible ray")
                ret.append(Ray(origin=PI, direction=ri.direction,
                               intensity=0,
                               wavelength=ri.wavelength, n=ni, label=ri.label,
                               orig_surf=self.id))
            else:
                oz = copysign(sqrt(oz2), rz)
                # Check for transmitted or and reflected orders. Here intensities
                # have no meaning
                if self.reflectivity != 1:
                    ret.append(Ray(origin=PI, direction=(ox, oy, oz),
                                   intensity=ri.intensity/len(self.M),
                                   wavelength=ri.wavelength, n=ni, label=ri.label,
                                   orig_surf=self.id))
                if self.reflectivity != 0:
                    ret.append(Ray(origin=PI, direction=(ox, oy, -oz),
                                   intensity=ri.intensity/len(self.M),
                                   wavelength=ri.wavelength, n=ni,
                                   label=ri.label, orig_surf=self.id))

        return ret
