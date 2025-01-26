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
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.plane cimport Plane
from pyoptools.misc.function_2d.poly_2d.poly_2d cimport Poly2D
from pyoptools.misc.cmisc.eigen cimport Vector3d, Vector2d
import cython

from libc.math cimport sqrt, copysign, M_PI, INFINITY

cdef class RPPMask(Plane):
    """Reflective plane phase mask surface.

    A class to define diffractive plane surfaces that can be either reflective,
    transmissive, or both. The surface's behavior is controlled by its reflectivity:
    1 for purely reflective, 0 for purely transmissive, or a value between 0 and 1
    for both behaviors.

    Parameters
    ----------
    phm : Poly2D, optional
        Polynomial describing the phase modulation of the grating.
        The X and Y input values of the polynomial are in microns.
        Default is a zero polynomial.
    M : list, optional
        List of diffraction orders to consider.
        Default is [1].
    *args, **kwargs
        Additional arguments passed to the parent Plane class.

    Attributes
    ----------
    phm : Poly2D
        Phase modulation polynomial
    phx : Poly2D
        X derivative of the phase modulation polynomial
    phy : Poly2D
        Y derivative of the phase modulation polynomial
    M : list
        Diffraction orders

    Notes
    -----
    The surface shape is given by the shape attribute inherited from the
    parent class.

    Examples
    --------
    Create a 10mm x 10mm linear grating with a fringe every 15 microns:

    >>> g = RPPMask(shape=Rectangular(size=(10,10)),
    ...             phm=Poly2D([0,2*pi/15.,0]),
    ...             M=[1])
    """

    cdef public Poly2D phx
    cdef public Poly2D phy
    cdef public Poly2D phm
    cdef public list M

    def __init__(self, Poly2D phm=None, M=None, *args, **kwargs):
        """Initialize the RPPMask.

        Parameters
        ----------
        phm : Poly2D, optional
            Phase modulation polynomial
        M : list, optional
            Diffraction orders
        *args, **kwargs
            Arguments passed to parent class
        """
        super().__init__(*args, **kwargs)

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

    @cython.cdivision(True)
    cpdef list propagate(self, Ray incident_ray, double ni, double nr):
        """Calculate the propagation of a ray through a diffraction grating.

        This method implements the physical behavior of light when it encounters
        a reflective plane phase mask, considering both transmitted and reflected
        diffracted rays based on the surface's reflectivity.
        Parameters
        ----------
        incident_ray : Ray
            The incoming ray to be propagated
        ni : float
            Index of refraction of the incident medium
        nr : float
            Index of refraction of the refractive medium

        Returns
        -------
        list
            List of Ray objects representing the propagated rays (transmitted
            and/or reflected) after interaction with the grating

        Notes
        -----
        - All units are in millimeters
        - The intensity of the output rays is divided equally among the diffraction
          orders
        - Physically impossible rays (where oz2 < 0) are handled by returning the
          incident ray direction with zero intensity
        """
        cdef double l, rx, ry, rz, k, d, kx, ky, ox, oy, oz, oz2, M
        cdef Vector3d intersection_point, normal_vector, result_dir

        # Calculate intersection point and normal vector
        self.intersection_cy(incident_ray, intersection_point)
        self.normal_cy(intersection_point, normal_vector)

        # Convert wavelength to millimeters
        l = incident_ray.wavelength * 1e-3

        # Get incident ray direction components
        rx, ry, rz = incident_ray.direction

        # Calculate phase gradient
        cdef Vector2d K = Vector2d(
            self.phx.eval_cy(intersection_point(0), intersection_point(1)),
            self.phy.eval_cy(intersection_point(0), intersection_point(1))
        )

        # Calculate grating parameters
        k = K.norm()
        if k != 0:
            d = 2. * M_PI / k
            kx = K(0) / k
            ky = K(1) / k
        else:
            kx = ky = 0
            d = INFINITY

        ret = []
        ray_intensity = incident_ray.intensity / len(self.M)

        # Process each diffraction order
        for M in self.M:
            ox = rx + M * l/d * kx
            oy = ry + M * l/d * ky
            oz2 = rz**2 - 2 * M * l/d * (rx*kx + ry*ky) - (M * l/d)**2

            if oz2 < 0:
                print("Warning: eliminating physically impossible ray")
                ret.append(Ray.fast_init(
                    intersection_point,
                    incident_ray._direction,
                    0,  # Zero intensity for impossible rays
                    incident_ray.wavelength,
                    ni,
                    incident_ray.label,
                    incident_ray.draw_color,
                    None,
                    0,
                    self.id,
                    0,
                    incident_ray._parent_cnt + 1
                ))
            else:
                oz = copysign(sqrt(oz2), rz)

                # Handle transmitted rays
                if self.reflectivity != 1:
                    result_dir = Vector3d(ox, oy, oz)
                    ret.append(Ray.fast_init(
                        intersection_point,
                        result_dir,
                        ray_intensity,
                        incident_ray.wavelength,
                        ni,
                        incident_ray.label,
                        incident_ray.draw_color,
                        None,
                        0,
                        self.id,
                        0,
                        incident_ray._parent_cnt + 1
                    ))

                # Handle reflected rays
                if self.reflectivity != 0:
                    result_dir = Vector3d(ox, oy, -oz)
                    ret.append(Ray.fast_init(
                        intersection_point,
                        result_dir,
                        ray_intensity,
                        incident_ray.wavelength,
                        ni,
                        incident_ray.label,
                        incident_ray.draw_color,
                        None,
                        0,
                        self.id,
                        0,
                        incident_ray._parent_cnt + 1
                    ))

        return ret
