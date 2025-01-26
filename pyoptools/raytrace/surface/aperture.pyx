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
# Description:     Stops definition module
# Symbols Defined: Stop, StopCircAp, StopRecAp
# ------------------------------------------------------------------------------

from pyoptools.raytrace.surface.plane cimport Plane

from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_nan_to_vector3d


cdef class Aperture(Plane):
    """
    Class to define a surface with an aperture.

    The `Aperture` class is used to define stops in an optical system. It
    represents a surface that has both an external shape and an internal
    shape, which defines the aperture (hole) in the surface.

    Parameters
    ----------
    shape : Shape, optional
        An instance of a `Shape` subclass that defines the external shape
        of the stop. This parameter is inherited from the `Plane` class
        and passed through to it.
    ap_shape : Shape, optional
        An instance of a `Shape` subclass that defines the internal shape
        of the stop (the aperture shape). This parameter specifies the
        shape of the hole in the surface.

    Examples
    --------
    Creating an aperture surface with a rectangular external shape and
    a circular aperture::

        ap = Aperture(shape=Rectangular(size=(60, 60)), ap_shape=Circular(radius=2.5))

    Notes
    -----
    To create a stop component, it is recommended to use the `Stop` class
    from the `comp_lib` module. The `Stop` class creates the aperture surface
    and encapsulates it in a component that can be used within a `System`.
    """

    cdef public Shape ap_shape

    def __init__(self, ap_shape=None, *args, **kwargs):
        Plane.__init__(self, *args, **kwargs)
        self.ap_shape = ap_shape

        # Add items to the state list
        self.addkey("ap_shape")

    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:
        """
        """

        Plane._calculate_intersection(self,
                                      incident_ray,
                                      intersection_point)
        if self.ap_shape.hit_cy(intersection_point):
            assign_nan_to_vector3d(intersection_point)

    cpdef  list propagate(self, Ray ri, double ni, double nr):
        # Calculate the intersection point and the surface normal
        cdef Vector3d PI

        self._calculate_intersection(ri, PI)

        ret_ray = Ray.fast_init(PI,
                                ri._direction,
                                0,
                                ri.wavelength,
                                ri.n,
                                ri.label,
                                ri.draw_color,
                                ri,
                                ri.pop,
                                self.id,
                                0,
                                ri._parent_cnt+1)

        return [ret_ray]

#    cpdef list propagate(self, Ray ri, double ni, double nr):
#        """
#        Determine whether a ray continues propagating after intersecting the surface.#
#
#        This method overrides `Surface.propagate` to decide if a ray should continue
#        propagating or not after it intersects the surface. It checks if the
#        intersection point lies within the aperture defined by `ap_shape`.
#
#        Parameters
#        ----------
#        ri : Ray
#            The incoming ray to be propagated. This ray is in the coordinate system
#            of the surface.
#        ni : double
#            The refractive index of the medium from which the ray is incoming.
#            This parameter is not used in this implementation.
#        nr : double
#            The refractive index of the medium into which the ray would propagate.
#            This parameter is not used in this implementation.
#
#        Returns
#        -------
#        list of Ray
#            A list containing the resulting ray. The intensity of the ray is set to
#            zero if the intersection point does not lie within the aperture, indicating
#            that the ray does not continue propagating. Otherwise, the ray continues
#            with its original intensity.
#
#        Warnings
#        --------
#        This surface only checks if the ray continues propagating or not based on
#        the aperture shape. It does not calculate refraction or reflection.
#        Therefore, it must not be used to model lenses or mirrors.
#
#        Notes
#        -----
#        - The method calculates the intersection point between the incoming ray
#        and the surface and checks if this point is within the defined aperture.
#        - If the intersection point is within the aperture, the ray continues
#        propagating with its original intensity. If not, the ray's intensity
#        is set to zero, effectively stopping its propagation.
#        """
#        # Calculate the intersection point and the surface normal
#        cdef Vector3d PI
#
#        self._calculate_intersection(ri, PI)
#
#        # Check if the intersection point is within the aperture
#        if self.ap_shape.hit_cy(PI):
#            i = ri.intensity
#        else:
#            i = 0.
#
#        ret_ray = Ray.fast_init(PI,
#                               ri._direction,
#                                i,
#                                ri.wavelength,
#                                ri.n,
#                                ri.label,
#                                ri.draw_color,
#                                ri,
#                                ri.pop,
#                                self.id,
#                                0,
#                                ri._parent_cnt+1)
#
#        return [ret_ray]

    def __repr__(self):
        """
        Return a string representation of the Aperture optical surface.

        This method provides a detailed string representation of the `Aperture`
        object, including its shape, aperture shape, and other relevant attributes
        inherited from the `Plane` class.

        Returns
        -------
        str
            A string that represents the current state of the `Aperture` object.
        """
        return (
            f"Aperture(shape={repr(self.shape)}, "
            f"ap_shape={repr(self.ap_shape)}, "
            f"id={self.id})"
        )
