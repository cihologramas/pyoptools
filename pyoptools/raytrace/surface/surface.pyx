# ------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amezquita Orozco
# Description:     Surface definition module
# Symbols Defined: surface
# ------------------------------------------------------------------------------

from pyoptools.misc.pmisc import hitlist2int_list, hitlist2int, interpolate_g
from numpy import array, dot, isinf as npisinf, power, any, linspace, \
    meshgrid, indices, argwhere, tan, polyfit, arange, where, inf

import cython
from matplotlib.tri import Triangulation
from pyoptools.misc.picklable.picklable cimport Picklable

from pyoptools.raytrace.shape.circular cimport Circular
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.misc.function_2d.poly_2d.poly_2d import ord2i
from pyoptools.misc.lsq import polyfit2d
from scipy import interpolate
from numpy.ma import array as ma_array
from warnings import warn

from libc.math cimport sqrt, acos, cos, isnan, M_PI, INFINITY
from libc.stdlib cimport malloc, free

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_to_vector3d, \
    convert_vector3d_to_tuple, assign_nan_to_vector3d


cdef class Surface(Picklable):
    """"
    `Surface` is an abstract superclass for all optical surface objects. This
    class defines an API that all subclasses must implement.

    All subclasses of `Surface` share the following attributes and properties:

    Attributes
    ----------
    reflectivity : float
        A floating-point number between 0 and 1 representing the reflectivity
        of the surface. A value of 0 indicates a completely transparent
        surface, while a value of 1 indicates a completely reflective surface.
        Values between 0 and 1 represent beam splitters.
    shape : Shape
        An instance of the `Shape` class that defines the aperture shape of
        the surface.
    id : list
        A list of surface identifiers or tags. (Clarification of its specific
        purpose may be needed.)

    Properties
    ----------
    hit_list : list
        A read-only property that provides access to the `_hit_list`
        attribute, which stores information about rays that have impacted the
        surface. This list tracks the points of impact in the surface's
        reference system and keeps a pointer to the hitting ray, allowing
        retrieval of information such as intensity, wavelength, and other ray
        properties. The `reset` method clears this list.

    Notes
    -----
    Subclasses of `Surface` should implement specific methods to define
    geometry and behavior.
    """

    def __init__(self,
                 reflectivity=0.,
                 shape=None,
                 filter_spec=None):
        self.reflectivity = reflectivity
        self.shape = shape if shape is not None else Circular(radius=10.)
        self._hit_list = []

        # The id of each surface gets registered when the component is created,
        # and when the system is created

        self.id = []

        if filter_spec is None:
            filter_spec = ("nofilter",)

        if not isinstance(filter_spec, tuple) or len(filter_spec) < 1:
            raise ValueError("filter_spec must be a non-empty tuple")

        filter_type = filter_spec[0].lower()
        if filter_type == "nofilter":
            self.reflectivity_function = \
                <ReflectivityFunctionPtr>self.constant_reflectivity
        elif filter_type == "shortpass":
            if len(filter_spec) != 2:
                raise ValueError("Short-pass filter requires a cutoff parameter")
            self.cutoff = filter_spec[1]
            self.reflectivity_function = \
                <ReflectivityFunctionPtr>self.shortpass_reflectivity
        elif filter_type == "longpass":
            if len(filter_spec) != 2:
                raise ValueError("Long-pass filter requires a cutoff parameter")
            self.cutoff = filter_spec[1]
            self.reflectivity_function = \
                <ReflectivityFunctionPtr>self.longpass_reflectivity
        elif filter_type == "bandpass":
            if len(filter_spec) != 3:
                raise ValueError("Bandpass filter requires lower and upper "
                                 "cutoff parameters")
            self.lower_cutoff = filter_spec[1]
            self.upper_cutoff = filter_spec[2]
            self.reflectivity_function = \
                <ReflectivityFunctionPtr>self.bandpass_reflectivity
        else:
            raise ValueError(f"Invalid filter type: {filter_type}")

        Picklable.__init__(self, "reflectivity", "shape", "_hit_list", "id")

    @property
    def hit_list(self):
        """
        A read-only property that provides an immutable view of the hit list.

        Returns
        -------
        tuple
            A tuple containing the elements of the hit list, ensuring
            immutability.
        """
        # Return the list as a tuple to prevent modification
        return tuple(self._hit_list)

    cpdef list topo(self, list X, list Y):
        """
        Calculate the Z values for given X and Y coordinates on the Surface object.

        This method computes the topography (Z values) of the Surface object
        corresponding to the provided X and Y coordinates. It is primarily used
        for plotting the surface.

        Parameters
        ----------
        X : list of float
            A list of X coordinates on the surface for which Z values are to be
            calculated.
        Y : list of float
            A list of Y coordinates on the surface for which Z values are to be
            calculated.

        Returns
        -------
        list of float
            A list of Z values corresponding to each (X, Y) pair, representing the
            height of the surface at those coordinates.

        Note
        ----
        This method uses the topo_cy method, which must be overloaded in all
        Surface subclasses.
        """
        cdef int n = min(len(X), len(Y))
        cdef double* Z = <double*>malloc(n * sizeof(double))
        cdef int i
        cdef double x, y
        for i in range(n):
            x = X[i]
            y = Y[i]
            Z[i] = self.topo_cy(x, y)

        output = [Z[i] for i in range(n)]
        free(Z)
        return output

    def get_z_at_point(self, double x, double y):
        """
        Python wrapper for calculating Z coordinate at a single point on the Surface.

        This method provides a Python-accessible interface to topo_cy, allowing
        calculation of the surface height at a single (X, Y) coordinate pair.

        Parameters
        ----------
        x : double
            The X coordinate of the point.
        y : double
            The Y coordinate of the point.

        Returns
        -------
        double
            The Z coordinate (height) of the surface at the specified point.

        Note
        ----
        This is a wrapper around topo_cy which must be implemented in Surface
        subclasses.
        """

        return self.topo_cy(x, y)

    cdef double topo_cy(self, double x, double y) noexcept nogil:
        """
        Core implementation for calculating Z coordinate at a point on the Surface.

        This is the fundamental method for computing the surface height at a given
        point. It must be overloaded by all Surface subclasses to define their
        specific surface geometry.

        Parameters
        ----------
        x : double
            The X coordinate of the point.
        y : double
            The Y coordinate of the point.

        Returns
        -------
        double
            The Z coordinate (height) of the surface at the specified point.

        Note
        ----
        This is a Cython method marked with noexcept and nogil for performance.
        All Surface subclasses must override this method to implement their specific
        surface geometry calculations.
        """

        with gil:
            warn("Method topo_cy, from class Surface, should be overloaded" +
                 " in class "+self.__class__.__name__)

    cdef void _calculate_normal(self, Vector3d& intersection_point,
                                Vector3d& normal) noexcept nogil:
        """
        Calculate the normal vector at the point of intersection between a ray
        and the surface.

        This method computes the normal vector at the given point of
        intersection on the surface. The intersection point should already be
        calculated using the _calculate_intersection` method, and this method
        should not call it again to avoid redundant computations.

        Parameters
        ----------
        intersection_point : Vector3d
            An eigen::Vector3d instance representing the point of intersection
            on the surface, as returned by the `_calculate_intersection`
            method.
        normal : Vector3d
            An eigen::Vector3d instance representing the vector where the
            normal vector at the intersection point will be stored.

        Notes
        -----
        - This method should be overloaded in all subclasses of `Surface` to
        implement the specific normal calculation for each surface type.
        - The normal vector is calculated based on the surface's geometry at
        the given intersection point and must be normalized such that
        |normal| = 1.
        - The rest of the `pyoptools` code assumes that the normal vector is
        correctly normalized; therefore, this method must ensure that
        |normal| = 1.
        - This method should not call `_calculate_intersection` to avoid
        redundant computations. The intersection point will be provided by the
        caller.
        - This method should not be called directly. Use the `intersection` or
        `intersection_cy` methods instead.

        Raises
        ------
        NotImplementedError
            If the method is not overloaded in a subclass.
        """
        with gil:
            raise NotImplementedError("The _calculate_normal method must be "
                                      "overloaded in the subclass " +
                                      self.__class__.__name__)

    cdef inline void normal_cy(self, Vector3d& intersection_point,
                               Vector3d& normal) noexcept nogil:
        """
        Calls the surface-specific `_calculate_normal` method of the subclass.

        This method is designed to call the `_calculate_normal` method of the
        subclass to compute the normal vector at the given intersection point.
        It introduces a slight redundancy, but helps maintain a consistent API
        design similar to the `intersection` method.

        This Cython method should be used exclusively by Cython code. For
        Python code, use the corresponding Python method `normal`, which
        internally calls this Cython method.

        Parameters
        ----------
        intersection_point : Vector3d
            An eigen::Vector3d instance representing the point of intersection
            on the surface.
        normal : Vector3d
            An eigen::Vector3d instance where the normal vector at the
            intersection point will be stored.

        Notes
        -----
        - This method provides a consistent interface for normal vector
        computation, similar to how the `intersection` method is structured.
        - The Cython method `normal_cy` should be used by all Cython code to
        ensure optimal performance and direct access to internal functionality.
        - The Python method `normal` is intended for use in Python code and
        calls this Cython method to perform the actual computation.
        """
        self._calculate_normal(intersection_point, normal)

    def normal(self, intersection_point):
        """
        Calculate the normal vector at a specific intersection point.

        This method returns the normalized normal vector at a given
        intersection point. The method must be overloaded in all `Surface`
        subclasses to provide surface-specific geometry calculations.

        Parameters
        ----------
        intersection_point : tuple or list
            The intersection point as a tuple or list (x, y, z) in the
            coordinate system of the surface.

        Returns
        -------
        tuple
            A tuple (dx, dy, dz) representing the normalized normal vector at
            the intersection point.

        Notes
        -----
        - The normal vector is normalized and perpendicular to the surface at
        the given intersection point.
        - This method should be overloaded in all subclasses of `Surface` to
        provide specific implementations based on the surface geometry.
        """
        cdef Vector3d intersection_point_vector, normal_vector

        assign_to_vector3d(intersection_point, intersection_point_vector)

        self.normal_cy(intersection_point_vector, normal_vector)

        return convert_vector3d_to_tuple(normal_vector)

    cdef void _calculate_intersection(self,
                                      Ray incident_ray,
                                      Vector3d& intersection_point) noexcept nogil:

        """
        Calculate the point of intersection between a ray and a surface.

        This method computes the point of intersection between the surface
        and the incident ray in the surface's coordinate system. It must be
        overloaded in all subclasses of `Surface` to implement the
        geometry-specific intersection calculation for each surface type.

        Parameters
        ----------
        incident_ray : Ray
            The incident ray, which must be in the coordinate system of the
            surface.
        intersection_point : Vector3d
            An eigen::Vector3d instance where the intersection point will be
            stored.

        Notes
        -----
        - This method does not perform any checks for the aperture; such checks
        are automatically handled by the `intersection_cy` method.
        - This method should not be called directly. Use the `intersection` or
        `intersection_cy` methods instead.
        - This function must be overloaded in all subclasses of `Surface` to
        provide the specific intersection calculation for each type of surface.

        Raises
        ------
        NotImplementedError
            If the method is not overloaded in a subclass.
        """
        with gil:
            raise NotImplementedError("The _calculate_intersection method must"
                                      " be overloaded in the subclass " +
                                      self.__class__.__name__)

    cdef inline void intersection_cy(self, Ray incident_ray,
                                     Vector3d& intersection_point):
        """
        Compute the intersection point of a ray with the surface.

        This method calculates the intersection point of the given incident ray
        with the surface and stores the result in the provided
        `intersection_point` vector. It performs multiple validation checks:
        1. Whether the ray intersects within the defined boundary of the surface
        2. Whether the intersection occurs in front of the ray origin (forward
           propagation)

        If any validation fails, the method sets `intersection_point` to a vector
        of NaN values.
        Parameters
        ----------
        incident_ray : Ray
            The incident ray that intersects with the surface.
        intersection_point : Vector3d&
            Reference to an eigen::Vector3d instance where the intersection point will
            be stored. If the intersection is invalid, this will be set to
            (NaN, NaN, NaN).
        Notes
        -----
        - This is a Cython-specific, low-level function intended for internal use.
        - The method doesn't return a value; instead, it modifies the passed
          `intersection_point`.
        - Forward propagation is enforced by ensuring the intersection occurs in the
          direction of the ray (positive distance).
        """
        # Calculate the intersection point between ray and surface
        self._calculate_intersection(incident_ray, intersection_point)

        # Validation 1: Check if the intersection point is inside the surface boundary
        if not self.shape.hit_cy(intersection_point):
            assign_nan_to_vector3d(intersection_point)
            return

        # Validation 2: Check if the intersection is in front of the ray
        # (forward propagation)
        # Calculate distance along ray direction from origin to intersection
        cdef double dist = (intersection_point -
                            incident_ray._origin).dot(incident_ray._direction)

        # If distance is negative or too small (numerical precision threshold),
        # the intersection is invalid
        if dist < 1.e-10:
            assign_nan_to_vector3d(intersection_point)

    def intersection(self, Ray incident_ray):
        """
        Calculate the point of intersection between a ray and this surface.

        This method returns the point of intersection between the surface and
        the incident ray. The intersection point is calculated in the
        coordinate system of the surface.

        This calculation is performed considering only the ray and this surface
        in isolation. It does not take into account any other optical elements
        that may be present in the system. For this method, the only existing
        elements are the ray and the surface.

        If there is no point of intersection, such as when the ray is outside
        the surface's boundary, the method returns a tuple (NaN, NaN, NaN).

        Parameters
        ----------
        incident_ray : Ray
            The incident ray, which must be in the coordinate system of the
            surface.

        Returns
        -------
        tuple
            A tuple (x, y, z) containing the coordinates of the point of
            intersection between the ray and the surface. If there is no
            intersection, it returns (NaN, NaN, NaN).
        """
        cdef Vector3d intersection_point
        self.intersection_cy(incident_ray, intersection_point)
        return convert_vector3d_to_tuple(intersection_point)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef double distance_cy(self, Ray incident_ray, Vector3d& intersection_point):
        """
        Calculate the distance from the ray's origin to the intersection point
        on the surface.

        This method computes the distance from the origin of the given incident
        ray to the intersection point on the surface. It calculates the
        intersection point and stores it in `intersection_point`. If there
        is no valid intersection, the method sets the intersection point to
        NaN and returns infinity for the distance.

        Parameters
        ----------
        incident_ray : Ray
            The incident ray for which the distance to the intersection point
            is calculated.
        intersection_point : Vector3d
            An eigen::Vector3d instance where the coordinates of the
            intersection point will be stored.

        Returns
        -------
        double
            The distance from the ray's origin to the intersection point. If
            there is no intersection, the method returns infinity.

        Notes
        -----
        - If the intersection point is behind the ray's origin or too close
        due to rounding errors, the distance is considered infinite.
        - The method handles cases where the surface is ahead of or behind the
        ray, or when there is no valid intersection.
        - Distances greater than 1e10 are considered as infinity to avoid
        numerical precision issues with extremely distant intersections.
        """
        # Calculate the intersection point and store it in intersection_point
        self.intersection_cy(incident_ray, intersection_point)

        cdef double dist

        # Check if any component of the intersection point is NaN, indicating
        # no intersection
        if  (isnan(intersection_point(0)) or
             isnan(intersection_point(1)) or
             isnan(intersection_point(2))):
            dist = INFINITY
        else:
            # Calculate the distance as the dot product of the vector from ray origin
            # to intersection point with the ray direction
            dist = (intersection_point - incident_ray._origin). \
                dot(incident_ray._direction)

        # Treat extremely large distances (>1e10) as infinity
        # This helps avoid numerical precision issues with very distant intersections
        if dist > 1e10:
            dist = INFINITY

        return dist

    def distance(self, incident_ray):
        """
        Calculate the distance propagated by a ray until it intersects with a
        surface.

        This method returns the distance traveled by a ray from its origin to
        its point of intersection with a surface.

        Parameters
        ----------
        incident_ray : Ray
            The incident ray, which must be in the coordinate system of the
            surface.

        Returns
        -------
        tuple
            A tuple `(distance, point_of_intersection, surface)`, where:

            distance : float
                The distance from the ray's origin to the point of intersection
                with the surface.
            point_of_intersection : tuple
                The coordinates of the intersection point in the coordinate
                system of the surface.
            surface : Surface
                A reference to the surface object.
        """
        cdef Vector3d intersection_point
        cdef double dist

        dist = self.distance_cy(incident_ray, intersection_point)
        return dist, convert_vector3d_to_tuple(intersection_point)

    @cython.cdivision(True)
    cpdef list propagate(self, Ray incident_ray, double ni, double nr):
        """
        Calculate the refracted ray on a surface.

        This method calculates the ray refracted (or both refracted and
        reflected rays) on a surface.

        Parameters
        ----------
        incident_ray : Ray
            The incident ray, which must be in the coordinate system of the
            surface.
        ni : double
            The refractive index of the incident medium.
        nr : double
            The refractive index of the refracted medium.

        Returns
        -------
        list
            A list containing the rays resulting from the ray-surface
            interaction.
        """

        cdef double incidence_angle, gamma, gamma1, ddot, reflect
        cdef Vector3d intersection_point, normal_vector, A, A1, S1, S2, S3

        self.intersection_cy(incident_ray, intersection_point)

        self.normal_cy(intersection_point, normal_vector)

        # The information about the hit list is added in the
        # system.propagate_ray method.
        # Do not uncomment the following line:
        # self._hit_list.append((intersection_point, incident_ray))

        # S1 = array(incident_ray.direction)*ni
        S1 = incident_ray._direction * ni

        # Calculate the incident angle
        # ddot = vector3_dot_product(&S1, &normal_vector)
        ddot = S1.dot(normal_vector)
        incidence_angle = (acos(ddot/ni))

        # some times because rounding errors |incident_ray.dir|>1. In this cases
        # I=nan
        if isnan(incidence_angle):
            incidence_angle = 0.

        # Take the correct normal
        if incidence_angle > M_PI/2:
            # P=-P
            normal_vector = -normal_vector
            ddot = S1.dot(normal_vector)
            incidence_angle = (acos(ddot/ni))

        gamma = nr*sqrt((ni/nr*cos(incidence_angle))**2 - (ni/nr)**2+1.) \
            - ni*cos(incidence_angle)

        # A=gamma*P
        # vector3_times_scalar(&normal_vector,gamma,&A)
        A = normal_vector * gamma
        # S2=S1+A
        # add_vector3(&S1,&A,&S2)
        S2 = S1 + A

        # If the refractive index (nr) is negative, an error is raised. This
        # should not happen.
        # If reflectivity is 0, the optical surface does not function as a beam
        # splitter.
        # If there is total internal reflection, the surface should also not
        # function as a beam splitter.

        # If the optical surface is a beam splitter, this method should return
        # a list containing both the transmitted and reflected rays:
        # [Transmitted ray, Reflected ray].

        # normalize_vector3(&S2)
        S2.normalize()

        # This could be used in the future to define spectral filters, such as
        # dichroic reflectors.
        # For example:
        # try:
        #     reflect = self.reflectivity(incident_ray.wavelength)
        # except TypeError:
        #     reflect = self.reflectivity

        reflect = self.reflectivity_function(self, incident_ray.wavelength)

        if reflect > 1 or reflect < 0:
            raise ValueError

        if ni<0 or nr<0:
            # This case should never happen
            raise ValueError, \
                f"Negative refractive index detected. ni={ni} nr={nr}"

        elif (reflect == 0) and not(isnan(S2(0)) or isnan(S2(1)) or isnan(S2(2))):
            # Normal refraction case
            return [Ray.fast_init(intersection_point,
                                  S2,
                                  incident_ray.intensity,
                                  incident_ray.wavelength,
                                  nr, incident_ray.label,
                                  incident_ray.draw_color,
                                  None,
                                  0.,
                                  self.id,
                                  0,
                                  incident_ray._parent_cnt+1)]
        elif (isnan(S2(0)) or isnan(S2(1)) or isnan(S2(2))):
            # Total internal reflection case
            gamma1 = -2.*ni*cos(incidence_angle)

            # A1=gamma1*P
            A1 = normal_vector * gamma1
            # S3=S1+A1
            S3 = S1+A1

            S3.normalize()  # S3=S3/sqrt(dot(S3,S3))

            return [Ray.fast_init(intersection_point,
                                  S3,
                                  incident_ray.intensity,
                                  incident_ray.wavelength,
                                  ni,
                                  incident_ray.label,
                                  incident_ray.draw_color,
                                  None,
                                  0.,
                                  self.id,
                                  0,
                                  incident_ray._parent_cnt+1)]

        else:
            # BeamSplitter case
            gamma1 = -2.*ni*cos(incidence_angle)
            # A1=gamma1*P

            A1 = normal_vector * gamma1

            # S3=S1+A1
            S3 = S1 + A1

            S3.normalize()  # S3=S3/sqrt(dot(S3,S3))

            if reflect != 1.:
                return [Ray.fast_init(intersection_point,
                                      S2,
                                      incident_ray.intensity*(1.-reflect),
                                      incident_ray.wavelength,
                                      nr,
                                      incident_ray.label,
                                      incident_ray.draw_color,
                                      None,
                                      0.,
                                      self.id,
                                      0,
                                      incident_ray._parent_cnt+1),
                        Ray.fast_init(intersection_point,
                                      S3,
                                      incident_ray.intensity*reflect,
                                      incident_ray.wavelength,
                                      ni,
                                      incident_ray.label,
                                      incident_ray.draw_color,
                                      None,
                                      0,
                                      self.id,
                                      0,
                                      incident_ray._parent_cnt+1)]
            else:
                return [Ray.fast_init(intersection_point,
                                      S3,
                                      incident_ray.intensity,
                                      incident_ray.wavelength,
                                      ni,
                                      incident_ray.label,
                                      incident_ray.draw_color,
                                      None,
                                      0.,
                                      self.id,
                                      0,
                                      incident_ray._parent_cnt+1)]

    cdef double constant_reflectivity(self, double wavelength) noexcept nogil:
        return self.reflectivity

    cdef double longpass_reflectivity(self, double wavelength) noexcept nogil:
        return self.reflectivity if wavelength < self.cutoff else 0.0

    cdef double shortpass_reflectivity(self, double wavelength) noexcept nogil:
        return self.reflectivity if wavelength > self.cutoff else 0.0

    cdef double bandpass_reflectivity(self, double wavelength) noexcept nogil:
        return (0.0 if self.lower_cutoff < wavelength < self.upper_cutoff
                else self.reflectivity)

    # Note: We could add a function that runs a custom python function so it
    # could be provided at runtime

    cpdef pw_propagate1(self, Ray ri, ni, nr, rsamples, isamples, knots):
        '''Method to calculate wavefront emerging from the surface

        .. Note::
            This method needs to be checked, because the incoming plane
            wave has not a constant intensity. The ray distribution is
            affected by the surface topography. it is better to use the
            slow version pw_propagate.

        This method calculates the wavefront emerging from the optical surface
        when illuminated by an unity amplitude plane wave.
        The k vector of the incoming PW is given by the direction of ri.
        The wavelength of the PW is given by the wavelength of ri.
        The origin of ri is not used at all.

        The input and output planes are located at z=0 in the surface
        coordinated system.

        Note: this method needs to be checked, because the incoming plane wave
        has not a constant intensity. The ray distribution is affected by the
        surface topography. it is better to use the slow version pw_propagate1.

        Note: The ray comes from the negative side. Need to change this

        Arguments:


             ri -- isurface.pyncident ray

             ni -- refraction index in the incident media n.

             nr -- refraction index in the refracted media

        ri must be in the coordinate system of the surface
        '''

        # Calculate the position maximum and minimum z values for the surface
        X, Y, Z = self.shape.mesh(topo=self.topo)

        # Calculate a mesh where the test rays will be shot from
        # ,  size=(xmin, xmax, ymin, ymax))
        X, Y, H = self.shape.mesh(ndat=rsamples)
        X, Y, Z = self.shape.mesh(ndat=rsamples, topo=self.topo)

        # The ecuation where the rays are shooted to is Z=0
        # N_ = array((0, 0, 1))

        # The ecuation where the rays are shooted from is
        # X*dx+Y*dy+Z*dz=0, where dx,dy,dx are the plane normal.
        N_p = ri.direction

        ix, iy = indices(X.shape)

        # Create a list of indexes to be able to select the points that are
        # going to be used as spline generators, and as control points
        idx = where((ix+iy) % 2 == 0, False, True)

        X = X.flatten()
        Y = Y.flatten()
        Z = Z.flatten()
        H = H.flatten()
        idx = idx.flatten()
        # Get the rays that hit the surface
        li = H.nonzero()[0]
        # rin=Ray( dir=L1, wavelength=wavelength)
        dir = ri.direction
        S1 = dir*ni
        # Intersection point and optical path for the incident rays
        xi = []
        yi = []
        zi = []
        ii = []

        for i in li:
            # Destination of the ray
            P1 = array((X[i], Y[i], Z[i]))

            # Normal to the destination point
            Np = self.normal(P1)

            # Calculate incident and refracted angle. To optimize the routine,
            # only normal diffracted ray is taken into account

            # Calculate the incident angle
            I = (acos(dot(S1, Np)/ni))

            # Take the correct normal
            if I > M_PI/2:
                Np = -Np
                I = (acos(dot(S1, Np)/ni))
            # Calculate the diffracted direction

            gamma = nr*sqrt(power(ni/nr*cos(I), 2) -
                            power(ni/nr, 2)+1.) - ni*cos(I)
            A = gamma*Np
            S2 = nr*(S1+A)
            L2 = (S2/sqrt(dot(S2, S2)))
            ###

            # Calculate the distance between P1 ant the plane X*dx+Y*dy+Z*dz=0,
            # along the N_p direction
            din = -dot(N_p, -P1)

            # Calculate the distance between P1 Z=0 along the L2 line
            ddi = -P1[2]/L2[2]

            # Calculate the intersection point
            PIR = P1+ddi*L2

            # Calculate optical path
            d = din*ni+ddi*nr
            if d != inf:
                x, y, _z = PIR
                xi.append(x)
                yi.append(y)
                zi.append(d)
                ii.append(idx[i])
        xi = array(xi)
        yi = array(yi)
        zi = array(zi)
        ii = array(ii)

        xmax = xi.max()  # =abs(array(xi+yi)).max()
        ymax = xmax
        xmin = xi.min()  # -ymax
        ymin = xmin  # -xmax

        xx = linspace(xmin, xmax, isamples[0])
        yy = linspace(ymin, ymax, isamples[1])

        # Use only half of the samples to create the Spline,
        isp = argwhere(ii)  # (ii == True)
        ich = argwhere(not ii)  # (ii == False)

        xsp = xi[isp]
        ysp = yi[isp]
        zsp = zi[isp]

        xch = xi[ich]
        ych = yi[ich]
        zch = zi[ich]

        # Distribute homogeneously the knots
        xk = linspace(xsp.min(), xsp.max(), knots)
        yk = linspace(ysp.min(), ysp.max(), knots)

        # LSQBivariateSpline using some knots gives smaller error than
        # SmoothBivariateSpline
        di = interpolate.LSQBivariateSpline(xsp, ysp, zsp, xk[1:-1], yk[1:-1])
        # di=interpolate.SmoothBivariateSpline(xsp, ysp, zsp)

        # Evaluate error
        zch1 = di.ev(xch, ych)
        er = (zch.flatten()-zch1).std()

        # I, xe, ye=histogram2d(yi, xi, (xx, yy))
        I = hitlist2int(xi, yi, xi, xx, yy)
        d = ma_array(di(xx, yy).transpose(), mask=I.mask)

        # XD, YD=meshgrid(xx, yy)

        # d1=griddata(xi,  yi,  zi,  xx, yy)
        return I, d, er

    cpdef pw_propagate(self, Ray ri, ni, nr, rsamples, shape, order, z):
        '''Method to calculate wavefront emerging from the surface


        This method calculates the wavefront emerging from the optical surface
        when illuminated by an unity amplitude plane wave.
        The k vector of the incoming PW is given by the direction of ri.
        The wavelength of the PW is given by the wavelength of ri.
        The origin of ri is not used at all.

        Arguments:


             ri -- surface incident ray

             ni -- refraction index in the incident media n.

             nr -- refraction index in the refracted media

             rsamples -- number of rays used to sample the surface

             shape -- shape of the output wavefront

             order -- order of the polynomial fit

             z -- Z position of the input and output plane. The origin is the
                  surface vertex

        ri must be in the coordinate system of the surface
        Note: The ray comes from the negative side. Need to change this
        '''
        from ray_trace.surface import Plane

        xi, yi, zi = self.pw_propagate_list(ri, ni, nr, rsamples, z)

        xmin, xmax, ymin, ymax = self.shape.limits()
        xx = linspace(xmin, xmax, shape[0])
        yy = linspace(ymin, ymax, shape[1])

        # I, xe, ye=histogram2d(yi, xi, (xx, yy))
        I = hitlist2int(xi, yi, xi, xx, yy)

        p, er = polyfit2d(xi, yi, zi, order=order)

        d = ma_array(p.evalvv(xx, yy), mask=I.mask)
        # d,er=interpolate_g(xi,yi,zi,xx,yy,knots=knots,error=True,mask=I.mask)

        # XD, YD=meshgrid(xx, yy)

        # d1=griddata(xi,  yi,  zi,  xx, yy)
        return I, d, er

    cpdef pw_propagate_list(self, Ray ri, ni, nr, rsamples, z):
        '''Method to calculate wavefront emerging from the surface

        This method calculates samples of the  wavefront emerging from
        the optical surface when illuminated by an unity amplitude plane
        wave.
        The k vector of the incoming PW is given by the direction of ri.
        The wavelength of the PW is given by the wavelength of ri.
        The origin of ri is not used at all.

        The returned value is a list containing the x,y coordinates of the ray
        list in the output surface, and the optical path at such point.

        To create an output matrix, this values must be interpolated.



        Arguments:


             ri -- incident ray

             ni -- refraction index in the incident media n.

             nr -- refraction index in the refracted media

             rsamples -- number of rays used to sample the surface (Tuple)

             z -- Z position of the input and output plane. The origin is the
                  surface vertex

        ri must be in the coordinate system of the surface
        Note: The ray comes from the negative side. Need to change this
        '''
        from pyoptools.raytrace.surface import Plane

        # plane where the aperture is located
        Za = z

        # Create an array of rays to simulate the plane wave
        dir = ri.direction
        xmin, xmax, ymin, ymax = self.shape.limits()

        # Calculate the position maximum and minimum z values for the surface
        # X, Y, H=self.shape.mesh(ndat=rsamples)
        # X, Y, H=self.shape.mesh(ndat=rsamples)

        # X;Y;Z coordinates of the points the rays are aiming to  in the
        # aperture plane. The 20% increase in size, is to assure that
        # all the surface is sampled even if the rays are tilted
        dx = float(xmax-xmin)*.1
        dy = float(ymax-ymin)*.1

        Xa = linspace(xmin-dx, xmax+dx, int(rsamples[0]*1.2))
        Ya = linspace(ymin-dy, ymax+dy, int(rsamples[1]*1.2))
        # print Xa,Ya
        Xa, Ya = meshgrid(Xa, Ya)

        #
        Xa = Xa.flatten()
        Ya = Ya.flatten()

        # The ecuation where the rays are shooted from is
        # X*dx+Y*dy+Z*dz=0, where dx,dy,dx are the plane normal.
        Z = -(Xa*dir[0]+Ya*dir[1])/dir[2]+Za
        # H=H.flatten()

        # li=H.nonzero()[0]
        pl = Plane()
        xi = []
        yi = []
        zi = []
        for i in range(len(Xa)):

            P0 = array((Xa[i], Ya[i], Z[i]))

            ri = Ray(origin=P0, direction=dir, wavelength=ri.wavelength)
            # if the ray do not intersect the surface, break current iteration
            if any(npisinf(self.intersection(ri))):
                continue

            # print self.intersection(ri)
            rd = self.propagate(ri, ni, nr)

            # take only the transmitted ray
            rd = rd[0]

            # Translate rd, to put it in the aperture reference plane
            rd.origin[2] = rd.origin[2]-Za
            di = self.distance_s(ri)[0]
            dr = pl.distance_s(rd)[0]

            PI = pl._intersection(rd)

            # Calculate optical path
            d = di*ni+dr*nr
            # print d,di,dr,PI
            if d != inf:
                x, y, z = PI
                xi.append(x)
                yi.append(y)
                zi.append(d)
        xi = array(xi)
        yi = array(yi)
        zi = array(zi)

        return xi, yi, zi

    cpdef wf_propagate(self, wf, ni, nr, samples, shape, knots):
        '''Method to calculate wavefront emerging from the surface


        This method calculates the wavefront emerging from the optical surface
        when illuminated by an arbitrary wavefront.

        The input and output planes are located at z=0 in the surface
        coordinated system.

        Arguments:

             wf -- Field instance containing the incoming wavefront

             ni -- refraction index in the incident media n.

             nr -- refraction index in the refracted media.

             samples -- Tuple containing the number of rays used to sample the field.

             shape -- Tuple containing the shape of the output field
        '''
        from ray_trace.surface import Plane
        # Get the wavefront ray representation
        # TODO: This representation only takes into account the phase, but not
        #       the intensity.
        #       This has to be fixed, because in practice this is not OK
        rays = wf.rayrep(samples[0], samples[1])

        # TODO: check which rays pass inside the aperture

        # rin=Ray( dir=L1, wavelength=wavelength)

        # Intersection point and optical path for the incident rays

        xi = []
        yi = []
        zi = []

        pl = Plane()
        for ri in rays:
            # Calculate the intersection point
            rd = self.propagate(ri, ni, nr)
            # Take only de transmitted ray
            rd = rd[0]

            # Incident ray propagation distance until the optical surface
            di = self.distance(ri)[0]

            # Refracted ray propagation until the output surface (plane Z=0)
            # The distance method is not used, because it eliminates the negative
            # propagation values

            PI = pl._intersection(rd)
            dr = dot(PI-rd.origin, rd.direction)

            d = di*ni+dr*nr+ri.optical_path_parent()

            if d != inf:
                x, y, _z = PI
                xi.append(x)
                yi.append(y)
                zi.append(d)
            else:
                print(ri)

        d = interpolate_g(xi, yi, zi, samples=shape,
                          knots=knots, error=False, mask=None)

        return d

    def __repr__(self):
        # ~ '''Return an string with the representation of the optical surface
        # ~
        # ~ It must be overloaded in all subclasses
        # ~ '''

        return "OptSurf(reflectivity="+str(self.reflectivity)+")"

    cpdef reset(self):
        ''' Remove information from the hit_list
        '''
        self._hit_list = []

    cpdef pw_cohef(self, ni, nr, ilimit, slimit, step, order, rsamples, zb):
        '''Method to generate the taylor polinomial coheficients that can be
        used to obtain the phase and intensity for different plane wave
        inclinations.

        Notes:
               - The pupil is normalized to the radius of the lens
               - This method assumes a circular shaped pupil
               - The phase is returned as an optical path
               - The intensity is normalized to 1

        Arguments:

        =========  ======================================================
        ni,nr      Refraction index from the incident and refracted sides
        ilimit     Inferior limit for incidence angle of the plane wave
                   in radians
        slimit     Superior limit for incidence angle of the plane wave
                   in radians
        step       Step to be used to generate the interpolation data
        order      Order of the Taylor interpolation used
        rsamples   Tuple containing the number of ray samples to be used
                   in each direction
        zb         Z position of the plane where the measurementas are
                   made. The origin is the vertex of the surface. If
                   None, an estimate position is taken.
        =========  ======================================================

        '''
        xm = self.shape.limits()[1]
        # print xm
        # Get the z position of the border of the surface

        if zb is None:
            zm = self.topo(xm, 0)
        else:
            zm = zb

        # Optical path coheficient list
        opcl = []

        # Intensity coheficient list
        icl = []

        xd = []
        for i in arange(ilimit, slimit, step):
            xd.append(i)

            x, y, d = self.pw_propagate_list(
                Ray(direction=(tan(i), 0, 1)), ni, nr, rsamples=rsamples, z=zm)

            # Normalize the pupil
            x = x/xm
            y = y/xm

            # Get the optical path polynomial coheficients
            pf, _ef = polyfit2d(x, y, d, order=order)

            # get the intensity data
            xi, yi, I = hitlist2int_list(x, y)

            # Get the intensity polynomial coheficients
            pi, _ei = polyfit2d(xi, yi, I, order=order)

            # TODO: Print something if error is too big
            # if ei>0.001: print ""
            # if ef>1e-6: print ""

            opcl.append(pf.cohef.flatten())

            icl.append(pi.cohef.flatten())

        dph = array(opcl)
        di = array(icl)

        # Lista de los coheficientes de los polinomios para generar los
        # coheficientes
        phcohef = []

        # get the number of coheficients the polinomial expansion has
        ncohef = ord2i(order)

        for i in range(ncohef):
            # figure()
            # plot(xd,d[:,i])
            pof = polyfit(xd, dph[:, i], 15)
            phcohef.append(pof)

        icohef = []

        for i in range(ncohef):
            # figure()
            # plot(xd,d[:,i])
            pof = polyfit(xd, di[:, i], 15)
            icohef.append(pof)

        return phcohef, icohef, zm

    def polylist(self):

        points = []

        X, Y = self.shape.pointlist()
        Z = self.topo(X, Y)

        for i in range(len(X)):
            points.append((X[i], Y[i], Z[i]))

        # Need to find a better way to do this not using delaunay
        # or maybe to generate all using triangulations????
        tri = Triangulation(X, Y)
        trip = tri.triangles
        return points, trip
