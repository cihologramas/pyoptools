from pyoptools.raytrace.component.component import Component
from pyoptools.raytrace.system.system import System
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.plane cimport Plane
from pyoptools.misc.cmisc.eigen cimport assign_tuple_to_vector3d, Vector3d
from math import isnan

import warnings
from functools import wraps


def deprecated_params(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        if "f" in kwargs:
            warnings.warn("The 'f' parameter is deprecated. Use 'focal_length' "
                          "instead.", DeprecationWarning, stacklevel=2)
            kwargs["focal_length"] = kwargs.pop("f")

        if "pupils" in kwargs:
            warnings.warn("The 'pupils' parameter is deprecated. Use "
                          "'pupil_config' instead.", DeprecationWarning,
                          stacklevel=2)
            kwargs["pupil_config"] = kwargs.pop("pupils")

        if "complete_trace" in kwargs:
            warnings.warn("The 'complete_trace' parameter is deprecated. Use "
                          "'show_internal_rays' instead.", DeprecationWarning,
                          stacklevel=2)
            kwargs["show_internal_rays"] = kwargs.pop("complete_trace")

        return func(*args, **kwargs)
    return wrapper


class IdealThickLens(System):
    """
    Class defined to create an Ideal Thick Lens.

    It is implemented as a subsystem because the propagation inside involves
    multiple rays: first surface -> first principal plane -> second principal
    plane -> second surface

    Parameters
    ----------
    lens_shape : object
        Shape of the lens (shape of both lens surfaces).
    lens_thickness : float
        Thickness of the lens (distance between the two surfaces).
    principal_plane_positions : tuple of float, optional
        Principal planes position. The position of each principal plane is
        measured from its corresponding surface. Default is (0., 0.).
    pupil_config : tuple or None, optional
        Configuration for the pupil. If None, no pupil is defined. Default is
        None. When provided, it should be a tuple containing:
        (pupil_distance, pupil_aperture_shape, pupil_reference_surface)

        - pupil_distance : float
            Distance of the pupil from the reference surface.
        - pupil_aperture_shape : object
            Shape of the pupil aperture.
        - pupil_reference_surface : bool
            Determines which lens surface is used as the reference for pupil
            positioning. If True, E1 is the reference surface. If False, E2 is
            the reference surface.
    focal_length : float, optional
        Focal length of the lens. Default is 100.
    show_internal_rays : bool, optional
        If True, the trace between the principal planes will be shown.
        If False, only the trace from and to the lens surfaces will be exact.
        Default is False.

    Notes
    -----
    The origin of the component is located at the midpoint between the two lens
    surfaces.

    E1 is positioned at -lens_thickness/2 along the z-axis.
    E2 is positioned at +lens_thickness/2 along the z-axis.

    Rays can travel in either direction through the lens:
    - For rays traveling in the positive z-direction, E1 is encountered first.
    - For rays traveling in the negative z-direction, E2 is encountered first.

    For the `pupil_config` parameter:
        The pupil_distance is always measured from the reference surface.

        When pupil_reference_surface is True (E1 is reference):
            - For rays entering through E1, it acts as an entrance pupil.
            - For rays exiting through E1, it acts as an exit pupil.
            - For rays entering through E2, it acts as an exit pupil.
            - For rays exiting through E2, it acts as an entrance pupil.
        When pupil_reference_surface is False (E2 is reference), the roles of
        E1 and E2 are reversed.

    The behavior of the pupil depends on the direction of ray propagation and
    the reference surface chosen.

    TODO: Convert from python to cython everything possible to increase
    performance
    """
    @deprecated_params
    def __init__(self,
                 lens_shape,
                 lens_thickness,
                 principal_plane_positions=(0., 0.),
                 pupil_config=None,
                 focal_length=100,
                 show_internal_rays=False):

        System.__init__(self, n=1)

        self.__E1__ = Plane(shape=lens_shape)
        self.__E2__ = Plane(shape=lens_shape)
        self.__H1__ = Plane(shape=lens_shape)
        self.__H2__ = Plane(shape=lens_shape)
        # Pupil setup
        if pupil_config is None:
            self.__PUP1__ = False
            self.__PUP2__ = False
        else:
            pupil_distance, pupil_aperture_shape, pupil_reference_surface = \
                pupil_config

            self.__PUP1__ = bool(pupil_reference_surface)
            self.__PUP2__ = not pupil_reference_surface

            if self.__PUP1__:
                self.__P1__ = Plane(shape=pupil_aperture_shape)
                self.__C5__ = Component(surflist=[(self.__P1__, (0, 0, 0), (0, 0, 0))])
                self.complist["P1"] = (self.__C5__,
                                       (0, 0, -lens_thickness/2+pupil_distance),
                                       (0, 0, 0))
            if self.__PUP2__:
                self.__P2__ = Plane(shape=pupil_aperture_shape)
                self.__C6__ = Component(surflist=[(self.__P2__,
                                                   (0, 0, 0),
                                                   (0, 0, 0))])
                self.complist["P2"] = (self.__C6__,
                                       (0, 0, lens_thickness/2+pupil_distance),
                                       (0, 0, 0))

        # Component setup
        self.__C1__ = Component(surflist=[(self.__E1__, (0, 0, 0), (0, 0, 0))])
        self.__C2__ = Component(surflist=[(self.__E2__, (0, 0, 0), (0, 0, 0))])
        self.__C3__ = Component(surflist=[(self.__H1__, (0, 0, 0), (0, 0, 0))])
        self.__C4__ = Component(surflist=[(self.__H2__, (0, 0, 0), (0, 0, 0))])

        self.complist["E1"]=(self.__C1__, (0, 0, -lens_thickness/2.), (0, 0, 0))
        self.complist["E2"]=(self.__C2__, (0, 0, lens_thickness/2.), (0, 0, 0))
        self.complist["H1"]=(self.__C3__,
                             (0, 0, -lens_thickness/2.+principal_plane_positions[0]),
                             (0, 0, 0))
        self.complist["H2"]=(self.__C4__,
                             (0, 0, lens_thickness/2.+principal_plane_positions[1]),
                             (0, 0, 0))

        self.focal_length = focal_length
        self.show_internal_rays = show_internal_rays

    def propagate_ray(self, Ray ri):
        """
        Propagate a ray through the ideal thick lens.

        Parameters
        ----------
        ri : Ray
            Incident ray to be propagated.

        Returns
        -------
        Ray
            The propagated ray with its child rays.
        """
        cdef double wav = ri.wavelength
        cdef double n = ri.n

        cdef Vector3d p_cy, d_cy
        cdef tuple[double, double, double] P, D

        cdef bint ST, ST2

        cdef Ray R, R_E1, R_E1_E2, R_E2, R_H1, R_H2

        _D, _C, S = self.distance(ri)

        # Determine if the ray enters through E1 or E2
        if S == self.__E1__:
            P1 = self.complist["P1"] if self.__PUP1__ else False
            E1 = self.complist["E1"]
            H1 = self.complist["H1"]
            H2 = self.complist["H2"]
            E2 = self.complist["E2"]
            P2 = self.complist["P2"] if self.__PUP2__ else False
        else:
            P1 = self.complist["P2"] if self.__PUP2__ else False
            E1 = self.complist["E2"]
            H1 = self.complist["H2"]
            H2 = self.complist["H1"]
            E2 = self.complist["E1"]
            P2 = self.complist["P1"] if self.__PUP1__ else False

        # Check if the ray passes through the entrance pupil
        if P1:
            C, P, D = P1
            assign_tuple_to_vector3d(P, p_cy)
            assign_tuple_to_vector3d(D, d_cy)
            R = ri.ch_coord_sys_f(p_cy, d_cy)
            PI = C["S0"][0].intersection(R)
            ST = isnan(PI[0])
        else:  # No entrance pupil check
            ST = False

        # Find the intersection point on the entrance surface E1
        C, P, D = E1
        assign_tuple_to_vector3d(P, p_cy)
        assign_tuple_to_vector3d(D, d_cy)
        R = ri.ch_coord_sys_f(p_cy, d_cy)
        PI = C["S0"][0].intersection(R)

        # Generate the ray from E1 to H1
        direction = R.direction if E1[1][2] < H1[1][2] else -R.direction
        R = Ray(origin=PI, direction=direction, wavelength=wav, n=n,
                intensity=ri.intensity)
        R_E1 = R.ch_coord_sys_inv_f(p_cy, d_cy, False)

        if ST:
            # The ray doesn't pass the entrance pupil. Draw it up to the entrance
            # surface
            R_E1.intensity = 0
            ri.add_child(R_E1)
            return ri

        # Propagate to H1
        C, P, D = H1
        assign_tuple_to_vector3d(P, p_cy)
        assign_tuple_to_vector3d(D, d_cy)
        R = R_E1.ch_coord_sys_f(p_cy, d_cy)
        PI = C["S0"][0].intersection(R)

        # Create the ray between H1 and H2
        direction = (0, 0, 1) if H1[1][2] < H2[1][2] else (0, 0, -1)
        R = Ray(origin=PI, direction=direction, wavelength=wav, n=n,
                intensity=ri.intensity)
        R_H1 = R.ch_coord_sys_inv_f(p_cy, d_cy, False)

        # Propagate to H2
        C, P, D = H2
        assign_tuple_to_vector3d(P, p_cy)
        assign_tuple_to_vector3d(D, d_cy)

        R = R_H1.ch_coord_sys_f(p_cy, d_cy)
        PI = C["S0"][0].intersection(R)

        # Calculate the refraction at H2
        _rx, _ry, rz = ri.direction
        FP = ri.direction * self.focal_length/ abs(rz)
        d = FP - PI
        if self.focal_length< 0:
            d = -d

        # Create the ray between H2 and E2
        direction = d if H2[1][2] < E2[1][2] else -d
        R = Ray(origin=PI,
                direction=direction,
                wavelength=wav,
                n=n,
                intensity=ri.intensity)
        R_H2 = R.ch_coord_sys_inv_f(p_cy, d_cy, False)

        # Propagate to E2
        C, P, D = E2
        assign_tuple_to_vector3d(P, p_cy)
        assign_tuple_to_vector3d(D, d_cy)
        R = R_H2.ch_coord_sys_f(p_cy, d_cy)
        PI = C["S0"][0].intersection(R)

        # Check the exit aperture
        PE2 = not isnan(PI[0])

        # Check the exit pupil
        if P2:
            C0, _P0, _D0 = P2
            R = R_H2.ch_coord_sys_f(p_cy, d_cy)
            PII = C0["S0"][0].intersection(R)
            ST2 = not (isnan(PII[0]) or isnan(PII[1]) or isnan(PII[2]))
        else:  # No exit pupil check
            ST2 = True

        ie2 = ri.intensity if ST2 and PE2 else 0

        direction = R.direction if H2[1][2] < E2[1][2] else -R.direction
        R = Ray(origin=PI, direction=direction, wavelength=wav, intensity=ie2, n=n)
        R_E2 = R.ch_coord_sys_inv_f(p_cy, d_cy, False)

        if self.show_internal_rays:
            ri.add_child(R_E1)
            R_E1.add_child(R_H1)
            R_H1.add_child(R_H2)
            R_H2.add_child(R_E2)
        else:
            R_E1_E2 = Ray(origin=R_E1.origin, direction=R_E2.origin-R_E1.origin,
                          wavelength=wav, n=n,
                          intensity=ri.intensity)
            ri.add_child(R_E1_E2)
            R_E1_E2.add_child(R_E2)

        return ri

    def distance(self, Ray ri):
        """
        Calculate the distance to the first intersection with the lens.

        Only C1 and C2, which are the entrance and exit surfaces of the system,
        are taken into account.

        Parameters
        ----------
        ri : Ray
            Incident ray.

        Returns
        -------
        tuple
            (distance, intersection point, intersected surface)
        """
        cdef float min_dist = float("inf")
        cdef Vector3d p_cy, d_cy
        cdef tuple[double, double, double] P, D
        min_pi = None
        min_surf = None

        for comp in [self.complist["E1"], self.complist["E2"]]:
            C, P, D = comp
            # Reorient the ray to the coordinate system of the element
            # and calculate the ray's path until it intersects with the element
            assign_tuple_to_vector3d(P, p_cy)
            assign_tuple_to_vector3d(D, d_cy)
            R = ri.ch_coord_sys_f(p_cy, d_cy)
            Dist = C.distance(R)

            if Dist[0] < min_dist:
                min_dist = Dist[0]
                min_pi = Dist[1]
                min_surf = Dist[2]

        return min_dist, min_pi, min_surf
