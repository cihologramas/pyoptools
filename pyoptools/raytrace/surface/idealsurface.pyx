
"""
Modulo con clases y funciones auxiliares.
"""


from numpy import array, float64
from pyoptools.misc.definitions import inf_vect

from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_nan_to_vector3d, \
                                        assign_doubles_to_vector3d
from libc.math cimport fabs, isnan

from pyoptools.misc.cmisc.eigen cimport Vector3d, assign_nan_to_vector3d, \
                                        assign_doubles_to_vector3d

cdef class IdealSurface(Surface):
    """
    A class that represents an ideal surface. It is used to create
    mathematically ideal lenses.

    Clase que representa una superficie ideal. Se utiliza para crear
    lentes ideales.
    """

    cdef public double f

    def __init__(self, f=100, *args, **kwargs):
        """
        f represents the focal length of the ideal surface
        """
        Surface.__init__(self, *args, **kwargs)

        self.f = f

        # Add attributes to the state list
        self.addkey("f")

    cdef inline double topo_cy(self, double x, double y) noexcept nogil:
        return 0


    cdef void _calculate_intersection(self,
                                    Ray incident_ray,
                                    Vector3d& intersection_point) noexcept nogil:
        """
        Calculate the intersection point between a ray and the XY plane (z = 0).

        Parameters
        ----------
        incident_ray : Ray
            The incoming ray for which the intersection with the plane is to be calculated.
        intersection_point : Vector3d&
            Output: intersection point in the surface coordinate system.
            If there is no valid intersection, this vector is set to NaN.

        Notes
        -----
        - The plane is z = 0 with normal (0, 0, 1).
        - If the ray is parallel to the plane (direction z component is zero),
        no intersection is reported (NaN).
        """

        cdef double x1, y1, z1
        cdef double dx, dy, dz
        cdef double u
        cdef double X, Y, Z

        x1 = incident_ray._origin(0)
        y1 = incident_ray._origin(1)
        z1 = incident_ray._origin(2)

        dx = incident_ray._direction(0)
        dy = incident_ray._direction(1)
        dz = incident_ray._direction(2)

        # Parallel to plane (or numerically extremely close)
        if dz == 0.0:
            assign_nan_to_vector3d(intersection_point)
            return

        # Solve z1 + u*dz = 0  ->  u = -z1/dz
        u = -z1 / dz

        X = x1 + u * dx
        Y = y1 + u * dy
        Z = z1 + u * dz  # should be 0 (up to floating error)

        assign_doubles_to_vector3d(X, Y, Z, intersection_point)


    cdef void _calculate_normal(self,
                                Vector3d& intersection_point,
                                Vector3d& normal) noexcept nogil:
        """
        Calculate the normal vector to the XY plane (z = 0).

        Parameters
        ----------
        intersection_point : Vector3d&
            Point on the surface (unused; normal is constant for a plane).
        normal : Vector3d&
            Output: the surface normal (unit length).
        """
        # Constant normal for z=0 plane
        assign_doubles_to_vector3d(0.0, 0.0, 1.0, normal)


    #@cython.cdivision(True)
    cpdef list propagate(self, Ray incident_ray, double ni, double nr):
        """
        Propagate a ray through an ideal focusing surface with focal length self.f.

        The surface is mathematically ideal: it redirects rays as if they were focused
        to (or came from) a focal point defined by the intersection of the principal
        ray with the focal plane.

        Parameters
        ----------
        incident_ray : Ray
            The incident ray in the coordinate system of the surface.
        ni : double
            Refractive index of the incident medium.
        nr : double
            Refractive index of the transmitted medium (used for transmitted ray).
        Returns
        -------
        list
            Rays resulting from transmission/reflection depending on reflectivity.
        """

        cdef Vector3d intersection_point, FP, d, dirv
        cdef double rz, scale, reflect

        # Compute intersection with the surface (XY plane in your case)
        self.intersection_cy(incident_ray, intersection_point)

        # If no intersection, return empty list (intersection_point set to NaN)
        if isnan(intersection_point(0)) or isnan(intersection_point(1)) or isnan(intersection_point(2)):
            return []

        dirv = incident_ray._direction
        rz = dirv(2)

        # Ray parallel to focal plane definition (old code used abs(rz) in denominator)
        if rz == 0.0:
            return []

        # Focal point as where the principal ray hits the focal plane:
        # FP = dir * f / abs(rz)
        scale = self.f / fabs(rz)
        FP = dirv * scale

        # Spectral / constant reflectivity handling (same convention as spherical example)
        reflect = self.reflectivity_function(self, incident_ray.wavelength)
        if reflect > 1.0 or reflect < 0.0:
            raise ValueError

        # Build rays: transmitted (ideal focusing) and/or reflected
        # Transmitted direction
        d = FP - intersection_point
        if self.f < 0.0:
            d = -d
        d.normalize()

        # Reflected direction
        cdef Vector3d rdir
        rdir = intersection_point - FP
        rdir.normalize()

        if (reflect == 0.0):
            # Pure transmission
            return [Ray.fast_init(intersection_point,
                                d,
                                incident_ray.intensity,
                                incident_ray.wavelength,
                                nr,
                                incident_ray.label,
                                incident_ray.draw_color,
                                None,
                                0.,
                                self.id,
                                0,
                                incident_ray._parent_cnt + 1)]

        elif (reflect == 1.0):
            # Pure reflection
            return [Ray.fast_init(intersection_point,
                                rdir,
                                incident_ray.intensity,
                                incident_ray.wavelength,
                                ni,
                                incident_ray.label,
                                incident_ray.draw_color,
                                None,
                                0.,
                                self.id,
                                0,
                                incident_ray._parent_cnt + 1)]
        else:
            # Beam-splitter behavior: split intensity
            return [Ray.fast_init(intersection_point,
                                d,
                                incident_ray.intensity * (1.0 - reflect),
                                incident_ray.wavelength,
                                nr,
                                incident_ray.label,
                                incident_ray.draw_color,
                                None,
                                0.,
                                self.id,
                                0,
                                incident_ray._parent_cnt + 1),
                    Ray.fast_init(intersection_point,
                                rdir,
                                incident_ray.intensity * reflect,
                                incident_ray.wavelength,
                                ni,
                                incident_ray.label,
                                incident_ray.draw_color,
                                None,
                                0.,
                                self.id,
                                0,
                                incident_ray._parent_cnt + 1)]

