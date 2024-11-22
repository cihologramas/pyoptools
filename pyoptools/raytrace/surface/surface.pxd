
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.picklable.picklable cimport Picklable
from pyoptools.misc.cmisc.eigen cimport Vector3d

ctypedef double (*ReflectivityFunctionPtr)(Surface, double)


cdef class Surface(Picklable):
    # Object representing the surface reflectivity
    cdef public double reflectivity

    # Cutoff values for filters, the values are given in um as the wavelenghts
    # for the rays

    cdef double cutoff, lower_cutoff, upper_cutoff

    # Pointer to the reflectivity function to use

    cdef ReflectivityFunctionPtr reflectivity_function

    # Object containing the Shape limiting the surface
    cdef public Shape shape

    # List of rays that have hit the surface
    cdef public list _hit_list

    # Surface identifiers or tags (clarify its purpose if needed)
    cdef public list id

    # Method representing the topography of the surface.
    # Used for rendering the surface in Jupyter Lab
    # Also to find the bounding box when using iterative methods to find the
    # ray intersection.
    cpdef list topo(self, list x, list y)

    # fast cython version of topo, for one point (x,y). Must be overloaded in all
    # subclasses of Surface
    cdef double topo_cy(self, double x, double y) noexcept nogil

    # Cython method to calculate the surface intersection with a ray
    # To be used by pyoptools Cython functions and methods
    # TODO: Check what needs to be done for this method to be nogil
    cdef void intersection_cy(self,
                              Ray incident_ray,
                              Vector3d& intersection_point)  # noexcept nogil

    # Cython method to calculate the normal at a given intersection point
    # To be used by pyoptools Cython functions and methods
    cdef void normal_cy(self,
                        Vector3d& intersection_point,
                        Vector3d &Normal) noexcept nogil

    # Method to be overloaded by subclasses that define the specific surface geometry
    cdef void _calculate_intersection(self, Ray, Vector3d&) noexcept nogil

    # Method to be overloaded by subclasses that define the specific surface geometry
    cdef void _calculate_normal(self, Vector3d&, Vector3d&) noexcept nogil

    # Method to calculate the distance from the ray to the surface intersection
    cdef double distance_cy(self, Ray incident_ray, Vector3d& intersection_point)

    # Method to clear all internal structures and data within a surface
    # after propagation. Prepares the surface to be ready for a new propagation
    cpdef reset(self)

    # Calculates the propagation of a ray through the surface
    cpdef list propagate(self, Ray ri, double ni, double nr)

    # Methods to be able to simulate optical filters. These are just for
    # internal use in the class.
    cdef double constant_reflectivity(self, double wavelength) noexcept nogil
    cdef double longpass_reflectivity(self, double wavelength) noexcept nogil
    cdef double shortpass_reflectivity(self, double wavelength) noexcept nogil
    cdef double bandpass_reflectivity(self, double wavelength) noexcept nogil

    # Legacy methods that need to be reviewed
    cpdef pw_propagate1(self, Ray ri, ni, nr, rsamples, isamples, knots)
    cpdef pw_propagate(self, Ray ri, ni, nr, rsamples, shape, order, z)
    cpdef pw_propagate_list(self, Ray ri, ni, nr, rsamples, z)
    cpdef wf_propagate(self, wf, ni, nr, samples, shape, knots)
    cpdef pw_cohef(self, ni, nr, ilimit, slimit, step, order, rsamples, zb)
