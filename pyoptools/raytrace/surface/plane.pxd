from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.misc.cmisc.linalg cimport Vector3

cdef class Plane(Surface):
    cdef void _calculate_intersection(self, Ray, Vector3*) noexcept nogil
    cdef void _calculate_normal(self, Vector3*, Vector3*) noexcept nogil
