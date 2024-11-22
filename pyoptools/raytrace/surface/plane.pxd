from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.misc.cmisc.eigen cimport Vector3d

cdef class Plane(Surface):
    cdef void _calculate_intersection(self, Ray, Vector3d&) noexcept nogil
    cdef void _calculate_normal(self, Vector3d&, Vector3d&) noexcept nogil
    cdef double topo_cy(self, double x, double y) noexcept nogil
