from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray
cimport numpy as np
cdef class Plane(Surface):
    cpdef _intersection(self,Ray A)
    cpdef np.ndarray normal(self,ri)
