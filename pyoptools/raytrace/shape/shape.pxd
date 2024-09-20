from pyoptools.misc.picklable.picklable cimport Picklable
from pyoptools.misc.cmisc.eigen cimport Vector3d

cdef class Shape(Picklable):
    cdef bint hit_cy(self, Vector3d& point)  noexcept nogil
    cpdef pointlist(self)
    cpdef limits(self)
