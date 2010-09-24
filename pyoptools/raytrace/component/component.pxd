from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.plist.plist cimport plist
from pyoptools.misc.picklable.picklable cimport Picklable

cdef class Component(Picklable):
    cdef public plist _surflist
    cdef public object _material

    cpdef distance(self,Ray ri_)
