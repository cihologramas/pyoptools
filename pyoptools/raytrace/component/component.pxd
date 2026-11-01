from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.plist.plist cimport plist
from pyoptools.misc.picklable.picklable cimport Picklable

cdef class Component(Picklable):
    cdef public plist _surflist
    cdef public object _material
    cdef public object render_color, render_opacity, render_smooth, render_thickness

    cpdef distance(self, Ray ri_)
