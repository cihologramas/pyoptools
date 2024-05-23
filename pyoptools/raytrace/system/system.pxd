from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.misc.plist.plist cimport plist
from pyoptools.misc.picklable.picklable cimport Picklable

cdef class System(Picklable):
    cdef public plist _complist
    cdef public double n  # This could be changed to material
    cdef public list _np_rays
    cdef public list _p_rays
    cdef public int _max_ray_parent_cnt
    cdef public float intensity_threshold

    cdef public int _exit_status_flag # 0: no error, 1: error

    cpdef distance(self, Ray ri)
    cpdef propagate_ray(self, Ray ri)
    cpdef propagate_ray_ns(self, Ray gr, dpath)
