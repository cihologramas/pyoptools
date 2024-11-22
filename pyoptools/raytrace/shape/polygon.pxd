
from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.cmisc.eigen cimport Vector2d


cdef class Polygon(Shape):
    cdef Vector2d point_a, point_b, point_c
    cdef public int samples
