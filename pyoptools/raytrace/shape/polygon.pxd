
from pyoptools.raytrace.shape.shape cimport Shape


cdef class Polygon(Shape):
    cdef public tuple coord
    cdef public int samples
    
