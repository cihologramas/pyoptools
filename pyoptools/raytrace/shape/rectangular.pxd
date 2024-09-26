
from pyoptools.raytrace.shape.shape cimport Shape


cdef class Rectangular(Shape):
    cdef public tuple[double, double] size, offset
    cdef public tuple[int, int] samples
