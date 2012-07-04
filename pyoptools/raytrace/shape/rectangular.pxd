
from pyoptools.raytrace.shape.shape cimport Shape


cdef class Rectangular(Shape):
    cdef public tuple size,samples,offset
    
