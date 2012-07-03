from pyoptools.raytrace.shape.shape cimport Shape

cdef class Circular(Shape):
    cdef public double radius
    cdef public tuple samples
