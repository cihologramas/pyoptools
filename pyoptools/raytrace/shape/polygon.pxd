from libcpp.vector cimport vector
from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.misc.cmisc.eigen cimport Vector2d, Vector3d


cdef class Polygon(Shape):
    cdef vector[Vector2d] poly_points
    cdef public int samples
    cdef tuple _coord
