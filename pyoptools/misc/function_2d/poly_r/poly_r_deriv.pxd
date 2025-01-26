from pyoptools.misc.cmisc.eigen cimport VectorXd
from ..function_2d cimport Function2D


cdef class PolyRDeriv(Function2D):
    cdef VectorXd _coeff
    cdef int order
    cdef int axis

    cdef double eval_cy(self, double x, double y) noexcept nogil
