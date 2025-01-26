from pyoptools.misc.cmisc.eigen cimport VectorXd, MatrixXd
from .poly_r_deriv cimport PolyRDeriv
from ..function_2d cimport Function2D


cdef class PolyR(Function2D):
    cdef VectorXd _coeff

    cdef public int order

    # Cache of the derivatives
    cdef public PolyRDeriv dx, dy

    cdef double eval_cy(self, double x, double y) noexcept nogil
