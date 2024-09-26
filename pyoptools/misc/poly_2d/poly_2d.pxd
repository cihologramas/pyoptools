from pyoptools.misc.cmisc.eigen cimport VectorXd, VectorXi, MatrixXd


cdef class Poly2D:

    cdef VectorXi px
    cdef VectorXi py
    cdef VectorXd _coeff
    cdef public int order
    cdef int _num_coeff

    # cache of the derivatives
    cdef Poly2D dx, dy

    cdef double eval_cy(self, double x, double y) noexcept nogil
    cdef void eval2d(self, MatrixXd& x, MatrixXd& y, MatrixXd& result) noexcept nogil

cpdef int pxpy2i(int px, int py)
cpdef int ord2i(int o)
cpdef tuple[int, int] index_to_powers(int i)
cpdef tuple indices_to_powers(int[:] indices)
