from pyoptools.misc.cmisc.eigen cimport VectorXd, VectorXi, MatrixXd
from ..function_2d cimport Function2D

cdef class Poly2D(Function2D):
    """
    Represents a 2D polynomial of the form:
      f(x, y) = Σ c[i] * x^(px[i]) * y^(py[i])
    where px and py are the powers of x and y respectively for each term.
    """

    # Polynomial coefficients and powers
    cdef VectorXi px          # Power of x for each term
    cdef VectorXi py          # Power of y for each term
    cdef VectorXd _coeff      # Coefficients for each term
    cdef public int order     # Order (degree) of the polynomial
    cdef int _num_coeff       # Number of coefficients/terms in the polynomial

    # Cached derivatives
    cdef Poly2D dx            # Cached derivative with respect to x (∂f/∂x)
    cdef Poly2D dy            # Cached derivative with respect to y (∂f/∂y)

    # Evaluation methods
    cdef double eval_cy(self, double x, double y) noexcept nogil

    cpdef tuple[Function2D, Function2D] dxdy(self)

# Utility functions
cpdef int pxpy2i(int px, int py)
cpdef int ord2i(int o)
cpdef tuple[int, int] index_to_powers(int i)
cpdef tuple indices_to_powers(int[:] indices)
