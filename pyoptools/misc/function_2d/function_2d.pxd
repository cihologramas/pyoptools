from pyoptools.misc.cmisc.eigen cimport VectorXd, VectorXi, MatrixXd


cdef class Function2D:
    """
    Base class declaration for 2D functions defined in Cartesian coordinates.
    This file is for importing the `Function2D` API in other `.pyx` files.
    """
    cdef double eval_cy(self, double x, double y) noexcept nogil
    cdef void eval2d(self, MatrixXd& x, MatrixXd& y, MatrixXd& result) noexcept nogil
    # def eval(self, double[:, ::1] x, double[:, ::1] y)
    cpdef tuple[Function2D, Function2D] dxdy(self)
