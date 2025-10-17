from libc.math cimport sqrt
from pyoptools.misc.cmisc.eigen cimport VectorXd, convert_vectorXd_to_list
from .poly_r_deriv import PolyRDeriv

import numpy as np

cdef class PolyR:
    """A 2D polynomial class for radially symmetric functions.

    Represents polynomials that depend only on the radial distance r = x**2 + y**2
    from the origin. While the polynomial is defined in terms of r, it can be
    evaluated at any point (x,y) in Cartesian coordinates.

    The polynomial is expressed in terms of r as:
    z = c_0 + c_1*r + c_2*r^2 + c_3*r^3 + ... + c_n*r^n

    where r = x**2 + y**2

    Parameters
    ----------
    coeff : array_like
        List or array of polynomial coefficients in ascending order of degree.
        The length of coeff determines the order of the polynomial.

    Attributes
    ----------
    order : int
        The order (degree) of the polynomial.
    _coeff : VectorXd
        Internal Eigen vector storing the polynomial coefficients.
    dx : PolyRDeriv or None
        Derivative with respect to x, computed lazily.
    dy : PolyRDeriv or None
        Derivative with respect to y, computed lazily.
    """

    def __init__(self, coeff):
        self.order = len(coeff)
        self._coeff = VectorXd(self.order)
        for i in range(self.order):
            (<double*>(&(self._coeff(i))))[0] = coeff[i]
        self.dx = None
        self.dy = None

    cdef double eval_cy(self, double x, double y) noexcept nogil:
        """Evaluate the polynomial at point (x,y) in Cython.

        Parameters
        ----------
        x : double
            X coordinate.
        y : double
            Y coordinate.
        Returns
        -------
        double
            Value of the polynomial at point (x,y).
        """
        cdef int k
        cdef double result = 0
        cdef double r = sqrt(x**2 + y**2)
        for k in range(self.order):
            result += self._coeff(k) * r**k
        return result

    def eval(self, double[:, ::1] x, double[:, ::1] y):
        """Evaluate the polynomial at arrays of coordinates (x,y).

        Parameters
        ----------
        x : ndarray
            2D array of x coordinates, shape (n_rows, n_cols).
        y : ndarray
            2D array of y coordinates, shape (n_rows, n_cols).

        Returns
        -------
        ndarray
            2D array containing the polynomial values at each (x,y) point,
            shape matches input arrays.
        """
        cdef int n_rows = x.shape[0]
        cdef int n_cols = x.shape[1]
        z_array = np.zeros((n_rows, n_cols), dtype=np.float64)
        cdef double[:, ::1] z = z_array
        cdef int i, j
        for i in range(n_rows):
            for j in range(n_cols):
                z[i, j] = self.eval_cy(x[i, j], y[i, j])
        return z_array

    cpdef tuple[Function2D, Function2D] dxdy(self):
        """Compute the partial derivatives with respect to x and y.

        Returns
        -------
        tuple of PolyRDeriv
            A tuple (dx, dy) containing the derivative polynomials:
                - dx: derivative with respect to x
                - dy: derivative with respect to y
        """

        coeff_list = convert_vectorXd_to_list(self._coeff)
        dx = PolyRDeriv(coeff_list, 0)
        dy = PolyRDeriv(coeff_list, 1)
        return (dx, dy)
