from libc.math cimport sqrt
from pyoptools.misc.cmisc.eigen cimport VectorXd
import numpy as np

cdef class PolyRDeriv:
    """Class that represents the derivative of a radial polynomial.

    This class implements the calculation of the partial derivatives of a radial
    polynomial with respect to x or y coordinates. A radial polynomial is a
    polynomial expressed in terms of the radial coordinate r = sqrt(x^2 + y^2).

    Parameters
    ----------
    coeff : list
        List of coefficients of the radial polynomial.
        The i-th element corresponds to the coefficient of r^i.
    axis : int
        Axis along which to take the derivative.
        Must be 0 (for x) or 1 (for y).

    Attributes
    ----------
    order : int
        The order of the polynomial (number of coefficients - 1)
    axis : int
        Axis along which to take the derivative (0 for x, 1 for y)
    _coeff : VectorXd
        Eigen vector storing the polynomial coefficients
    """

    def __init__(self, coeff, axis):
        """Initialize a radial polynomial derivative.

        Parameters
        ----------
        coeff : list
            List of coefficients of the radial polynomial.
            The i-th element corresponds to the coefficient of r^i.
        axis : int
            Axis along which to take the derivative.
            Must be 0 (for x) or 1 (for y).

        Raises
        ------
        ValueError
            If the axis parameter is not 0 or 1.
        """
        self.order = len(coeff)

        self._coeff = VectorXd(self.order)

        for i in range(self.order):
            (<double*>(&(self._coeff(i))))[0] = coeff[i]

        if not (axis == 0 or axis == 1):
            raise ValueError(str(axis)+" is not a valid axis for derivative,"
                             " use 0 or 1 for x/y.")
        self.axis = axis

    cdef double eval_cy(self, double x, double y) noexcept nogil:
        """Evaluate the derivative of the radial polynomial at point (x, y).

        This is a C-level function that computes the value of the derivative
        at the given coordinates. The calculation is performed in a nogil context
        for better performance.

        Parameters
        ----------
        x : double
            X coordinate of the evaluation point
        y : double
            Y coordinate of the evaluation point
        Returns
        -------
        double
            Value of the polynomial derivative at point (x, y)

        Notes
        -----
        Returns 0 for points too close to the origin (r < 1e-10) to avoid
        numerical instability.
        """
        cdef int k
        cdef double result = 0
        cdef double r = sqrt(x**2 + y**2)

        if r < 1e-10:
            return 0

        for k in range(1, self.order):
            result += k*self._coeff(k) * r**(k-2)

        if self.axis == 0:
            result = result*x
        else:
            result = result*y
        return result

    def eval(self, double[:, ::1] x, double[:, ::1] y):
        """Evaluate the derivative of the radial polynomial over a grid of points.

        Parameters
        ----------
        x : ndarray of shape (n, m)
            Array of x coordinates where to evaluate the derivative.
            Must be a 2D array of doubles.
        y : ndarray of shape (n, m)
            Array of y coordinates where to evaluate the derivative.
            Must be a 2D array of doubles.

        Returns
        -------
        ndarray of shape (n, m)
            Array containing the values of the polynomial derivative
            evaluated at each point (x[i,j], y[i,j]).

        Notes
        -----
        This method calls the C-level eval_cy function for each point
        in the input arrays. Input arrays must have the same shape.
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
