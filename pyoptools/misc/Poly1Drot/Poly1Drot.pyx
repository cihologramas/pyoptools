import numpy as np

# cimport cython

cdef extern from "math.h":
    double pow(double, double)

cdef class poly1DrotDeriv:
    """ Class defining the derivative of a rotationally symmetric 1D polynomial.
    Generally not instantiated indivually, but as part of poly1Drot class
    """

    cdef public object coef
    # cdef np.float64_t *coef_c
    cdef int clen

    cdef public int wrt

    def __init__(self, coef, wrt):
        """Initilizer for a rotationally symmetric 1D polynomial derivative

        Parameters
        ----------
        coef : tuple of float
            All the coefficients starting from index 0
        wrt: int
            Selector for the axis to find the derivative with respect to.
            Either 0 or 1 for the x or y axis respectively. Other values
            will raise a ValueError.
        """
        self.coef = np.array(coef, dtype=np.float64)

        self.coef_c = np.PyArray_DATA(self.coef)
        self.clen = np.uint32(len(self.coef))

        if not (wrt == 0 or wrt ==1):
            raise ValueError(wrt+" is not a valid axis for derivative,"
                                 " use 0 or 1 for x/y.")
        self.wrt = wrt

    cpdef double peval(self, double x, double y):
        """
        Evaluate the derivative at the point x, y

        Parameters
        ---
        x : float
            x value at which to evaluate
        y : float
            y value at which to evaluate
        """

        cdef double r, s, a
        cdef int i = 0

        r = pow(x**2 + y**2, 0.5)
        if r == 0:
            return 0

        s = 0
        for i in range(self.clen):
            a = self.coef_c[i]
            if self.wrt == 0:
                s += (i+1) * a * x * r**(i-1)
            elif self.wrt == 1:
                s += (i+1) * a * y * r**(i-1)
        return s

    def eval(self, x, y):
        raise NotImplementedError

cdef class poly1Drot:
    """ Class defining a rotationally symmetric 1D polynomial.
    """

    cdef public object coef

    # cdef np.float64_t *coef_c
    cdef int clen

    def __init__(self, coef):
        """Initilizer for a rotationally symmetric 1D polynomial
        This can be evaluated at any point x, y in the same way as a 2D
        polynomial.

        Parameters
        ----------
        coef : tuple of float
            All the coefficients starting from index 0
        """

        self.coef = np.array(coef, dtype=np.float64)
        self.coef_c = np.PyArray_DATA(self.coef)
        self.clen = np.uint32(len(self.coef))

    def eval(self, x, y):
        """Generic numpy evaluation method for the polynomial surface value
        """
        r = np.sqrt(x**2 + y**2)
        s = np.zeros(np.shape(r))
        for i, a in enumerate(self.coef):
            s += a * r**(i)
        return s

    cpdef double peval(self, double x, double y):
        """Cython implementation of polynomial evaluation.
        """
        cdef double r, s, a
        cdef int i = 0

        r = pow(x**2 + y**2, 0.5)
        s = 0
        for i in range(self.clen):
            a = self.coef_c[i]
            s += a * pow(r, i)
        return s

    def meval(self):
        raise ValueError

    def dxdy(self):
        """Method to return a tuple of poly1DrotDeriv objects, for derivatives
        along a x and y axis respectively.
        """

        dx = poly1DrotDeriv(self.coef, 0)
        dy = poly1DrotDeriv(self.coef, 1)

        return (dx, dy)
