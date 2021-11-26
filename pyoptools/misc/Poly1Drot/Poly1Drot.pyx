#cython: profile=True

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "math.h":
    double pow(double,double)

cdef class poly1DrotDeriv:

    #fn_selector = {
    #    'x' : lambda x, y, r, i, a : (i+1) * a * x * r**(i-1),
    #    'y' : lambda x, y, r, i, a : (i+1) * a * y * r**(i-1)
    #}

    cdef public object coef
    cdef public int wrt

    def __init__(self, coef, wrt):
        self.coef = np.array(coef, dtype=np.float64)
        if not (wrt == 0 or wrt ==1):
            raise ValueError(wrt+' is not a valid axis for derivative, use 0 or 1 for x/y.')
        self.wrt = wrt

    cpdef double peval(self, double x, double y):

        cdef double r, s, a
        cdef int i

        r = pow(x**2 + y**2, 0.5)
        if r == 0:
            return 0

        s = 0
        for i, a in enumerate(self.coef):
            if self.wrt == 0:
                s += (i+1) * a * x * r**(i-1)
            elif self.wrt == 1:
                s += (i+1) * a * y * r**(i-1)
        return s

    def eval(self, x, y):
        raise NotImplementedError

cdef class poly1Drot:

    cdef public object coef

    def __init__(self, coef):
        self.coef = np.array(coef, dtype=np.float64)

    def eval(self, x, y):
        r = np.sqrt(x**2 + y**2)
        s = np.zeros(np.shape(r))
        for i, a in enumerate(self.coef):
            s += a * r**(i)
        return s

    cpdef double peval(self, double x, double y):
        cdef float r, s
        cdef int i

        r = pow(x**2 + y**2, 0.5)
        s = 0
        for i, a in enumerate(self.coef):
            s += a * r**(i)
        return s

    def meval(self):
        raise ValueError

    def dxdy(self):

        dx = poly1DrotDeriv(self.coef, 0)
        dy = poly1DrotDeriv(self.coef, 1)

        return (dx, dy)
