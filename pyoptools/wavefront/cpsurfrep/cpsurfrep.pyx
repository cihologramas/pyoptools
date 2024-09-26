from pyoptools.misc import *

import numpy as np

cimport cython

# from numpy.ma import is_masked, MaskedArray
cdef extern from "math.h":
    double sqrt(double) nogil
    double atan2(double, double) nogil


@cython.boundscheck(False)  # turn of bounds-checking for entire function
@cython.wraparound(False)
def cpw_evaluate_c(self, double [:] k, samples=(512, 512)):
    """Plane wave evaluate, return the 2d polynomials

    Arguments:

    k -- Propagation direction of the incident plane wave

    """
    cdef double r, iang
    cdef int i

    # knorm=k/sqrt(dot(k,k))

    r=sqrt(k[0]**2+k[1]**2)

    iang =atan2(r, k[2])

    cdef double [:] cf
    cdef double [:] ci

    cf=np.zeros((self.nc,))
    ci=np.zeros((self.nc,))

    cdef double [:] pof
    for i in range(self.nc):
        pof=self.pf[i]
        cf[i]=np.polyval(pof, iang)
        poi=self.pi[i]
        ci[i]=np.polyval(poi, iang)

    df=poly2d(cf)
    di=poly2d(ci)
    return df, di
