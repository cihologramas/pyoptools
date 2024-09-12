cdef class poly2d:
    cdef public int[::1] cohef

    cdef int[::1] px
    cdef int[::1] py
    cdef double[::1] coeff
    cdef public int order
    cdef int clen

    # cache of the derivatives
    cdef poly2d dx, dy

    # TODO: This name has to be changed to something more meaningfull
    cpdef double peval(self, double x, double y)
    
    cpdef eval(self, x, y)
    cdef double[:, :] eval2d(self, double[:, :] x, double[:, :] y)
    cdef double[:, :, :] eval3d(self, double[:, :, :] x, double[:, :, :] y)
    cdef double[:, :, :, :] eval4d(self, double[:, :, :, :] x, double[:, :, :, :] y)
    cdef double[:, :, :, :, :] eval5d(self, double[:, :, :, :, :] x,
                                      double[:, :, :, :, :] y)

cpdef int pxpy2i(int px, int py)
cpdef int ord2i(int o)
cpdef tuple[int, int] index_to_powers(int i)
cpdef tuple indices_to_powers(int[:] indices)

