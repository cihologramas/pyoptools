cimport numpy as np

cdef class poly2d:
    cdef public object cohef
    # Slow lists to be used for python and for garbage colection
    # if the slow lists are not attributes of the class they can get destroyed
    # so px_c and py_c will get erased
    cdef object px64,py64,px,py
    #internal Fast C list
    cdef np.float64_t *px_c
    cdef np.float64_t *py_c
    cdef np.float64_t *cohef_c
    cdef public int order
    cdef int clen  
    cdef object ctx0,prg0,prg1,queue0
    
    #cache of the derivatives
    cdef poly2d dx,dy
    
    cpdef double peval(self, double x, double y)
    cpdef eval(self, x, y)
    cpdef eval_2(self, x, y, key)

cpdef int pxpy2i(int px,int py)
cpdef ord2i(o)
cpdef i2pxpy(i)
