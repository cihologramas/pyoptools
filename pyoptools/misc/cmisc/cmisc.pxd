# cython: profile=True

cimport numpy as np


cdef extern from "Python.h":

        cdef void Py_INCREF( object )


cdef extern from "numpy/arrayobject.h":

    cdef void import_array()
    
    cdef object PyArray_ZEROS( int nd, np.npy_intp *dims, int typenum, int fortran )
    cdef object PyArray_SimpleNew( int nd, np.npy_intp *dims, int typenum )
    cdef object PyArray_EMPTY( int nd, np.npy_intp *dims, int typenum, int fortran )
    
    int PyArray_ISCARRAY( np.ndarray instance ) # I can't get this one to work?!?

    int PyArray_FLOAT
    int PyArray_DOUBLE
    

cdef inline np.ndarray norm_vect(np.ndarray v)
cdef inline double vector_length(np.ndarray v)
    
cdef inline np.ndarray rot_mat_i(np.ndarray r)
cdef inline np.ndarray rot_mat(np.ndarray r)
cdef np.ndarray mvdot(np.ndarray mat, np.ndarray vec)
cdef np.ndarray mvdot1(np.ndarray mat, np.ndarray vec)

#cdef rot_mat_i(np.ndarray[np.double_t, ndim=1] r):
#cdef rot_mat(np.ndarray[np.double_t, ndim=1] r):
#def unwrap(inph,in_p=(), double uv=2*np.pi, int nn =1):
#cdef mvdot(np.ndarray[np.double_t, ndim=2, mode="c"] mat, np.ndarray[np.double_t, ndim=1, mode="c"] vec):

cdef void mvdotf(np.float64_t *ret, np.float64_t * mat,np.float64_t * vec )



cdef inline np.ndarray empty_mat( int M, int N )
cdef inline np.ndarray empty_vec( int M )
cdef inline np.ndarray zero_mat( int M, int N )
cdef inline np.ndarray zero_vec( int M )
cdef inline void clear_mat( np.ndarray A )
cdef inline void clear_vec( np.ndarray x )
