#!/usr/bin/env python
# cython: profile=True

import numpy as np
cimport numpy as np

import_array()

cimport cython


#from numpy.ma import is_masked, MaskedArray

cdef extern from "math.h":
    double fabs(double) nogil
    double ceil(double) nogil
    bint isnan(double x) nogil
    double sin(double) nogil
    double cos(double) nogil
    double sqrt(double) nogil

DTYPE = np.double
ctypedef np.double_t DTYPE_t


cdef inline np.ndarray norm_vect(np.ndarray v):
    '''Normalize a vector

    Normalizes a vector, and return a pointer to itself
    '''

    cdef np.float64_t * vd = <np.float64_t*>(np.PyArray_DATA(v))
    cdef double norm = sqrt(vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2])
    vd[0] = vd[0]/norm
    vd[1] = vd[1]/norm
    vd[2] = vd[2]/norm
    return v

cdef inline double vector_length(np.ndarray v):
    cdef np.float64_t * vd = <np.float64_t*>(np.PyArray_DATA(v))
    return sqrt(vd[0]*vd[0]+vd[1]*vd[1]+vd[2]*vd[2])


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.ndarray rot_mat_i(np.ndarray r):
    '''Returns the inverse transformation matrix for a rotation around the Z,Y,X axes

    Parameters

     r= (rx,ry,rz)
     '''

    # c=np.cos(r)
    # s=np.sin(r)
    cdef np.float64_t * rd = <np.float64_t*>(np.PyArray_DATA(r))

    cdef double c0 = cos(rd[0])
    cdef double c1 = cos(rd[1])
    cdef double c2 = cos(rd[2])

    cdef double s0 = sin(rd[0])
    cdef double s1 = sin(rd[1])
    cdef double s2 = sin(rd[2])

    # Slow Way
    # ~ rx=np.array([[ 1., 0., 0.],
    #~ [ 0., c0, s0],
    # ~ [ 0.,-s0, c0]])

    # ~ ry=np.array([[ c1, 0.,-s1],
    #~ [ 0., 1., 0.],
    # ~ [ s1, 0., c1]])

    # ~ rz=np.array([[ c2, s2, 0.],
    #~ [-s2, c2, 0.],
    # ~ [ 0., 0., 1.]])

    # ~ return np.dot(rx,np.dot(ry,rz))

    # Faster way
    # return np.array([[c1*c2,c1*s2,-s1],
    #                [c2*s0*s1-c0*s2,s0*s1*s2+c0*c2,c1*s0],
    #                [s0*s2+c0*c2*s1,c0*s1*s2-c2*s0,c0*c1]])

    # Fastest until now but ugly
    # TODO: Find a way to make this less ugly
    # cdef np.ndarray[np.double_t, ndim=2, mode="c"] rv=np.empty((3,3),dtype=np.double)
    cdef np.ndarray[np.double_t, ndim = 2, mode = "c"] rv = empty_mat(3, 3)
    rv[0, 0] = c1*c2
    rv[0, 1] = c1*s2
    rv[0, 2] = -s1

    rv[1, 0] = c2*s0*s1-c0*s2
    rv[1, 1] = s0*s1*s2+c0*c2
    rv[1, 2] = c1*s0

    rv[2, 0] = s0*s2+c0*c2*s1
    rv[2, 1] = c0*s1*s2-c2*s0
    rv[2, 2] = c0*c1

    return rv


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.ndarray rot_mat(np.ndarray r):
    '''Returns the transformation matrix for a rotation around the Z,Y,X axes

    The rotation is made first around the Z axis, then around the Y axis, and
    finally around the X axis.

    Parameters

    r= (rx,ry,rz)
    '''
    cdef np.float64_t * rd = <np.float64_t*>(np.PyArray_DATA(r))
    # c=cos(r)
    # s=sin(r)

    cdef double c0 = cos(rd[0])
    cdef double c1 = cos(rd[1])
    cdef double c2 = cos(rd[2])

    cdef double s0 = sin(rd[0])
    cdef double s1 = sin(rd[1])
    cdef double s2 = sin(rd[2])

    # ~ rx=array([[1. , 0., 0.],
    #~ [0. , c0,-s0],
    # ~ [0. , s0, c0]])
    # ~ ry=array([[ c1, 0., s1],
    #~ [ 0., 1., 0.],
    # ~ [-s1, 0., c1]])
    # ~ rz=array([[ c2,-s2, 0.],
    #~ [ s2, c2, 0.],
    # ~ [ 0., 0., 1.]])
    # ~ tm=dot(rz,dot(ry,rx))
    # ~ return tm
    # cdef np.ndarray[np.double_t, ndim=2, mode="c"] rv=np.empty((3,3),dtype=np.double)
    cdef np.ndarray[np.double_t, ndim = 2, mode = "c"] rv = empty_mat(3, 3)
    rv[0, 0] = c1*c2
    rv[0, 1] = c2*s0*s1-c0*s2
    rv[0, 2] = s0*s2+c0*c2*s1

    rv[1, 0] = c1*s2
    rv[1, 1] = s0*s1*s2+c0*c2
    rv[1, 2] = c0*s1*s2-c2*s0

    rv[2, 0] = -s1
    rv[2, 1] = c1*s0
    rv[2, 2] = c0*c1

    return rv


@cython.boundscheck(False)
@cython.wraparound(False)
def unwrap(inph, in_p=(), double uv=2*np.pi, int nn=1):
    """Return the input matrix unwrapped using the value given in uv.
    nn indicates how many nearest neighbours are going to be used in the unwrap
    procedure.

    The same as unwrapv, but using for-s
    """

    cdef np.ndarray[np.double_t, ndim = 2, mode = "c"] faseo
    if np.ma.isMaskedArray(inph) == False:
        faseo = inph.copy()
    else:
        faseo = inph.copy()
        faseo = np.where(np.ma.getmaskarray(inph), np.nan,  faseo)

    cdef int nx, ny
    nx, ny = (faseo.shape[0], faseo.shape[1])

    # If the initial unwraping point is not given, take the center of the image
    # as initial coordinate
    if in_p == ():
        in_p = (int(nx/2), int(ny/2))

    # Create a temporal space to mark if the points are already unwrapped
    # 0 the point has not been unwrapped
    # 1 the point has not been unwrapped, but it is in the unwrapping list
    # 2 the point was already unwrapped

    cdef np.ndarray[np.int_t, ndim = 2,  mode = "c"] fl = np.zeros((nx, ny), dtype=np.int)

    # List containing the points to unwrap

    l_un = [in_p]
    fl[in_p] = 1

    cdef int cx, cy
    cdef int nv, wvi
    cdef int i, j
    cdef int cx0, cx1, cy0, cy1
    cdef double wf, wv
    while len(l_un) > 0:
        # remove the first value from the list
        cx, cy = l_un.pop(0)

        # Put the coordinates of unwrapped the neighbors in the list
        # And check for wrapping
        nv = 0
        wv = 0
        cx0 = cx-nn
        cx1 = cx+nn+1
        cy0 = cy-nn
        cy1 = cy+nn+1
        if cx0 < 0:
            cx0 = 0
        if cy0 < 0:
            cy0 = 0

        if cx1 > nx:
            cx1 = nx
        if cy1 > ny:
            cy1 = ny

        for j in range(cy0, cy1):
            for i in range(cx0, cx1):

                if (fl[i, j] == 0) & (~isnan(faseo[i, j])):
                    fl[i, j] = 1
                    l_un.append((i, j))
                elif fl[i, j] == 2:
                    wf = (faseo[i, j]-faseo[cx, cy])/uv
                    wv = wv+wf
                    nv = nv+1
        if nv != 0:
            if wv > 0:
                wv = <int > (0.5+wv/nv)
            else:
                wv = <int > (-0.5+wv/nv)

        fl[cx, cy] = 2
        faseo[cx, cy] = faseo[cx, cy]+wv*uv

    return np.ma.masked_equal(faseo, np.nan)


@cython.boundscheck(False)
@cython.wraparound(False)
# cdef mvdot(np.ndarray[np.double_t, ndim=2, mode="c"] mat, np.ndarray[np.double_t, ndim=1, mode="c"] vec)
cdef np.ndarray mvdot(np.ndarray mat, np.ndarray vec):
    """
    Dot product between a 3x3 matrix and a 3 element vector
    """
    cdef np.ndarray [np.double_t, ndim= 1, mode = "c"] ret = zero_vec(3)
    cdef int i, j
    for j in range(3):
        for i in range(3):
            ret[j] = ret[j]+mat[j, i]*vec[i]
    return ret


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray mvdot1(np.ndarray mat, np.ndarray vec):
    """ 
    Dot product between a 3x3 matrix and a 3 element vector
    """
    cdef np.ndarray ret = zero_vec(3)
    cdef int i, j

    cdef np.float64_t * retp = <np.float64_t * >np.PyArray_DATA(ret)
    cdef np.float64_t * matp = <np.float64_t * >np.PyArray_DATA(mat)
    cdef np.float64_t * vecp = <np.float64_t * >np.PyArray_DATA(vec)

    for j in range(3):
        for i in range(3):
            retp[j] = retp[j]+matp[j*3+i]*vecp[i]

    return ret


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void mvdotf(np.float64_t * ret, np.float64_t * mat, np.float64_t * vec):
    cdef int i, j
    for j in range(3):
        for i in range(3):
            ret[j] = ret[j]+mat[j*3+i]*vec[i]


def dot_test(mat, vec):
    cdef np.ndarray[np.double_t, ndim= 1, mode = "c"] res = np.empty((3,), dtype=np.float64)
    mvdotf( < np.float64_t*>np.PyArray_DATA(res), 
           < np.float64_t*>np.PyArray_DATA(mat), < np.float64_t*>np.PyArray_DATA(vec))

    return res


def test_1(matrix):
    cdef np.float64_t * c = <np.float64_t*>np.PyArray_DATA(matrix)
    for i in range(9):
        print c[i]


def test_2(vector):
    cdef np.float64_t * c = <np.float64_t*>np.PyArray_DATA(vector)
    for i in range(3):
        print c[i]


#########################################################################
#
# Auxiliary functions taken from tokyo.pyx
# http://www.vetta.org/2009/09/tokyo-a-cython-blas-wrapper-for-fast-matrix-math/
#
#########################################################################


# Create a new empty double precision matrix
cdef inline np.ndarray empty_mat(int M, int N):
    cdef np.npy_intp length[2]
    length[0] = M
    length[1] = N
    Py_INCREF(np.NPY_DOUBLE)  # This is apparently necessary
    return PyArray_EMPTY(2, length, np.NPY_DOUBLE, 0)


# Create a new empty double precision vector
cdef inline np.ndarray empty_vec(int M):
    cdef np.npy_intp length[1]
    length[0] = M
    Py_INCREF(np.NPY_DOUBLE)  # This is apparently necessary
    return PyArray_EMPTY(1, length, np.NPY_DOUBLE, 0)


# Create a new zeroed double precision matrix
cdef inline np.ndarray zero_mat(int M, int N):
    cdef np.npy_intp length[2]
    length[0] = M
    length[1] = N
    Py_INCREF(np.NPY_DOUBLE)  # This is apparently necessary
    return PyArray_ZEROS(2, length, np.NPY_DOUBLE, 0)


# Create a new zeroed double precision vector
cdef inline np.ndarray zero_vec(int M):
    cdef np.npy_intp length[1]
    length[0] = M
    Py_INCREF(np.NPY_DOUBLE)  # This is apparently necessary
    return PyArray_ZEROS(1, length, np.NPY_DOUBLE, 0)


# Set a matrix to all zeros: must be doubles in contiguous memory.
cdef inline void clear_mat(np.ndarray A):

    if A.ndim != 2:
        raise ValueError("A is not a matrix")
    if A.descr.type_num != PyArray_DOUBLE:
        raise ValueError("A is not of type double")

    cdef double * ptr = <double*>A.data
    cdef unsigned int i
    for i in range(A.shape[0]*A.shape[1]):
        ptr[0] = 0.0
        ptr += 1


# Set a vector to all zeros: ust be doubles in contiguous memory.
cdef inline void clear_vec(np.ndarray x):

    if x.ndim != 1:
        raise ValueError("A is not a vector")
    if x.descr.type_num != PyArray_DOUBLE:
        raise ValueError("x is not of type double")

    cdef double * ptr = <double*>x.data
    cdef unsigned int i
    for i in range(x.shape[0]):
        ptr[0] = 0.0
        ptr += 1
