cimport cython

from cython.view cimport array
from libc.stdlib cimport malloc, free


cdef extern from "math.h":
    double fabs(double) nogil
    double ceil(double) nogil
    bint isnan(double x) nogil
    double sin(double) nogil
    double cos(double) nogil
    double sqrt(double) nogil


@cython.boundscheck(False)   # Disable bounds checking for array accesses
@cython.wraparound(False)    # Disable negative index wraparound
cdef inline double[::1] norm_vect(double[::1] v):
    """
    Normalize a vector of fixed length (3 elements).

    Parameters
    ----------
    v : double[::1]
        A 1-dimensional memoryview of type `double` with shape (3,) representing
        the vector to be normalized.

    Returns
    -------
    double[::1]
        A memoryview of the normalized vector (of the same shape as input).

    Notes
    -----
    This function normalizes the vector `v` in-place, modifying the input.

    The input vector `v` is expected to have exactly 3 elements.
    """

    assert v.shape[0] == 3, "Input vector must have exactly 3 elements."

    cdef double norm

    # Calculate the norm (Euclidean length) of the vector for 3 elements
    norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

    # Normalize the vector
    v[0] /= norm
    v[1] /= norm
    v[2] /= norm

    # Return the modified vector as a memoryview
    return v

@cython.boundscheck(False)   # Disable bounds checking for array accesses
@cython.wraparound(False)    # Disable negative index wraparound
cdef inline double vector_length(double[::1] v):
    """
    Compute the length (magnitude) of a vector.

    This function calculates the Euclidean norm (magnitude) of a 3-dimensional
    vector provided as a memoryview.

    Parameters
    ----------
    v : double[::1]
        A memoryview of type `double` with shape (3,) representing the vector
        whose length is to be computed.

    Returns
    -------
    double
        The length of the vector `v`.

    Notes
    -----
    The function computes the length of the vector using the formula:

    .. math::

        \text{length} = \sqrt{v[0]^2 + v[1]^2 + v[2]^2}

    The input vector `v` is expected to have exactly 3 elements. This function
    is implemented in Cython for efficient computation.
    """
    assert v.shape[0] == 3, "Input vector must have exactly 3 elements."

    cdef double length

    # Calculate the Euclidean norm of the vector
    length = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

    return length

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double[:, ::1] rot_mat_i(double[::1] r):
    """
    Returns the inverse transformation matrix for a rotation around the Z, Y, X axes.

    Parameters
    ----------
    r : double[::1]
        A 1-dimensional contiguous memoryview of type `double` with shape (3,)
        containing the rotation angles `(rx, ry, rz)` in radians for the rotation
        around the X, Y, and Z axes.

    Returns
    -------
    double[:, ::1]
        A 2-dimensional contiguous memoryview representing a 3x3 inverse rotation matrix.

    Notes
    -----
    The function calculates the inverse rotation matrix using trigonometric functions.
    """
    # Ensure the input has exactly 3 elements
    assert r.shape[0] == 3, "Input array 'r' must have exactly 3 elements."

    # Pre-calculate trigonometric functions
    cdef double cos_rx = cos(r[0])
    cdef double cos_ry = cos(r[1])
    cdef double cos_rz = cos(r[2])
    
    cdef double sin_rx = sin(r[0])
    cdef double sin_ry = sin(r[1])
    cdef double sin_rz = sin(r[2])

    # Allocate a contiguous memoryview on the heap
    cdef double[:, ::1] rv = \
        cython.view.array(shape=(3, 3), itemsize=sizeof(double), format="d", 
                          mode="c")

    rv[0, 0] = cos_ry * cos_rz
    rv[0, 1] = cos_ry * sin_rz
    rv[0, 2] = -sin_ry

    rv[1, 0] = cos_rz * sin_rx * sin_ry - cos_rx * sin_rz
    rv[1, 1] = sin_rx * sin_ry * sin_rz + cos_rx * cos_rz
    rv[1, 2] = cos_ry * sin_rx

    rv[2, 0] = sin_rx * sin_rz + cos_rx * cos_rz * sin_ry
    rv[2, 1] = cos_rx * sin_ry * sin_rz - cos_rz * sin_rx
    rv[2, 2] = cos_rx * cos_ry

    return rv  # Cast to a contiguous memoryview

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double[:, ::1] rot_mat(double[::1] r):
    """
    Returns the transformation matrix for a rotation around the Z, Y, X axes.

    The rotation is performed first around the Z axis, then around the Y axis,
    and finally around the X axis.

    Parameters
    ----------
    r : double[::1]
        A contiguous 1-dimensional memoryview of type `double` with shape (3,)
        containing the rotation angles `(rx, ry, rz)` in radians.

    Returns
    -------
    double[:, ::1]
        A 2-dimensional contiguous memoryview representing a 3x3 rotation matrix.

    Notes
    -----
    The function computes the rotation matrix using trigonometric functions.
    """
    # Ensure the input has exactly 3 elements
    assert r.shape[0] == 3, "Input array 'r' must have exactly 3 elements."

    # Pre-calculate trigonometric functions
    cdef double cos_rx = cos(r[0])
    cdef double cos_ry = cos(r[1])
    cdef double cos_rz = cos(r[2])

    cdef double sin_rx = sin(r[0])
    cdef double sin_ry = sin(r[1])
    cdef double sin_rz = sin(r[2])

    # Allocate a contiguous memoryview for the rotation matrix
     # Allocate a contiguous memoryview on the heap
    cdef double[:, ::1] rv = \
        cython.view.array(shape=(3, 3), itemsize=sizeof(double), format="d",
                          mode="c")

    # Compute the rotation matrix
    rv[0, 0] = cos_ry * cos_rz
    rv[0, 1] = cos_rz * sin_rx * sin_ry - cos_rx * sin_rz
    rv[0, 2] = sin_rx * sin_rz + cos_rx * cos_rz * sin_ry

    rv[1, 0] = cos_ry * sin_rz
    rv[1, 1] = sin_rx * sin_ry * sin_rz + cos_rx * cos_rz
    rv[1, 2] = cos_rx * sin_ry * sin_rz - cos_rz * sin_rx

    rv[2, 0] = -sin_ry
    rv[2, 1] = cos_ry * sin_rx
    rv[2, 2] = cos_rx * cos_ry

    # Cast to a contiguous memoryview and return
    return rv

cdef double[::1] empty_vec(int length):
        """
        Create an empty 1D array of doubles using Cython memoryviews.

        This function initializes an uninitialized 1D array (vector) of length `length` 
        with double-precision floating-point numbers. The array is created using 
        Cython's memory allocation features and returned as a memoryview.

        Parameters
        ----------
        length : int
            The length of the array (number of elements) to create.

        Returns
        -------
        vec : double[::1]
            A C-contiguous memoryview of the uninitialized array of doubles.

        Notes
        -----
        The array is created uninitialized, meaning it will contain arbitrary values.
        Use this function when you intend to fill the array immediately afterward, 
        as accessing uninitialized data can lead to undefined behavior.

        Examples
        --------
        >>> cdef double[::1] my_array = empty_vec(10)
        >>> print(len(my_array))
        10

        """
        # Allocate memory for 'length' elements of type double
        cdef double[::1] vec = array(shape=(length,), itemsize=sizeof(double), format="d", mode="c")
        
        return vec


cdef double[::1] zero_vec(int length):
    """
    Create a zero-initialized memoryview of doubles.

    Parameters
    ----------
    length : int
        The number of elements in the vector.

    Returns
    -------
    double[::1]
        A 1-dimensional contiguous memoryview of type double.
    """
    cdef int i
    # Allocate memory for 'length' elements of type double
    cdef double[::1] vec = array(shape=(length,), itemsize=sizeof(double), format="d", mode="c")

    # Zero-initialize the allocated memory
    for i in range(length):
        vec[i] = 0.

    return vec


cpdef zero_memarray(tuple shape):
    """
    Create a zero-initialized memoryview of doubles with the given shape.

    Parameters
    ----------
    shape : tuple
        The shape of the memoryview to be created. It determines the number
        of dimensions.

    Returns
    -------
    double[:, :], double[:, :, :], double[:, :, :, :], or double[:, :, :, :, :]
        A memoryview of type double initialized to zero with the given shape.
    """
    cdef int ndim = len(shape)

    if ndim == 2:
        return zero_memarray_2d(shape[0], shape[1])
    elif ndim == 3:
        return zero_memarray_3d(shape[0], shape[1], shape[2])
    elif ndim == 4:
        return zero_memarray_4d(shape[0], shape[1], shape[2], shape[3])
    elif ndim == 5:
        return zero_memarray_5d(shape[0], shape[1], shape[2], shape[3], shape[4])
    else:
        raise ValueError("Unsupported number of dimensions: %d. Supported dimensions are 2, 3, 4, or 5." % ndim)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef double[:, :] zero_memarray_2d(Py_ssize_t nx, Py_ssize_t ny):
    cdef double[:, :] arr = cython.view.array(shape=(nx, ny), itemsize=cython.sizeof(double), format="d", mode="c")
    cdef Py_ssize_t i, j
    for i in range(nx):
        for j in range(ny):
            arr[i, j] = 0.0
    return arr


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef double[:, :, :] zero_memarray_3d(Py_ssize_t nx, Py_ssize_t ny, Py_ssize_t nz):
    cdef double[:, :, :] arr = cython.view.array(shape=(nx, ny, nz), itemsize=cython.sizeof(double), format="d", mode="c")
    cdef Py_ssize_t i, j, k
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                arr[i, j, k] = 0.0
    return arr


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef double[:, :, :, :] zero_memarray_4d(Py_ssize_t nx, Py_ssize_t ny, Py_ssize_t nz, Py_ssize_t nw):
    cdef double[:, :, :, :] arr = cython.view.array(shape=(nx, ny, nz, nw), itemsize=cython.sizeof(double), format="d", mode="c")
    cdef Py_ssize_t i, j, k, l
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for l in range(nw):
                    arr[i, j, k, l] = 0.0
    return arr


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef double[:, :, :, :, :] zero_memarray_5d(Py_ssize_t nx, Py_ssize_t ny, Py_ssize_t nz, Py_ssize_t nw, Py_ssize_t nv):
    cdef double[:, :, :, :, :] arr = cython.view.array(shape=(nx, ny, nz, nw, nv), itemsize=cython.sizeof(double), format="d", mode="c")
    cdef Py_ssize_t i, j, k, l, m
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for l in range(nw):
                    for m in range(nv):
                        arr[i, j, k, l, m] = 0.0
    return arr
