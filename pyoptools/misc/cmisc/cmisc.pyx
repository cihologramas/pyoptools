cimport cython

from cython.view cimport array


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

    cdef double norm, sq_norm

    # Calculate the norm (Euclidean length) of the vector for 3 elements
    sq_norm = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2])

    # Normalize the vector
    if sq_norm != 1.:
        norm = sqrt(sq_norm)
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

        \text{length} = \\sqrt{v[0]^2 + v[1]^2 + v[2]^2}

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
        A 2-dimensional contiguous memoryview representing a 3x3 inverse
        rotation matrix.

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

    This function initializes an uninitialized 1D array (vector) of length
    `length` with double-precision floating-point numbers. The array is
    created using Cython's memory allocation features and returned as a
    memoryview.

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
    cdef double[::1] vec = array(shape=(length,),
                                 itemsize=sizeof(double),
                                 format="d", mode="c")

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
    cdef double[::1] vec = array(shape=(length,),
                                 itemsize=sizeof(double),
                                 format="d", mode="c")

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
    else:
        raise ValueError("Unsupported number of dimensions: %d. Supported "
                         "dimensions are 2, 3, 4, or 5." % ndim)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef double[:, ::1] zero_memarray_2d(Py_ssize_t nx, Py_ssize_t ny):
    cdef double[:, ::1] arr = cython.view.array(shape=(nx, ny),
                                                itemsize=cython.sizeof(double),
                                                format="d", mode="c")
    cdef Py_ssize_t i, j
    for i in range(nx):
        for j in range(ny):
            arr[i, j] = 0.0
    return arr


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cdef double[:, :, ::1] zero_memarray_3d(Py_ssize_t nx, Py_ssize_t ny, Py_ssize_t nz):
    cdef double[:, :, ::1] arr = cython.view.array(shape=(nx, ny, nz),
                                                   itemsize=cython.sizeof(double),
                                                   format="d", mode="c")
    cdef Py_ssize_t i, j, k
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                arr[i, j, k] = 0.0
    return arr


@cython.boundscheck(False)  # Disable bounds checking for performance
@cython.wraparound(False)   # Disable negative indexing for performance
cdef double[::1]to_vector(object py_obj):
    """
    Convert a Python object (list, tuple, numpy array) to a memory view of
    double values.

    Parameters:
    py_obj (list or tuple or array): The input Python object to convert.

    Returns:
    memoryview: A Cython memory view of the input object.
    """
    cdef Py_ssize_t size
    cdef double[::1] memview

    # Get the size of the input
    size = len(py_obj)

    # Create a memory view array to hold the converted values
    memview = empty_vec(size)

    # Copy elements from the Python object to the memory view
    for i in range(size):
        memview[i] = py_obj[i]

    return memview


@cython.boundscheck(False)  # Disable bounds checking for performance
@cython.wraparound(False)   # Disable negative indexing for performance
cpdef double[::1] dot_product_3x3_matrix_vector(double[:, :] matrix, double[:] vector):
    """
    Compute the dot product of a 3x3 matrix and a 3-element vector.

    Parameters
    ----------
    matrix : memoryview of shape (3, 3) and dtype double
        A memory view representing a 3x3 matrix.

    vector : memoryview of shape (3,) and dtype double
        A memory view representing a 3-element vector.

    Returns
    -------
    memoryview
        A memory view representing the resulting 3-element vector after the dot product.

    """
    cdef int i, j
    cdef double[::1] result = cython.view.array(shape=(3,),
                                                itemsize=cython.sizeof(cython.double),
                                                format="d")

    # Compute the dot product
    for i in range(3):
        result[i] = 0.0
        for j in range(3):
            result[i] += matrix[i, j] * vector[j]

    return result


@cython.boundscheck(False)  # Disable bounds checking for performance
@cython.wraparound(False)   # Disable negative indexing for performance
cdef bint allclose_cython(double[::1] a, double[::1] b, double atol):
    """
    Simplified Cython version of numpy.allclose for 1D memory views of length 3.

    Parameters
    ----------
    a : memoryview of shape (3,) and dtype double
        First input memory view representing a 3-element vector.

    b : memoryview of shape (3,) and dtype double
        Second input memory view representing a 3-element vector.

    atol : double, optional
        Absolute tolerance. Default is 1e-08.

    Returns
    -------
    bool
        True if all corresponding elements of `a` and `b` are approximately
        equal within `atol`.
        Otherwise, False.
    """
    cdef int i

    for i in range(3):
        if fabs(a[i] - b[i]) > atol:
            return False  # Elements are not approximately equal
    return True  # All elements are approximately equal


@cython.boundscheck(False)  # Disable bounds checking for performance
@cython.wraparound(False)   # Disable negative indexing for performance
cdef double norm_3d_vector(double[::1] vec):
    """
    Calculate the Euclidean norm of a 3D vector.

    Parameters
    ----------
    vec : memoryview of shape (3,) and dtype double, contiguous
        The 3D vector for which to calculate the norm. Must be a contiguous memory view.

    Returns
    -------
    double
        The Euclidean norm (length) of the 3D vector.
    """
    # Ensure the memoryview is of the correct shape
    assert vec.shape[0] == 3, "Input vector must have 3 elements."

    # Calculate the Euclidean norm
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])
