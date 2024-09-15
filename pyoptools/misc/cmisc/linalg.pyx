from cython cimport inline
from .linalg cimport Vector3, Matrix3x3

from libc.math cimport fabs, sin, cos, sqrt, nan


cdef inline Vector3 create_zero_vector3() noexcept nogil:
    """
    Create a new Vector3 initialized to zero.

    Returns
    -------
    Vector3
        A new Vector3 instance with all components set to zero.
    """
    cdef Vector3 result
    set_zero_vector3(&result)  # Pass the address of result to set_zero_vector3
    return result

cdef inline void set_zero_vector3(Vector3* v) noexcept nogil:
    """
    Set all components of the Vector3 to zero.

    Parameters
    ----------
    v : Vector3*
        Pointer to the Vector3 to be zeroed.
    """
    v[0].data[0] = 0.0
    v[0].data[1] = 0.0
    v[0].data[2] = 0.0

cdef inline void set_nan_vector3(Vector3* v) noexcept nogil:
    """
    Set all components of the Vector3 to nan.

    Parameters
    ----------
    v : Vector3*
        Pointer to the Vector3 to be zeroed.
    """
    v[0].data[0] = nan("")
    v[0].data[1] = nan("")
    v[0].data[2] = nan("")

cdef inline void add_vector3(const Vector3* a, const Vector3* b, Vector3* result) noexcept nogil:
    """
    Add two Vector3 instances and store the result in a third Vector3.

    Parameters
    ----------
    a : const Vector3*
        Pointer to the first Vector3.
    b : const Vector3*
        Pointer to the second Vector3.
    result : Vector3*
        Pointer to the Vector3 where the result will be stored.
    """
    result[0].data[0] = a[0].data[0] + b[0].data[0]
    result[0].data[1] = a[0].data[1] + b[0].data[1]
    result[0].data[2] = a[0].data[2] + b[0].data[2]

cdef inline void substract_vector3(const Vector3* a, const Vector3* b, Vector3* result) noexcept nogil:
    """
    Substract two Vector3 instances and store the result in a third Vector3.
    result = a-b

    Parameters
    ----------
    a : const Vector3*
        Pointer to the first Vector3.
    b : const Vector3*
        Pointer to the second Vector3.
    result : Vector3*
        Pointer to the Vector3 where the result will be stored.
    """
    result[0].data[0] = a[0].data[0] - b[0].data[0]
    result[0].data[1] = a[0].data[1] - b[0].data[1]
    result[0].data[2] = a[0].data[2] - b[0].data[2]

cdef inline void negate_vector3_inplace(Vector3* v) noexcept nogil:
    """
    Negate a Vector3 instance in place by reversing the sign of each component.

    Parameters
    ----------
    v : Vector3*
        Pointer to the Vector3 instance to be negated.

    Returns
    -------
    None
        The input Vector3 is modified in place.
    """
    v[0].data[0] = -v[0].data[0]
    v[0].data[1] = -v[0].data[1]
    v[0].data[2] = -v[0].data[2]

cdef inline void vector3_times_scalar(const Vector3* v, double scalar, Vector3* result) noexcept nogil:
    """
    Multiply a Vector3 instance by a constant scalar.

    Parameters
    ----------
    v : const Vector3*
        Pointer to the input Vector3 instance.
    scalar : double
        The constant scalar by which to multiply the vector.
    result : Vector3*
        Pointer to the output Vector3 instance where the result will be stored.

    Returns
    -------
    None
        The result is stored in the 'result' Vector3 pointer.
    """
    result[0].data[0] = v[0].data[0] * scalar
    result[0].data[1] = v[0].data[1] * scalar
    result[0].data[2] = v[0].data[2] * scalar



cdef inline bint vector3_equals(const Vector3* v1, const Vector3* v2, double tol=1e-9) noexcept nogil:
    """
    Evaluate the equality of two Vector3 instances within a specified tolerance.

    Parameters
    ----------
    v1 : const Vector3*
        Pointer to the first Vector3 instance.
    v2 : const Vector3*
        Pointer to the second Vector3 instance.
    tol : double, optional
        Tolerance for floating-point comparison. Default is 1e-9.

    Returns
    -------
    bint
        True if the two vectors are equal within the specified tolerance, False otherwise.
    """
    return (fabs(v1[0].data[0] - v2[0].data[0]) < tol and
            fabs(v1[0].data[1] - v2[0].data[1]) < tol and
            fabs(v1[0].data[2] - v2[0].data[2]) < tol)

cdef inline double vector3_dot_product(const Vector3* v1, const Vector3* v2) noexcept nogil:
    """
    Compute the dot product of two 3D vectors (Vector3).

    Parameters
    ----------
    v1 : const Vector3*
        Pointer to the first Vector3 instance.
    v2 : const Vector3*
        Pointer to the second Vector3 instance.

    Returns
    -------
    double
        The dot product of the two vectors.
    """
    return v1[0].data[0] * v2[0].data[0] + v1[0].data[1] * v2[0].data[1] + v1[0].data[2] * v2[0].data[2]


cdef inline void matrix3x3_vector3_dot(const Matrix3x3* m, const Vector3* v, Vector3* result) noexcept nogil:
    """
    Compute the dot product of a 3x3 matrix and a 3D vector.

    Parameters
    ----------
    m : const Matrix3x3*
        Pointer to the 3x3 matrix.
    v : const Vector3*
        Pointer to the 3D vector.
    result : Vector3*
        Pointer to the Vector3 where the result will be stored.
    """
    result[0].data[0] = m[0].data[0][0] * v[0].data[0] + m[0].data[0][1] * v[0].data[1] + m[0].data[0][2] * v[0].data[2]
    result[0].data[1] = m[0].data[1][0] * v[0].data[0] + m[0].data[1][1] * v[0].data[1] + m[0].data[1][2] * v[0].data[2]
    result[0].data[2] = m[0].data[2][0] * v[0].data[0] + m[0].data[2][1] * v[0].data[1] + m[0].data[2][2] * v[0].data[2]

cdef inline Vector3 vector3_from_python_object(obj):
    """
    Convert a Python object (list, tuple, or numpy-like array) to a Vector3.

    Parameters
    ----------
    obj : object
        A Python object that is either a list, tuple, or NumPy-like array with 3 elements.

    Returns
    -------
    Vector3
        A new Vector3 struct populated with the data from the Python object.

    Raises
    ------
    ValueError
        If the input object is not a valid type or does not have exactly 3 elements.
    """
    cdef Vector3 v
    cdef int i

    if len(obj) != 3:
        raise ValueError("Input must have exactly 3 elements.")
    try:
        v.data[0] = float(obj[0])
        v.data[1] = float(obj[1])
        v.data[2] = float(obj[2])
    except (TypeError, ValueError):
        raise ValueError("All elements must be convertible to float.")

    return v

cdef inline tuple vector3_to_tuple(Vector3 *v):
    """
    Convert a Vector3 struct to a Python tuple.

    Parameters
    ----------
    v : const Vector3*
        Pointer to the 3D vector.

    Returns
    -------
    tuple
        A tuple containing the three components of the Vector3: (x, y, z).
    """
    return (v[0].data[0], v[0].data[1], v[0].data[2])

cdef inline double vector3_magnitude(const Vector3* v) noexcept nogil:
    """
    Compute the magnitude (norm) of a 3D vector (Vector3).

    Parameters
    ----------
    v : const Vector3*
        Pointer to the Vector3 instance whose magnitude is to be calculated.

    Returns
    -------
    double
        The magnitude of the vector, calculated as the square root of the sum
        of the squares of its components.

    Notes
    -----
    - The magnitude is defined as:
      ||v|| = sqrt(v_x^2 + v_y^2 + v_z^2)
      where v_x, v_y, and v_z are the x, y, and z components of the vector, respectively.
    - This function does not modify the input vector and is safe to use in a `nogil` context.
    """
    # Compute the magnitude of the vector
    return sqrt(v[0].data[0] * v[0].data[0] +
                v[0].data[1] * v[0].data[1] +
                v[0].data[2] * v[0].data[2])



cdef inline void normalize_vector3(Vector3* v) noexcept nogil:
    """
    Normalize a Vector3 to have a unit length.

    Parameters
    ----------
    v : Vector3*
        Pointer to the Vector3 to be normalized. The components of `v` will be modified in place.

    Raises
    ------
    ValueError
        If the magnitude of the vector is zero, indicating that the vector cannot be normalized.
    """
    cdef double magnitude

    # Compute the magnitude of the vector
    magnitude = vector3_magnitude(v)
    
    if magnitude == 0.0:
        raise ValueError("Cannot normalize a vector with zero magnitude.")

    # Normalize the vector components
    v[0].data[0] /= magnitude
    v[0].data[1] /= magnitude
    v[0].data[2] /= magnitude

cdef inline void compute_rotation_matrix(const Vector3* rotation, Matrix3x3* result) noexcept nogil:
    """
    Compute the 3x3 rotation matrix from the given rotation angles (Rx, Ry, Rz) in radians
    and store it in a Matrix3x3 struct.

    The rotation matrix is calculated using the ZYX Euler angles convention:
    1. **First**, a rotation around the Z-axis by angle Rz (roll).
    2. **Second**, a rotation around the Y-axis by angle Ry (yaw).
    3. **Third**, a rotation around the X-axis by angle Rx (pitch).

    The resulting rotation matrix represents the combined effect of these rotations
    applied in the specified order (Z [0]. Y [0]. X).

    Parameters
    ----------
    rotation : const Vector3*
        Pointer to the Vector3 containing the rotation angles (Rx, Ry, Rz) in radians.
        - Rx: Rotation around the X-axis (pitch).
        - Ry: Rotation around the Y-axis (yaw).
        - Rz: Rotation around the Z-axis (roll).
    result : Matrix3x3*
        Pointer to the Matrix3x3 struct where the computed rotation matrix will be stored.

    Notes
    -----
    The ZYX Euler angles convention is commonly used in aerospace, robotics, and computer 
    graphics to describe the orientation of an object or coordinate system relative to a 
    fixed reference frame. In this convention:
    - The first rotation (Rz) is around the fixed Z-axis.
    - The second rotation (Ry) is around the fixed Y-axis.
    - The third rotation (Rx) is around the fixed X-axis.
    Due to the non-commutative nature of 3D rotations, the order in which these rotations
    are applied is crucial.

    Raises
    ------
    None
    """
    # Pre-calculate trigonometric functions
    cdef double cos_rx = cos(rotation[0].data[0])
    cdef double cos_ry = cos(rotation[0].data[1])
    cdef double cos_rz = cos(rotation[0].data[2])

    cdef double sin_rx = sin(rotation[0].data[0])
    cdef double sin_ry = sin(rotation[0].data[1])
    cdef double sin_rz = sin(rotation[0].data[2])

    # Compute the rotation matrix and store it in the result
    result[0].data[0][0] = cos_ry * cos_rz
    result[0].data[0][1] = cos_rz * sin_rx * sin_ry - cos_rx * sin_rz
    result[0].data[0][2] = sin_rx * sin_rz + cos_rx * cos_rz * sin_ry

    result[0].data[1][0] = cos_ry * sin_rz
    result[0].data[1][1] = sin_rx * sin_ry * sin_rz + cos_rx * cos_rz
    result[0].data[1][2] = cos_rx * sin_rz * sin_ry - cos_rz * sin_rx

    result[0].data[2][0] = -sin_ry
    result[0].data[2][1] = cos_ry * sin_rx
    result[0].data[2][2] = cos_rx * cos_ry

cdef void compute_rotation_matrix_i(const Vector3* rotation, Matrix3x3* result) noexcept nogil:
    # Pre-calculate trigonometric functions
    cdef double cos_rx = cos(rotation[0].data[0])
    cdef double cos_ry = cos(rotation[0].data[1])
    cdef double cos_rz = cos(rotation[0].data[2])

    cdef double sin_rx = sin(rotation[0].data[0])
    cdef double sin_ry = sin(rotation[0].data[1])
    cdef double sin_rz = sin(rotation[0].data[2])

    result[0].data[0][0] = cos_ry * cos_rz
    result[0].data[0][1] = cos_ry * sin_rz
    result[0].data[0][2] = -sin_ry

    result[0].data[1][0] = cos_rz * sin_rx * sin_ry - cos_rx * sin_rz
    result[0].data[1][1] = sin_rx * sin_ry * sin_rz + cos_rx * cos_rz
    result[0].data[1][2] = cos_ry * sin_rx

    result[0].data[2][0] = sin_rx * sin_rz + cos_rx * cos_rz * sin_ry
    result[0].data[2][1] = cos_rx * sin_ry * sin_rz - cos_rz * sin_rx
    result[0].data[2][2] = cos_rx * cos_ry

