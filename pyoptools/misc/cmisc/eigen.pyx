# distutils: language=c++

from .eigen cimport Vector3d
from libc.math cimport sin, cos, NAN
from numpy import array


cdef tuple convert_vector3d_to_tuple(Vector3d& v):
    """
    Convert an Eigen::Vector3d to a Python tuple.

    Parameters
    ----------
    v : Vector3d
        An Eigen::Vector3d instance to be converted.

    Returns
    -------
    tuple
        A tuple containing the x, y, and z components of the vector.
    """
    # Access vector components and return them as a Python tuple
    return (v(0), v(1), v(2))

cdef list convert_vectorXd_to_list(VectorXd& v):

    cdef int vect_lenght = v.size()
    cdef int i
    ret_val = []
    for i in range(vect_lenght):
        ret_val.append(v(i))
    return ret_val


cdef inline void assign_to_vector3d(object obj, Vector3d& v):
    """
    Assign values from a Python object (such as a list, tuple, or NumPy array)
    to an existing Eigen::Vector3d.

    Parameters
    ----------
    obj : object
        A Python object that supports indexing (e.g., list, tuple, or NumPy
        array) with at least 3 elements.

    v : Vector3d
        An Eigen::Vector3d instance that will be modified in place.

    Raises
    ------
    ValueError
        If the input object does not have at least 3 elements or if its
        elements cannot be converted to doubles.
    """

    # Ensure that the object has at least 3 elements
    if not hasattr(obj, "__getitem__") or len(obj) < 3:
        raise ValueError("Input object must have at least 3 elements.")

    # Convert and assign values to the Eigen::Vector3d
    cdef double* ptr

    # Assign x-coordinate
    ptr = &v(0)
    ptr[0] = <double>obj[0]  # Convert to double if necessary

    # Assign y-coordinate
    ptr = &v(1)
    ptr[0] = <double>obj[1]  # Convert to double if necessary

    # Assign z-coordinate
    ptr = &v(2)
    ptr[0] = <double>obj[2]  # Convert to double if necessary


cdef inline void assign_tuple_to_vector3d(tuple[double, double, double] t,
                                          Vector3d& v) noexcept nogil:

    """
    Assign values from a tuple to an existing Eigen::Vector3d.

    Parameters
    ----------
    t : tuple[double, double, double]
        A tuple of doubles with 3 elements.

    v : Vector3d
        An Eigen::Vector3d instance that will be modified in place.

    Raises
    ------
    ValueError
        If the input object does not have at least 3 elements or if its
        elements cannot be converted to doubles.
    """

    # Convert and assign values to the Eigen::Vector3d
    cdef double* ptr

    # Assign x-coordinate
    ptr = &v(0)
    ptr[0] = t[0]

    # Assign y-coordinate
    ptr = &v(1)
    ptr[0] = t[1]

    # Assign z-coordinate
    ptr = &v(2)
    ptr[0] = t[2]


cdef inline void assign_tuple_to_vector2d(tuple[double, double] t,
                                          Vector2d& v) noexcept nogil:

    """
    Assign values from a tuple to an existing Eigen::Vector2d.

    Parameters
    ----------
    t : tuple[double, double]
        A tuple of doubles with 2 elements.

    v : Vector2d
        An Eigen::Vector2d instance that will be modified in place.

    Raises
    ------
    ValueError
        If the input object does not have at least 3 elements or if its
        elements cannot be converted to doubles.
    """

    # Convert and assign values to the Eigen::Vector3d
    cdef double* ptr

    # Assign x-coordinate
    ptr = &v(0)
    ptr[0] = t[0]

    # Assign y-coordinate
    ptr = &v(1)
    ptr[0] = t[1]

cdef inline void compute_rotation_matrix(const Vector3d& rotation,
                                         Matrix3d& result) noexcept nogil:
    """
    Compute the 3x3 rotation matrix from the given rotation angles (Rx, Ry,
    Rz) in radians and store it in an Eigen::Matrix3d instance.

    The rotation matrix is calculated using the XYZ Euler angles convention:
    1. **First**, a rotation around the X-axis by angle Rx (pitch).
    2. **Second**, a rotation around the Y-axis by angle Ry (yaw).
    3. **Third**, a rotation around the Z-axis by angle Rz (roll).

    The resulting rotation matrix represents the combined effect of these
    rotations applied in the specified order Rz*Ry*Rx.

    Parameters
    ----------
    rotation : const Eigen::Vector3d
        An Eigen::Vector3d instance containing the rotation angles (Rx, Ry,
        Rz) in radians.
        - Rx: Rotation around the X-axis (pitch).
        - Ry: Rotation around the Y-axis (yaw).
        - Rz: Rotation around the Z-axis (roll).
    result : Eigen::Matrix3d
        An Eigen::Matrix3d instance where the computed rotation matrix
        will be stored.

    Notes
    -----
    The original axes are rotated counterclockwise by the given angle at each
    step.
    """

    # Pre-calculate trigonometric functions
    cdef double cos_rx = cos(rotation(0))
    cdef double cos_ry = cos(rotation(1))
    cdef double cos_rz = cos(rotation(2))

    cdef double sin_rx = sin(rotation(0))
    cdef double sin_ry = sin(rotation(1))
    cdef double sin_rz = sin(rotation(2))

    # Compute the rotation matrix and store it in the result
    (<double *>(&(result(0, 0))))[0] = cos_ry * cos_rz
    (<double *>(&(result(0, 1))))[0] = cos_rz * sin_rx * sin_ry - cos_rx * sin_rz
    (<double *>(&(result(0, 2))))[0] = sin_rx * sin_rz + cos_rx * cos_rz * sin_ry

    (<double *>(&(result(1, 0))))[0] = cos_ry * sin_rz
    (<double *>(&(result(1, 1))))[0] = sin_rx * sin_ry * sin_rz + cos_rx * cos_rz
    (<double *>(&(result(1, 2))))[0] = cos_rx * sin_rz * sin_ry - cos_rz * sin_rx

    (<double *>(&(result(2, 0))))[0] = -sin_ry
    (<double *>(&(result(2, 1))))[0] = cos_ry * sin_rx
    (<double *>(&(result(2, 2))))[0] = cos_rx * cos_ry

cdef inline void compute_rotation_matrix_i(const Vector3d& rotation,
                                           Matrix3d& result) noexcept nogil:
    """
    Compute the 3x3 rotation matrix from the given rotation angles (Rx, Ry,
    Rz) in radians and store it in an Eigen::Matrix3d instance.

    The rotation matrix is calculated using the ZYX Euler angles convention:
    1. **First**, a rotation around the Z-axis by angle Rz (roll).
    2. **Second**, a rotation around the Y-axis by angle Ry (yaw).
    3. **Third**, a rotation around the X-axis by angle Rx (pitch).

    The resulting rotation matrix represents the combined effect of these rotations
    applied in the specified order Rx*Ry*Rz.

    Parameters
    ----------
    rotation : const Eigen::Vector3
        Reference to the Vector3 containing the rotation angles (Rx, Ry, Rz) in radians.
        - Rx: Rotation around the X-axis (pitch).
        - Ry: Rotation around the Y-axis (yaw).
        - Rz: Rotation around the Z-axis (roll).
    result : Eigen::Matrix3d
        An Eigen::Matrix3d instance where the computed rotation matrix
        will be stored.

    Note
    ----
    The original axes are rotated counterclockwise by the given angle at each step.
    """

    # Pre-calculate trigonometric functions
    cdef double cos_rx = cos(rotation(0))
    cdef double cos_ry = cos(rotation(1))
    cdef double cos_rz = cos(rotation(2))

    cdef double sin_rx = sin(rotation(0))
    cdef double sin_ry = sin(rotation(1))
    cdef double sin_rz = sin(rotation(2))

    (<double *>(&(result(0, 0))))[0] = cos_ry * cos_rz
    (<double *>(&(result(0, 1))))[0] = cos_ry * sin_rz
    (<double *>(&(result(0, 2))))[0] = -sin_ry

    (<double *>(&(result(1, 0))))[0] = cos_rz * sin_rx * sin_ry - cos_rx * sin_rz
    (<double *>(&(result(1, 1))))[0] = sin_rx * sin_ry * sin_rz + cos_rx * cos_rz
    (<double *>(&(result(1, 2))))[0] = cos_ry * sin_rx

    (<double *>(&(result(2, 0))))[0] = sin_rx * sin_rz + cos_rx * cos_rz * sin_ry
    (<double *>(&(result(2, 1))))[0] = cos_rx * sin_ry * sin_rz - cos_rz * sin_rx
    (<double *>(&(result(2, 2))))[0] = cos_rx * cos_ry


cdef inline void assign_nan_to_vector3d(Vector3d& v) noexcept nogil:
    """
    Assigns NaN to all components of an Eigen::Vector3d.

    Parameters
    -----------
    v : Eigen::Vector3d
        An Eigen::Vector3d instance to which NaN values will be assigned.
    """

    # Assign NaN to each component of the vector
    (<double*>(&(v(0))))[0] = NAN
    (<double*>(&(v(1))))[0] = NAN
    (<double*>(&(v(2))))[0] = NAN

cdef inline void assign_doubles_to_vector3d(double x, double y, double z,
                                            Vector3d& v) noexcept nogil:
    """
    Assign values to a Vector3d object from individual double components.

    This function assigns the given `x`, `y`, and `z` double values to the
    corresponding components of a `Vector3d` object. The assignment is
    done directly and efficiently, without GIL (Global Interpreter Lock)
    interference, making it suitable for use in performance-critical
    Cython code.

    Parameters
    ----------
    x : double
        The value to assign to the x-component of the vector.
    y : double
        The value to assign to the y-component of the vector.
    z : double
        The value to assign to the z-component of the vector.
    v : Vector3d&
        A reference to the `Vector3d` object to which the values will be assigned.

    Notes
    -----
    - The function uses direct memory access to assign the values, ensuring
      that the operation is performed efficiently.
    - This function is declared as `inline` to minimize function call overhead
      and is marked `noexcept` to indicate that it does not throw exceptions.
    - The `nogil` keyword allows this function to be used in `nogil` sections
      of Cython code, which is important for performance in multi-threaded
      applications.
    """
    (<double*>(&(v(0))))[0] = x
    (<double*>(&(v(1))))[0] = y
    (<double*>(&(v(2))))[0] = z

cdef bint is_approx(Vector3d& a, Vector3d& b, double tol) noexcept nogil:
    """
    Determine if two 3D vectors are approximately equal within a given tolerance.

    This function computes the Euclidean distance between vectors `a` and `b`
    and checks if it is smaller than the specified tolerance `tol`.

    Parameters
    ----------
    a : Eigen::Vector3d
        The first 3D vector.
    b : Eigen::Vector3d
        The second 3D vector.
    tol : double
        The tolerance within which the vectors are considered approximately
        equal.

    Returns
    -------
    bool
        True if the Euclidean distance between the two vectors is less than
        `tol`, False otherwise.
    """
    return (a-b).norm()<=tol


cdef object convert_vector3d_to_array(Vector3d& v):
    """
    Convert an Eigen::Vector3d to a numpy array.

    Parameters
    ----------
    v : Vector3d
        An Eigen::Vector3d instance to be converted.

    Returns
    -------
    array
        An array containing the x, y, and z components of the vector.
    """
    # Access vector components and return them as a Python tuple
    return array((v(0), v(1), v(2)))
