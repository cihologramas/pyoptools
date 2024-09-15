from cython cimport inline

# Expose the Vector3 and Matrix3x3 structs
cdef struct Vector3:
    double data[3]

cdef struct Matrix3x3:
    double data[3][3]

# Expose the functions to other Cython modules
cdef Vector3 create_zero_vector3() noexcept nogil
cdef void set_zero_vector3(Vector3* v) noexcept nogil
cdef void set_nan_vector3(Vector3* v) noexcept nogil
cdef void add_vector3(const Vector3* a, const Vector3* b, Vector3* result) noexcept nogil
cdef void substract_vector3(const Vector3* a, const Vector3* b, Vector3* result) noexcept nogil
cdef void negate_vector3_inplace(Vector3* v) noexcept nogil
cdef void vector3_times_scalar(const Vector3* v, double scalar, Vector3* result) noexcept nogil
cdef bint vector3_equals(const Vector3* v1, const Vector3* v2, double tol=*) noexcept nogil
cdef double vector3_dot_product(const Vector3* v1, const Vector3* v2) noexcept nogil
cdef void matrix3x3_vector3_dot(const Matrix3x3* m, const Vector3* v, Vector3* result) noexcept nogil
cdef Vector3 vector3_from_python_object(object obj)
cdef tuple vector3_to_tuple(Vector3 *v)
cdef double vector3_magnitude(const Vector3* v) noexcept nogil
cdef void normalize_vector3(Vector3* v) noexcept nogil
cdef void compute_rotation_matrix(const Vector3* rotation, Matrix3x3* result) noexcept nogil
cdef void compute_rotation_matrix_i(const Vector3* rotation, Matrix3x3* result) noexcept nogil
