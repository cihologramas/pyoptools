cdef double[::1] norm_vect(double[::1] v)
cdef double vector_length(double[::1] v)

cdef double[::, ::1] rot_mat_i(double[::1] r)
cdef double[::, ::1] rot_mat(double[::1] r)

cdef double[::1] empty_vec(int length)
cdef double[::1] zero_vec(int length)
cdef double[::1] to_vector(object py_obj)

cpdef zero_memarray(tuple shape)
cdef double[:, ::1] zero_memarray_2d(Py_ssize_t nx, Py_ssize_t ny)
cdef double[:, :, ::1] zero_memarray_3d(Py_ssize_t nx, Py_ssize_t ny, Py_ssize_t nz)


cpdef double[::1] dot_product_3x3_matrix_vector(double[:, :] matrix, double[:] vector)

cdef bint allclose_cython(double[::1] a, double[::1] b, double atol)

cdef double norm_3d_vector(double[::1] vec)
