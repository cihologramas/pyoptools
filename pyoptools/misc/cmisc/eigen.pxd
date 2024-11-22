# distutils: language=c++

from libc.stddef cimport size_t  # Import size_t

cdef extern from "<Eigen/Dense>" namespace "Eigen":
    cdef cppclass Vector3d nogil:
        Vector3d()
        Vector3d(double, double, double)
        double x() const
        double y() const
        double z() const
        double& operator()(int)
        Vector3d operator+(const Vector3d&) const
        Vector3d operator-(const Vector3d&) const
        Vector3d operator*(double) const
        Vector3d operator/(double) const
        Vector3d operator-() const
        void normalize()
        Vector3d normalized() const
        double norm() const
        double dot(const Vector3d&) const

    # Declaration for Matrix3d
    cdef cppclass Matrix3d nogil:
        Matrix3d()
        Matrix3d(const Matrix3d&)
        double& operator()(int, int)
        Matrix3d operator+(const Matrix3d&) const
        Matrix3d operator-(const Matrix3d&) const
        Matrix3d operator*(const Matrix3d&) const
        Vector3d operator*(const Vector3d&) const
        Matrix3d transpose() const
        Vector3d dot(const Vector3d&) const

    cdef cppclass Vector2d nogil:
        Vector2d()
        Vector2d(double, double)
        double x() const
        double y() const
        double& operator()(int)
        Vector2d operator+(const Vector2d&) const
        Vector2d operator-(const Vector2d&) const
        Vector2d operator*(double) const
        Vector2d operator/(double) const
        Vector2d operator-() const
        void normalize()
        Vector2d normalized() const
        double norm() const
        double dot(const Vector2d&) const

    cdef cppclass MatrixXd nogil:
        MatrixXd()
        MatrixXd(int rows, int cols)
        MatrixXd(const MatrixXd&)
        double& operator()(int, int)
        MatrixXd operator+(const MatrixXd&) const
        MatrixXd operator-(const MatrixXd&) const
        MatrixXd operator*(const MatrixXd&) const
        MatrixXd transpose() const
        double determinant() const
        int rows() const
        int cols() const
        void resize(int newRows, int newCols)

    cdef cppclass VectorXd nogil:
        VectorXd()
        VectorXd(int size)
        VectorXd(const VectorXd&)
        double& operator()(int index)
        VectorXd operator+(const VectorXd&) const
        VectorXd operator-(const VectorXd&) const
        double dot(const VectorXd&) const
        VectorXd normalized() const
        double norm() const
        int size() const
        void setZero()
        void resize(int newSize)

    cdef cppclass VectorXi nogil:
        VectorXi()
        VectorXi(int size)
        VectorXi(const VectorXi&)
        int& operator()(int index)
        int size() const
        VectorXi operator+(const VectorXi&) const
        VectorXi operator-(const VectorXi&) const
        int dot(const VectorXi&) const
        VectorXi normalized() const
        double norm() const
        void setZero()
        void resize(int newSize)

cdef tuple convert_vector3d_to_tuple(Vector3d& v)
cdef list convert_vectorXd_to_list(VectorXd& v)


cdef void assign_to_vector3d(object obj, Vector3d& v)

cdef void assign_tuple_to_vector3d(tuple[double, double, double] t,
                                   Vector3d& v) noexcept nogil
cdef void assign_tuple_to_vector2d(tuple[double, double] t,
                                   Vector2d& v) noexcept nogil


cdef void compute_rotation_matrix_i(const Vector3d& rotation,
                                    Matrix3d& result) noexcept nogil
cdef void compute_rotation_matrix(const Vector3d& rotation,
                                  Matrix3d& result) noexcept nogil

cdef void assign_nan_to_vector3d(Vector3d& v) noexcept nogil

cdef void assign_doubles_to_vector3d(double x,
                                     double y,
                                     double z,
                                     Vector3d& v) noexcept nogil

cdef bint is_approx(Vector3d& a, Vector3d& b, double tol) noexcept nogil

cdef object convert_vector3d_to_array(Vector3d& v)
