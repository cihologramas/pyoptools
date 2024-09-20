# distutils: language=c++

cdef extern from "Eigen/Dense" namespace "Eigen":
    cdef cppclass Vector3d nogil:
        Vector3d()  # Default constructor
        Vector3d(double, double, double)  # Constructor with three double arguments
        double x() const  # Getter for x-coordinate
        double y() const  # Getter for y-coordinate
        double z() const  # Getter for z-coordinate
        double& operator()(int)  # Access elements with index
        # const double& operator()(int) const  # Access elements with index (const version)
        # Vector3d& operator+=(const Vector3d&)  # Addition-assignment operator
        Vector3d operator+(const Vector3d&) const  # Addition operator
        Vector3d operator-(const Vector3d&) const  # Subtraction operator
        Vector3d operator*(double) const  # Scalar multiplication
        Vector3d operator-() const # Unary negation operator
        void normalize()  # Normalize the vector in place, returns reference to itself <- not sure about it returning a reference
        Vector3d normalized() const  # Returns a normalized version of the vector without modifying the original
        double norm() const  # Compute the Euclidean norm (length) of the vector
        double dot(const Vector3d&) const

    # Declaration for Matrix3d
    cdef cppclass Matrix3d nogil:
        Matrix3d()  # Default constructor
        Matrix3d(const Matrix3d&)  # Copy constructor
        double& operator()(int, int)  # Access elements by (row, col)
        #const double& operator()(int, int) const  # Access elements by (row, col) (const version)
        Matrix3d operator+(const Matrix3d&) const  # Addition operator
        Matrix3d operator-(const Matrix3d&) const  # Subtraction operator
        Matrix3d operator*(const Matrix3d&) const  # Multiplication operator (Matrix-Matrix)
        Vector3d operator*(const Vector3d&) const  # Multiplication operator (Matrix-Vector)
        Matrix3d transpose() const  # Transpose function
        Vector3d dot(const Vector3d&) const

    cdef cppclass Vector2d nogil:
        Vector2d()  # Default constructor
        Vector2d(double, double)  # Constructor with three double arguments
        double x() const  # Getter for x-coordinate
        double y() const  # Getter for y-coordinate
        double& operator()(int)  # Access elements with index
        # const double& operator()(int) const  # Access elements with index (const version)
        # Vector3d& operator+=(const Vector2d&)  # Addition-assignment operator
        Vector2d operator+(const Vector2d&) const  # Addition operator
        Vector2d operator-(const Vector2d&) const  # Subtraction operator
        Vector2d operator*(double) const  # Scalar multiplication
        Vector2d operator/(double) const  # Scalar division
        Vector2d operator-() const # Unary negation operator
        void normalize()  # Normalize the vector in place, returns reference to itself <- not sure about it returning a reference
        Vector2d normalized() const  # Returns a normalized version of the vector without modifying the original
        double norm() const  # Compute the Euclidean norm (length) of the vector
        double dot(const Vector2d&) const

cdef tuple convert_vector3d_to_tuple(Vector3d& v)

cdef void assign_to_vector3d(object obj, Vector3d& v)

cdef void assign_tuple_to_vector3d(tuple[double, double, double] t, Vector3d& v) noexcept nogil
cdef void assign_tuple_to_vector2d(tuple[double, double] t, Vector2d& v) noexcept nogil


cdef void compute_rotation_matrix_i(const Vector3d& rotation, Matrix3d& result) noexcept nogil
cdef void compute_rotation_matrix(const Vector3d& rotation, Matrix3d& result) noexcept nogil

cdef void assign_nan_to_vector3d(Vector3d& v) noexcept nogil

cdef void assign_doubles_to_vector3d(double x, double y, double z, Vector3d& v) noexcept nogil

cdef bint is_approx(Vector3d& a, Vector3d& b, double tol) noexcept nogil