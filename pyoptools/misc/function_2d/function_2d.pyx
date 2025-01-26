from pyoptools.misc.cmisc.eigen cimport MatrixXd
from numpy import empty, float64
cimport cython

cdef class Function2D:
    """
    Abstract base class for 2D functions defined in Cartesian coordinates.

    This class provides a framework for representing and evaluating 2D mathematical
    functions, including their partial derivatives. All concrete implementations
    should inherit from this class and implement the required methods.

    Notes
    -----
    The class provides two main evaluation methods:
    - eval2d: A nogil Cython method using Eigen matrices
    - eval: A Python-accessible method using NumPy arrays

    Subclasses must implement:
    - eval_cy: Core evaluation function
    - dx: x-derivative
    - dy: y-derivative
    """

    cdef void eval2d(self, MatrixXd& x, MatrixXd& y, MatrixXd& result) noexcept nogil:
        """
        Evaluate the function for matrices of x and y coordinates using Eigen.

        Parameters
        ----------
        x : MatrixXd&
            Matrix of x coordinates
        y : MatrixXd&
            Matrix of y coordinates
        result : MatrixXd&
            Output matrix to store results

        Notes
        -----
        This is a nogil implementation using Eigen matrices for high performance.
        """
        cdef Py_ssize_t i, j, n_cols, n_rows

        n_cols = x.cols()
        n_rows = x.rows()
        result.resize(n_rows, n_cols)

        for i in range(n_rows):
            for j in range(n_cols):
                (<double*>(&(result(i, j))))[0] = self.eval_cy(x(i, j), y(i, j))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def eval(self, double[:, ::1] x, double[:, ::1] y):
        """
        Evaluate the function for arrays of x and y coordinates.

        Parameters
        ----------
        x : ndarray
            2D array of x coordinates
        y : ndarray
            2D array of y coordinates

        Returns
        -------
        ndarray
            2D array containing function values at each (x,y) point

        Notes
        -----
        This is the Python-accessible implementation using NumPy arrays.
        Decorated for maximum performance with bounds checking disabled.
        """
        cdef int n_rows = x.shape[0]
        cdef int n_cols = x.shape[1]

        z_array = empty((n_rows, n_cols), dtype=float64)
        cdef double[:, ::1] z = z_array

        cdef int i, j
        for i in range(n_rows):
            for j in range(n_cols):
                z[i, j] = self.eval_cy(x[i, j], y[i, j])

        return z_array

    cdef double eval_cy(self, double x, double y) noexcept nogil:
        """
        Core evaluation method for a single point.

        Parameters
        ----------
        x : double
            x coordinate
        y : double
            y coordinate

        Returns
        -------
        double
            Function value at (x,y)

        Notes
        -----
        This method must be implemented by all subclasses.
        It is the fundamental evaluation method used by both eval() and eval2d().
        """
        with gil:
            raise NotImplementedError("Subclasses must implement eval_cy()")

    cpdef tuple[Function2D, Function2D] dxdy(self):
        """
        Compute both partial derivatives simultaneously.

        Returns
        -------
        tuple[Function2D, Function2D]
            A tuple containing (∂f/∂x, ∂f/∂y)

        Notes
        -----
        This method must be implemented by all subclasses.
        """
        raise NotImplementedError("Subclasses must implement dxdy()")
