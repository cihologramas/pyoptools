cdef extern from "math.h":
    double pow(double, double)
    double sin(double)
    double cos(double)
    double sqrt(double)
    double floor(double)

cimport cython


# Import C stdlib
from libc.stdlib cimport *
from pyoptools.misc.cmisc.cmisc cimport zero_vec, zero_memarray, zero_memarray_2d, zero_memarray_3d, zero_memarray_4d, zero_memarray_5d


cimport cython

from cython.view cimport array

cdef class poly2d:
    '''Class to define a 2D polynomial

        .. math::
            z=c0+
            c1*x+c2*y+
            c3*x^2+c4*x*y+c5*y^2+
            c6*x^3+c7*x^2*y+c8*x*y^2+c9*y^3+...
    '''

    def __cinit__(self, double[::1] coeff):
        """
        """
        cdef int i
        cdef tuple[int, int] powers

        self.coeff = coeff

        # get the lenght of the coheficient vector
        self.clen = coeff.shape[0]

        # Save the powers so thy can be used in a python way or a C fast way
        # Remember if they are not saved as a class attribute, they will get garbage
        # collected

        self.px = array(shape=(self.clen,), itemsize=sizeof(double),
                                    format="d", mode="c")

        self.py = array(shape=(self.clen,), itemsize=sizeof(double),
                                    format="d", mode="c")

        self.order = 0
        for i in range( self.clen):
            powers = index_to_powers(i)
            self.px[i] = powers[0]
            self.py[i] = powers[1]

            # Get the order of the polynomial as the highest power
            if self.order < powers[0]:
                self.order = powers[0]

            if self.order < powers[1]:
                self.order = powers[1]

        self.dx=None
        self.dy=None

    # def __reduce__(self):
    #    """
    #    Support pickling of the poly2d object.
    #
    #    Returns
    #    -------
    #    tuple
    #        A tuple that contains the class and the arguments needed to
    #        recreate the object.
    #    """
    #    # Convert memoryview to bytes for serialization
    #    coeff_bytes = self.coeff.tobytes()
    #    coeff_shape = self.coeff.shape
    #    return (self.__class__, (memoryview(coeff_bytes).cast('d', shape=coeff_shape),))


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    def __add__(self, other):
        """
        Overload the addition operator for poly2d objects.

        Parameters
        ----------
        other : poly2d or scalar
            The polynomial or scalar to add.

        Returns
        -------
        poly2d
            A new poly2d object representing the sum of the two polynomials or
            NotImplemented if addition is not possible.
        """
        cdef int o1, o2, o, i
        cdef int new_length
        cdef double[::1] ncohef

        # If both self and other are instances of poly2d
        if isinstance(self, poly2d) and isinstance(other, poly2d):
            # Determine the maximum order of the two polynomials
            o1 = self.order
            o2 = other.order
            o = max(o1, o2)

            # Calculate the new length for the coefficient vector
            new_length = pxpy2i(0, o) + 1

            # Create a zero-initialized coefficient vector with the new length
            ncohef = zero_vec(new_length)

            # Copy coefficients from both polynomials to the new vector
            ncohef[:self.clen] = self.coeff

            for i in range(self.clen):
                ncohef[i] += other.coeff[i]

            return poly2d(ncohef)

        # If other is a scalar, add it to the constant term
        elif isinstance(other, (int, float)):
            # Create a copy of the coefficients to modify
            ncohef = zero_vec(self.clen)
            ncohef[:] = self.coeff
            ncohef[0] += other
            return poly2d(ncohef)

        # If addition is not possible, return NotImplemented
        return NotImplemented

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    def __sub__(self, other):
        """
        Overload the subtraction operator for poly2d objects.

        Parameters
        ----------
        other : poly2d or scalar
            The polynomial or scalar to subtract.

        Returns
        -------
        poly2d
            A new poly2d object representing the result of the subtraction or
            NotImplemented if subtraction is not possible.
        """
        cdef int o1, o2, o, i
        cdef int new_length
        cdef double[::1] ncohef

        # If both self and other are instances of poly2d
        if isinstance(self, poly2d) and isinstance(other, poly2d):
            # Determine the maximum order of the two polynomials
            o1 = self.order
            o2 = other.order
            o = max(o1, o2)

            # Calculate the new length for the coefficient vector
            new_length = pxpy2i(0, o) + 1

            # Create a zero-initialized coefficient vector with the new length
            ncohef = zero_vec(new_length)

            # Copy coefficients from both polynomials to the new vector
            ncohef[:self.clen] = self.coeff

            for i in range(self.clen):
                ncohef[i] -= other.coeff[i]

            return poly2d(ncohef)

        # If other is a scalar, subtract it from the constant term
        elif isinstance(other, (int, float)):
            # Create a copy of the coefficients to modify
            ncohef = zero_vec(self.clen)
            ncohef[:] = self.coeff
            ncohef[0] -= other
            return poly2d(ncohef)

        # If subtraction is not possible, return NotImplemented
        return NotImplemented

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    def __neg__(self):
        """
        Overload the unary negation operator for poly2d objects.

        Returns
        -------
        poly2d
            A new poly2d object representing the negative of the polynomial.
        """
        cdef int i
        cdef double[::1] neg_coeff

        # Create a zero-initialized vector for the negative coefficients
        neg_coeff = zero_vec(self.clen)

        # Negate each coefficient
        for i in range(self.clen):
            neg_coeff[i] = -self.coeff[i]

        return poly2d(neg_coeff)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    def __mul__(self, other):
        """
        Overload the multiplication operator for poly2d objects.

        Parameters
        ----------
        other : poly2d or scalar
            The polynomial or scalar to multiply.

        Returns
        -------
        poly2d
            A new poly2d object representing the result of the multiplication or
            NotImplemented if multiplication is not possible.
        """
        cdef int rxp, ryp, axp, ayp, pxp, pyp
        cdef int o1, o2
        cdef unsigned int i, j, ir
        cdef double[::1] ncohef

        # Multiplication with a scalar
        if isinstance(other, (float, int)):
            ncohef = zero_vec(self.clen)
            for i in range(self.clen):
                ncohef[i] = self.coeff[i] * other
            return poly2d(ncohef)

        # Multiplication with another poly2d object
        elif isinstance(other, poly2d):
            o1 = self.order
            o2 = other.order

            # Create a zero-initialized coefficient vector with the new length
            ncohef = zero_vec(pxpy2i(0, o1 + o2) + 1)

            # Multiply the coefficients of both polynomials
            for i in range(self.clen):
                for j in range(other.clen):
                    axp = other.px[j]
                    ayp = other.py[j]
                    pxp = self.px[i]
                    pyp = self.py[i]

                    rxp = axp + pxp
                    ryp = ayp + pyp

                    ir = pxpy2i(rxp, ryp)
                    ncohef[ir] += self.coeff[i] * other.coeff[j]

            return poly2d(ncohef)

        # If multiplication is not possible, return NotImplemented
        return NotImplemented

    def __str__(self):
        """
        Generate a human-readable string representation of the polynomial.

        Returns
        -------
        str
            A string representing the polynomial.
        """
        cdef int i, px, py
        cdef double c
        cdef list terms = []
        cdef str term, retval = ""

        for i in range(self.clen):
            # Get the powers of x and y for the current term
            px, py = index_to_powers(i)
            c = self.coeff[i]

            # Skip terms with a coefficient of 0
            if c == 0:
                continue

            # Construct the term string
            term = ""

            # Handle coefficient part
            if c == -1 and (px != 0 or py != 0):  # Add minus for -1, but skip for the constant term
                term += "-"
            elif c != 1 or (px == 0 and py == 0):  # Add the coefficient except when it's 1 (unless it's the constant term)
                term += str(c)

            # Handle the x part
            if px != 0:
                term += "x"
                if px != 1:
                    term += "^" + str(px)

            # Handle the y part
            if py != 0:
                term += "y"
                if py != 1:
                    term += "^" + str(py)

            # Append the constructed term to the list
            terms.append(term)

        # Join all terms with a plus sign, handle empty cases gracefully
        retval = " + ".join(terms).replace("+-", "- ")

        return retval if retval else "0"

    def __repr__(self):
        """
        Generate a detailed string representation of the poly2d object.

        Returns
        -------
        str
            A string representing the poly2d object in detail.
        """
        return f"poly2d(coeff={list(self.coeff)}, order={self.order})"

    @cython.boundscheck(False)   # Disable bounds checking for array accesses
    @cython.wraparound(False)    # Disable negative index wraparound
    @cython.cdivision(True)      # Enable C-like division, skips division by zero checks
    @cython.nonecheck(False)     # Disable checks for None when accessing object attributes
    @cython.infer_types(True)    # Allow Cython to infer types automatically
    def dxdy(self):
        """
        Calculate the derivative with respect to X and Y for the polynomial.

        Returns
        -------
        tuple
            A tuple (dx, dy) containing two poly2d objects representing the
            derivatives with respect to X and Y, respectively.

        Note
        ----
        This method caches the derivative results the first time it is called.
        If you modify the coeff attribute, the cache will not be updated.
        """
        cdef double[::1] Dx, Dy
        cdef int i, px, py, dxi, dyi

        # Check if derivatives are already cached
        if (self.dx is None) or (self.dy is None):
            # Create zero-initialized vectors for the derivatives
            Dx = zero_vec(self.clen)
            Dy = zero_vec(self.clen)

            # Calculate the derivatives
            for i in range(1, self.clen):
                px = self.px[i]
                py = self.py[i]

                if px > 0:  # Only compute Dx if px > 0
                    dxi = pxpy2i(px - 1, py)
                    Dx[dxi] += px * self.coeff[i]

                if py > 0:  # Only compute Dy if py > 0
                    dyi = pxpy2i(px, py - 1)
                    Dy[dyi] += py * self.coeff[i]

            # Cache the results
            self.dx = poly2d(Dx)
            self.dy = poly2d(Dy)

        return self.dx, self.dy

    
    cpdef eval(self, x, y):
        """
        Evaluate the polynomial at x, y for input arrays with 2 to 5 dimensions.

        Parameters
        ----------
        x : double[:, :] or double[:, :, :] or double[:, :, :, :] or double[:, :, :, :, :]
            A memoryview containing the x values where the polynomial is to be evaluated.
        y : double[:, :] or double[:, :, :] or double[:, :, :, :] or double[:, :, :, :, :]
            A memoryview containing the y values where the polynomial is to be evaluated.

        Returns
        -------
        double[:, :], double[:, :, :], double[:, :, :, :], or double[:, :, :, :, :]
            A memoryview of type double with the same shape as x and y, representing
            the polynomial evaluated at (x, y).
        """
        cdef int ndim = len(x.shape)

        if ndim == 2:
            return self.eval2d(x, y)
        elif ndim == 3:
            return self.eval3d(x, y)
        elif ndim == 4:
            return self.eval4d(x, y)
        elif ndim == 5:
            return self.eval5d(x, y)
        else:
            raise ValueError("Unsupported number of dimensions: %d. Supported dimensions are 2, 3, 4, or 5." % ndim)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    cdef double[:, :] eval2d(self, double[:, :] x, double[:, :] y):
        """
        Evaluate the polynomial for 2D input arrays.

        Parameters
        ----------
        x : double[:, :]
            A 2D memoryview containing the x values.
        y : double[:, :]
            A 2D memoryview containing the y values.

        Returns
        -------
        double[:, :]
            A 2D memoryview of doubles representing the polynomial evaluated at (x, y).
        """
        cdef Py_ssize_t i, j, nx, ny, k
        cdef double[:, :] result

        # Get dimensions
        nx = x.shape[0]
        ny = x.shape[1]

        # Initialize result memoryview with the same shape
        result = zero_memarray_2d(nx, ny)

        # Evaluate the polynomial for each element in the matrices x and y
        for i in range(nx):
            for j in range(ny):
                for k in range(self.clen):
                    result[i, j] += self.coeff[k] * (x[i, j] ** self.px[k]) * (y[i, j] ** self.py[k])

        return result


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    cdef double[:, :, :] eval3d(self, double[:, :, :] x, double[:, :, :] y):
        """
        Evaluate the polynomial for 3D input arrays.

        Parameters
        ----------
        x : double[:, :, :]
            A 3D memoryview containing the x values.
        y : double[:, :, :]
            A 3D memoryview containing the y values.

        Returns
        -------
        double[:, :, :]
            A 3D memoryview of doubles representing the polynomial evaluated at (x, y).
        """
        cdef Py_ssize_t i, j, k, nx, ny, nz, m
        cdef double[:, :, :] result

        # Get dimensions
        nx = x.shape[0]
        ny = x.shape[1]
        nz = x.shape[2]

        # Initialize result memoryview with the same shape
        result = zero_memarray_3d(nx, ny, nz)

        # Evaluate the polynomial for each element in the matrices x and y
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for m in range(self.clen):
                        result[i, j, k] += self.coeff[m] * (x[i, j, k] ** self.px[m]) * (y[i, j, k] ** self.py[m])

        return result


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    cdef double[:, :, :, :] eval4d(self, double[:, :, :, :] x, double[:, :, :, :] y):
        """
        Evaluate the polynomial for 4D input arrays.

        Parameters
        ----------
        x : double[:, :, :, :]
            A 4D memoryview containing the x values.
        y : double[:, :, :, :]
            A 4D memoryview containing the y values.

        Returns
        -------
        double[:, :, :, :]
            A 4D memoryview of doubles representing the polynomial evaluated at (x, y).
        """
        cdef Py_ssize_t i, j, k, l, nx, ny, nz, nw, m
        cdef double[:, :, :, :] result

        # Get dimensions
        nx = x.shape[0]
        ny = x.shape[1]
        nz = x.shape[2]
        nw = x.shape[3]

        # Initialize result memoryview with the same shape
        result = zero_memarray_4d(nx, ny, nz, nw)

        # Evaluate the polynomial for each element in the matrices x and y
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for l in range(nw):
                        for m in range(self.clen):
                            result[i, j, k, l] += self.coeff[m] * (x[i, j, k, l] ** self.px[m]) * (y[i, j, k, l] ** self.py[m])

        return result


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    cdef double[:, :, :, :, :] eval5d(self, double[:, :, :, :, :] x, double[:, :, :, :, :] y):
        """
        Evaluate the polynomial for 5D input arrays.

        Parameters
        ----------
        x : double[:, :, :, :, :]
            A 5D memoryview containing the x values.
        y : double[:, :, :, :, :]
            A 5D memoryview containing the y values.

        Returns
        -------
        double[:, :, :, :, :]
            A 5D memoryview of doubles representing the polynomial evaluated at (x, y).
        """
        cdef Py_ssize_t i, j, k, l, m, nx, ny, nz, nw, nv, p
        cdef double[:, :, :, :, :] result

        # Get dimensions
        nx = x.shape[0]
        ny = x.shape[1]
        nz = x.shape[2]
        nw = x.shape[3]
        nv = x.shape[4]

        # Initialize result memoryview with the same shape
        result = zero_memarray_5d(nx, ny, nz, nw, nv)

        # Evaluate the polynomial for each element in the matrices x and y
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    for l in range(nw):
                        for m in range(nv):
                            for p in range(self.clen):
                                result[i, j, k, l, m] += self.coeff[p] * (x[i, j, k, l, m] ** self.px[p]) * (y[i, j, k, l, m] ** self.py[p])

        return result


    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    @cython.infer_types(True)
    cpdef double peval(self, double x, double y):
        """
        Evaluate the polynomial at a single point (x, y).

        Parameters
        ----------
        x : double
            The x value where the polynomial is to be evaluated.
        y : double
            The y value where the polynomial is to be evaluated.

        Returns
        -------
        double
            The result of the polynomial evaluated at (x, y).
        """
        cdef Py_ssize_t i
        cdef double result = 0.0
        cdef double x_power, y_power

        cdef double[::1] coeff = self.coeff
        cdef int[::1] px = self.px
        cdef int[::1] py = self.py

        # Evaluate the polynomial using a manual loop
        for i in range(self.clen):
            if coeff[i] != 0:
                # Calculate x^px[i] and y^py[i]
                x_power = x ** px[i] if px[i] != 0 else 1.0
                y_power = y ** py[i] if py[i] != 0 else 1.0
                
                # Add the term to the result
                result += coeff[i] * x_power * y_power

        return result




    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    @cython.nonecheck(False)
    def meval(self, double[:, :] x, double[:, :] y):
        """
        Evaluate the polynomial for the values given in the 2D matrices x, y.

        Parameters
        ----------
        x : double[:, :]
            A 2D matrix containing the x values where the polynomial is to be evaluated.
        y : double[:, :]
            A 2D matrix containing the y values where the polynomial is to be evaluated.

        Returns
        -------
        double[:, :]
            A 2D matrix with the same shape as x and y, containing the evaluated polynomial.
        """
        cdef Py_ssize_t i, j, nx, ny, k
        cdef double[:, :] Result

        # Get dimensions of input matrices
        nx = x.shape[0]
        ny = x.shape[1]

        # Initialize result matrix
        Result = array(shape=(nx, ny), itemsize=cython.sizeof(double), format="d",
                       mode="c")

        # Evaluate the polynomial for each element in the matrices x and y
        for i in range(nx):
            for j in range(ny):
                for k in range(self.clen):
                    Result[i, j] += self.coeff[k] * (x[i, j] ** self.px[k]) * (y[i, j] ** self.py[k])

        return Result


    @cython.boundscheck(False)  # Disable bounds checking
    @cython.wraparound(False)   # Disable negative index wraparound
    @cython.cdivision(True)     # Enable C-like division
    @cython.nonecheck(False)    # Disable checks for None
    def mevalr(self, double[:, :] x, double[:, :] y, double rot=0):
        """
        Evaluate the polynomial for the values given in the 2D matrices x, y, rotated.

        Parameters
        ----------
        x : double[:, :]
            A 2D matrix containing the x values where the polynomial is to be evaluated.
        y : double[:, :]
            A 2D matrix containing the y values where the polynomial is to be evaluated.
        rot : double
            The angle to rotate the coordinate points x, y before evaluating
            the polynomial. Given in radians.

        Returns
        -------
        double[:, :]
            A 2D matrix with the same shape as x and y, containing the evaluated polynomial.
            The coordinates x' and y' are calculated from rotating x and y.
        """
        cdef Py_ssize_t i, j, nx, ny, k
        cdef double[:, :] Result

        # Calculate the rotation matrix
        cdef double cosr = cos(rot)
        cdef double sinr = sin(rot)

        cdef double rx, ry

        # Get dimensions of input matrices
        nx = x.shape[0]
        ny = x.shape[1]

        # Initialize result matrix
        Result = array(shape=(nx, ny), itemsize=cython.sizeof(double),
                       format="d", mode="c")

        # Evaluate the polynomial for each element in the matrices x and y
        for i in range(nx):
            for j in range(ny):
                # Calculate rotated coordinates
                rx = x[i, j] * cosr - y[i, j] * sinr
                ry = x[i, j] * sinr + y[i, j] * cosr

                # Compute the polynomial value at each point
                for k in range(self.clen):
                    Result[i, j] += self.coeff[k] * (rx ** self.px[k]) * (ry ** self.py[k])

        return Result


    @cython.boundscheck(False)  # Disable bounds checking
    @cython.wraparound(False)   # Disable negative index wraparound
    @cython.cdivision(True)     # Enable C-like division
    @cython.nonecheck(False)    # Disable checks for None
    def vveval(self, double[::1] x, double[::1] y):
        """
        Evaluate the polynomial in a 2D mesh defined by the vectors x and y.

        Parameters
        ----------
        x : double[::1]
            A 1D vector with the values of x where the polynomial is to be evaluated.
        y : double[::1]
            A 1D vector with the values of y where the polynomial is to be evaluated.

        Returns
        -------
        double[:, :]
            A 2D matrix with the shape (nx, ny), where nx is the length of the x vector,
            and ny is the length of the y vector, containing the evaluated polynomial.
        """
        cdef Py_ssize_t i, j, nx, ny, k
        cdef double[:, :] Result

        # Get dimensions of input vectors
        nx = x.shape[0]
        ny = y.shape[0]

        # Initialize result matrix
        Result = array(shape=(nx, ny), itemsize=cython.sizeof(double),
                       format="d", mode="c")

        # Evaluate the polynomial for each element in the mesh grid
        for j in range(nx):
            for i in range(ny):
                # Initialize polynomial result at (x, y)
                for k in range(self.clen):
                    Result[j, i] += self.coeff[k] * (x[j] ** self.px[k]) * (y[i] ** self.py[k])

        return Result




cpdef tuple[int, int] index_to_powers(int i):
    """
    Convert index of coefficient vector to x and y powers in a polynomial.

    ========= = == == == == == == == == == == == == == == ===
    index     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 ...
    xpower    0  1  0  2  1  0  3  2  1  0  4  3  2  1  0 ...
    ypower    0  0  1  0  1  2  0  1  2  3  0  1  2  3  4 ...
    pol_order 0  1  1  2  2  2  3  3  3  3  4  4  4  4  4 ...
    ========= = == == == == == == == == == == == == == == ===

    Parameters
    ----------
    i : int
        Index in the coefficient vector.

    Returns
    -------
    (int, int)
        A tuple (x_power, y_power) representing the powers of x and y for the term.
    """
    cdef int d, index_within_degree
    cdef double d_double

    # Calculate degree d directly using the formula
    d_double = 1+(-3.0 + sqrt(1.0 + 8.0 * i)) / 2.0
    d = <int> floor(d_double)  # Convert to integer by flooring

    # Calculate the index within the terms of degree d
    index_within_degree = i - (d * (d + 1)) // 2

    # Compute the powers of x and y
    x_power = d - index_within_degree
    y_power = index_within_degree

    return x_power, y_power

from cython cimport view

cpdef tuple indices_to_powers(int[:] indices):
    """
    Convert a vector of coefficient indices into the corresponding x and y 
    exponents for a polynomial.
    """
    cdef int n = indices.shape[0]
    
    # Create 1D arrays for x and y powers
    cdef int[:] x_powers = view.array(shape=(n,), 
                                      itemsize=cython.sizeof(cython.int),
                                      format="i", mode="c")
    cdef int[:] y_powers = view.array(shape=(n,), 
                                      itemsize=cython.sizeof(cython.int),
                                      format="i", mode="c")
    
    cdef int idx

    # Loop to compute x and y powers based on the index
    for idx in range(n):
        # Directly unpack the tuple returned by index_to_powers
        x_powers[idx], y_powers[idx] = index_to_powers(indices[idx])

    return x_powers, y_powers


cpdef int pxpy2i(int px, int py):
    """Method that returns the index given power in x and the power in y
    """
    cdef int po=px+py
    cdef int i0=(((po-1)*2+3)**2-1)//8
    return i0+py

cpdef int ord2i(int o):
    """
    Method that returns the number of coefficients of a polynomial of order o

    ===== == == == == == ===
    order  0  1  2  3  4 ...
    i      1  3  6 10 15 ...
    ===== == == == == == ===
    """
    return (o+2)*(o+1)//2
