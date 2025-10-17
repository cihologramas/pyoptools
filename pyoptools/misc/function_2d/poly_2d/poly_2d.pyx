
from libc.math cimport sqrt, floor
from pyoptools.misc.cmisc.eigen cimport VectorXd, \
    convert_vectorXd_to_list
import numpy as np

cimport cython

# Import C stdlib
from libc.stdlib cimport *

cdef class Poly2D:
    '''Class to define a 2D polynomial

        .. math::
            z=c0+
            c1*x+c2*y+
            c3*x^2+c4*x*y+c5*y^2+
            c6*x^3+c7*x^2*y+c8*x*y^2+c9*y^3+...
    '''

    def __init__(self, coeff_py):
        """
        """
        cdef int i
        cdef tuple[int, int] powers

        # get the lenght of the coheficient vector
        # Maibe is better to remove

        self._num_coeff = len(coeff_py)

        self._coeff = VectorXd(self._num_coeff)

        # copy the python coefficient list to eigen vector
        for i in range(self._num_coeff):
            (<double*>(&(self._coeff(i))))[0] = coeff_py[i]

        # Save the powers so thy can be usedlatter

        self.px.resize(self._num_coeff)
        self.py.resize(self._num_coeff)

        # Fill the power vectors. This is to speed up the calcularion ahead

        self.order = 0
        for i in range(self._num_coeff):
            powers = index_to_powers(i)
            (<int*>(&(self.px(i))))[0] = powers[0]
            (<int*>(&(self.py(i))))[0] = powers[1]

            # Get the order of the polynomial as the highest power
            if self.order < powers[0]:
                self.order = powers[0]

            if self.order < powers[1]:
                self.order = powers[1]

        self.dx=None
        self.dy=None

    def __add__(self, other):
        """
        Overload the addition operator for Poly2D objects.

        Parameters
        ----------
        other : Poly2D or scalar
            The polynomial or scalar to add.

        Returns
        -------
        Poly2D
            A new Poly2D object representing the sum of the two polynomials or
            NotImplemented if addition is not possible.
        """
        cdef int o1, o2, o, i
        cdef int new_length
        cdef VectorXd ncoeff

        cdef Poly2D __other = <Poly2D> other

        # If both self and other are instances of Poly2D
        if isinstance(self, Poly2D) and isinstance(other, Poly2D):
            # Determine the maximum order of the two polynomials
            o1 = self.order
            o2 = __other.order
            o = max(o1, o2)

            # Calculate the new length for the coefficient vector
            new_length = pxpy2i(0, o) + 1

            # Create a zero-initialized coefficient vector with the new length
            ncoeff.resize(new_length)
            ncoeff.setZero()
            # Copy coefficients from both polynomials to the new vector
            for i in range(self._coeff.size()):
                (<double*>(&(ncoeff(i))))[0]=self._coeff(i)

            for i in range(__other._coeff.size()):
                (<double*>(&(ncoeff(i))))[0]=ncoeff(i)+__other._coeff(i)

            return Poly2D(convert_vectorXd_to_list(ncoeff))

        # If other is a scalar, add it to the constant term
        elif isinstance(other, (int, float)):
            # Create a copy of the coefficients to modify
            ncohef = self._coeff
            (<double*>(&(ncoeff(0))))[0] = ncoeff(0) + other
            return Poly2D(convert_vectorXd_to_list(ncohef))

        # If addition is not possible, return NotImplemented
        return NotImplemented

    def __sub__(self, other):
        """
        Overload the subtraction operator for Poly2D objects.

        Parameters
        ----------
        other : Poly2D or scalar
            The polynomial or scalar to subtract.

        Returns
        -------
        Poly2D
            A new Poly2D object representing the result of the subtraction or
            NotImplemented if subtraction is not possible.
        """
        cdef int o1, o2, o, i
        cdef int new_length
        cdef VectorXd ncoeff

        cdef Poly2D __other = <Poly2D> other

        # If both self and other are instances of Poly2D
        if isinstance(self, Poly2D) and isinstance(other, Poly2D):
            # Determine the maximum order of the two polynomials
            o1 = self.order
            o2 = other.order
            o = max(o1, o2)

            # Calculate the new length for the coefficient vector
            new_length = pxpy2i(0, o) + 1

            # Create a zero-initialized coefficient vector with the new length
            ncoeff.resize(new_length)
            ncoeff.setZero()
            # Copy coefficients from both polynomials to the new vector
            for i in range(self._coeff.size()):
                (<double*>(&(ncoeff(i))))[0]=self._coeff(i)

            for i in range(__other._coeff.size()):
                (<double*>(&(ncoeff(i))))[0]=ncoeff(i)- __other._coeff(i)

            return Poly2D(convert_vectorXd_to_list(ncoeff))

        # If other is a scalar, add it to the constant term
        elif isinstance(other, (int, float)):
            # Create a copy of the coefficients to modify
            ncoeff = self._coeff
            (<double*>(&(ncoeff(0))))[0] = ncoeff(0) - other
            return Poly2D(convert_vectorXd_to_list(ncoeff))

        # If addition is not possible, return NotImplemented
        return NotImplemented

    def __neg__(self):
        """
        Overload the unary negation operator for Poly2D objects.

        Returns
        -------
        Poly2D
            A new Poly2D object representing the negative of the polynomial.
        """
        cdef int i
        cdef VectorXd neg_coeff

        # Create a zero-initialized vector for the negative coefficients
        neg_coeff.resize(self._num_coeff)

        # Negate each coefficient
        for i in range(self._num_coeff):
            (<double*>(&(neg_coeff(i))))[0] = -self._coeff(i)
        return Poly2D(convert_vectorXd_to_list(neg_coeff))

    def __mul__(self, other):
        """
        Overload the multiplication operator for Poly2D objects.

        Parameters
        ----------
        other : Poly2D or scalar
            The polynomial or scalar to multiply.

        Returns
        -------
        Poly2D
            A new Poly2D object representing the result of the multiplication or
            NotImplemented if multiplication is not possible.
        """
        cdef int rxp, ryp, axp, ayp, pxp, pyp
        cdef int o1, o2
        cdef unsigned int i, j, ir
        cdef VectorXd ncoeff

        # Multiplication with a scalar
        if isinstance(other, (float, int)):
            ncoeff.resize(self._num_coeff)
            for i in range(self._num_coeff):
                (<double*>(&(ncoeff(i))))[0] = self._coeff(i) * other
            return Poly2D(convert_vectorXd_to_list(ncoeff))

        # Multiplication with another Poly2D object
        elif isinstance(other, Poly2D):
            o1 = self.order
            o2 = other.order

            # Create a zero-initialized coefficient vector with the new length
            ncoeff.resize(pxpy2i(0, o1 + o2) + 1)
            ncoeff.setZero()

            # Multiply the coefficients of both polynomials
            for i in range(self._num_coeff):
                for j in range(other._num_coeff):
                    axp = other.px(j)
                    ayp = other.py(j)
                    pxp = self.px(i)
                    pyp = self.py(i)

                    rxp = axp + pxp
                    ryp = ayp + pyp

                    ir = pxpy2i(rxp, ryp)
                    # ncohef[ir] += self._coeff[i] * other._coeff[j]
                    (<double*>(&(ncoeff(ir))))[0] = \
                        (<double*>(&(ncoeff(ir))))[0] + \
                        self._coeff(i) * other._coeff(j)
            return Poly2D(convert_vectorXd_to_list(ncoeff))

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

        for i in range(self._num_coeff):
            # Get the powers of x and y for the current term
            px, py = index_to_powers(i)
            c = self._coeff(i)

            # Skip terms with a coefficient of 0
            if c == 0:
                continue

            # Construct the term string
            term = ""

            # Handle coefficient part

            # Add minus for -1, but skip for the constant term
            if c == -1 and (px != 0 or py != 0):
                term += "-"
            # Add the coefficient except when it's 1 (unless it's the constant term)
            elif c != 1 or (px == 0 and py == 0):
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
        retval = " + ".join(terms).replace("+ -", "- ")

        return retval if retval else "0"

    def __repr__(self):
        """
        Generate a detailed string representation of the Poly2D object.

        Returns
        -------
        str
            A string representing the Poly2D object in detail.
        """
        return f"Poly2D(coeff={convert_vectorXd_to_list(self._coeff)}," \
            f" order={self.order})"

    cpdef tuple[Function2D, Function2D] dxdy(self):
        """
        Calculate the derivative with respect to X and Y for the polynomial.

        Returns
        -------
        tuple
            A tuple (dx, dy) containing two Poly2D objects representing the
            derivatives with respect to X and Y, respectively.

        Note
        ----
        This method caches the derivative results the first time it is called.
        If you modify the coeff attribute, the cache will not be updated.
        """
        cdef VectorXd Dx, Dy
        cdef int i, px, py, dxi, dyi

        # Check if derivatives are already cached
        if (self.dx is None) or (self.dy is None):
            # Create zero-initialized vectors for the derivatives
            Dx.resize(self._num_coeff)
            Dx.setZero()
            Dy.resize(self._num_coeff)
            Dy.setZero()

            # Calculate the derivatives
            for i in range(1, self._num_coeff):
                px = self.px(i)
                py = self.py(i)

                if px > 0:  # Only compute Dx if px > 0
                    dxi = pxpy2i(px - 1, py)
                    (<double*>(&(Dx(dxi))))[0] = \
                        (<double*>(&(Dx(dxi))))[0] + px * self._coeff(i)

                if py > 0:  # Only compute Dy if py > 0
                    dyi = pxpy2i(px, py - 1)
                    (<double*>(&(Dy(dyi))))[0] = \
                        (<double*>(&(Dy(dyi))))[0] + py * self._coeff(i)

            # Cache the results
            self.dx = Poly2D(convert_vectorXd_to_list(Dx))
            self.dy = Poly2D(convert_vectorXd_to_list(Dy))

        return self.dx, self.dy

    cdef double eval_cy(self, double x, double y) noexcept nogil:

        cdef int k
        cdef double result = 0
        # TODO: Cache the powers by creating a vector with the powers of x,
        # and power of y

        for k in range(self._num_coeff):
            result += self._coeff(k) * (x ** self.px(k)) * (y ** self.py(k))
        return result

    def eval_1d(self, double[:]x, double[:]y):
        cdef int n_rows = x.shape[0]

        # Allocate a NumPy array for the result
        z_array = np.zeros((n_rows, ), dtype=np.float64)

        cdef double[:] z = z_array
        cdef int i
        for i in range(n_rows):
            z[i] = self.eval_cy(x[i], y[i])

        return z_array

    @property
    def num_coefficients(self):
        return self._num_coeff

    @property
    def coefficients(self):
        return np.array(convert_vectorXd_to_list(self._coeff))

    # @cython.boundscheck(False)  # Disable bounds checking
    # @cython.wraparound(False)   # Disable negative index wraparound
    # @cython.cdivision(True)     # Enable C-like division
    # @cython.nonecheck(False)    # Disable checks for None
    # def mevalr(self, double[:, :] x, double[:, :] y, double rot=0):
    #     """
    #     Evaluate the polynomial for the values given in the 2D matrices x, y, rotated.

    #     Parameters
    #     ----------
    #     x : double[:, :]
    #         A 2D matrix containing the x values where the polynomial is to be
    #         evaluated.
    #     y : double[:, :]
    #         A 2D matrix containing the y values where the polynomial is to be
    #         evaluated.
    #     rot : double
    #         The angle to rotate the coordinate points x, y before evaluating
    #         the polynomial. Given in radians.

    #     Returns
    #     -------
    #     double[:, :]
    #         A 2D matrix with the same shape as x and y, containing the evaluated
    #         polynomial. The coordinates x' and y' are calculated from rotating
    #         x and y.
    #     """
    #     cdef Py_ssize_t i, j, nx, ny, k
    #     cdef double[:, :] Result

    #     # Calculate the rotation matrix
    #     cdef double cosr = cos(rot)
    #     cdef double sinr = sin(rot)

    #     cdef double rx, ry

    #     # Get dimensions of input matrices
    #     nx = x.shape[0]
    #     ny = x.shape[1]

    #     # Initialize result matrix
    #     Result = array(shape=(nx, ny), itemsize=cython.sizeof(double),
    #                    format="d", mode="c")

    #     # Evaluate the polynomial for each element in the matrices x and y
    #     for i in range(nx):
    #         for j in range(ny):
    #             # Calculate rotated coordinates
    #             rx = x[i, j] * cosr - y[i, j] * sinr
    #             ry = x[i, j] * sinr + y[i, j] * cosr

    #             # Compute the polynomial value at each point
    #             for k in range(self._num_coeff):
    #                 Result[i, j] += self._coeff[k] * (rx ** self.px[k]) *
    #                 (ry ** self.py[k])

    #     return Result

    # @cython.boundscheck(False)  # Disable bounds checking
    # @cython.wraparound(False)   # Disable negative index wraparound
    # @cython.cdivision(True)     # Enable C-like division
    # @cython.nonecheck(False)    # Disable checks for None
    # def vveval(self, double[::1] x, double[::1] y):
    #     """
    #     Evaluate the polynomial in a 2D mesh defined by the vectors x and y.

    #     Parameters
    #     ----------
    #     x : double[::1]
    #         A 1D vector with the values of x where the polynomial is to be evaluated.
    #     y : double[::1]
    #         A 1D vector with the values of y where the polynomial is to be evaluated.

    #     Returns
    #     -------
    #     double[:, :]
    #         A 2D matrix with the shape (nx, ny), where nx is the length of
    #         the x vector, and ny is the length of the y vector, containing
    #         the evaluated polynomial.
    #     """
    #     cdef Py_ssize_t i, j, nx, ny, k
    #     cdef double[:, :] Result

    #     # Get dimensions of input vectors
    #     nx = x.shape[0]
    #     ny = y.shape[0]

    #     # Initialize result matrix
    #     Result = array(shape=(nx, ny), itemsize=cython.sizeof(double),
    #                    format="d", mode="c")

    #     # Evaluate the polynomial for each element in the mesh grid
    #     for j in range(nx):
    #         for i in range(ny):
    #             # Initialize polynomial result at (x, y)
    #             for k in range(self._num_coeff):
    #                 Result[j, i] += self._coeff[k] * (x[j] ** self.px[k]) *
    #                  (y[i] ** self.py[k])
    #     return Result

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

    return np.array(x_powers), np.array(y_powers)


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
