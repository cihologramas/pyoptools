"""Module defining the Zernike polynomials
"""

# The location of some scipy functions has changed

try:
    from scipy import factorial, comb as binomial
except ImportError:
    # both moved to special in scipy1.0.0
    from scipy.special import comb as binomial
    from scipy.special import factorial

from numpy import (
    mgrid,
    sqrt,
    arccos,
    zeros,
    transpose,
    pi,
    cos,
    sin,
    ones,
    array,
    where,
)

from numpy.ma import masked_array
from pyoptools.misc.function_2d.poly_2d.poly_2d import ord2i, indices_to_powers as i2pxpy, Poly2D


def polar_array(Rmax=1.0, DS=0.1, pr=1.0):
    """
    Function that generates 2 square matrices one with the rho coordinate,
    and the other with the theta coordinate, to be able to calculate
    functions using polar coordinates.

    It is similar to the mgrid function.

    **Arguments**

    ==== ===================================================
    Rmax Limit the pupil area -Rmax<=X<=Rmax -Rmax<=Y<=Rmax
    DS   Step between pixels
    pr   Pupil radius. Used to normalize the pupil.
    ==== ===================================================

    ..  TODO:: This function should be moved to a auxiliary functions module
    """

    X, Y = mgrid[-Rmax:Rmax + DS:DS, -Rmax:Rmax + DS:DS] / pr
    r = sqrt(X**2 + Y**2)
    th = arccos(transpose(X * 1.0 / r))
    th = where(th < 2.0 * pi, th, 0)
    th = where(X < 0, 2.0 * pi - th, th)
    return r, th


def rnm(n, m, rho):
    """
    Return an array with the zernike Rnm polynomial calculated at rho points.


    **ARGUMENTS:**

        === ==========================================
        n    n order of the Zernike polynomial
        m    m order of the Zernike polynomial
        rho  Matrix containing the radial coordinates.
        === ==========================================

    .. note:: For rho>1 the returned value is 0

    .. note:: Values for rho<0 are silently returned as rho=0

    """

    if type(n) is not int:
        raise Exception("n must be integer")
    if type(m) is not int:
        raise Exception("m must be integer")
    if (n - m) % 2 != 0:
        raise Exception("n-m must be even")
    if abs(m) > n:
        raise Exception("The following must be true |m|<=n")
    mask = where(rho <= 1, False, True)

    if n == 0 and m == 0:
        return masked_array(data=ones(rho.shape), mask=mask)
    rho = where(rho < 0, 0, rho)
    Rnm = zeros(rho.shape)
    S = (n - abs(m)) / 2
    for s in range(0, S + 1):
        CR = (
            pow(-1, s)
            * factorial(n - s)
            / (
                factorial(s)
                * factorial(-s + (n + abs(m)) / 2)
                * factorial(-s + (n - abs(m)) / 2)
            )
        )
        p = CR * pow(rho, n - 2 * s)
        Rnm = Rnm + p
    return masked_array(data=Rnm, mask=mask)


def zernike(n, m, rho, theta):
    """
    Returns the an array with the Zernike polynomial evaluated in the rho and
    theta.

    **ARGUMENTS:**

    ===== ==========================================
    n     n order of the Zernike polynomial
    m     m order of the Zernike polynomial
    rho   Matrix containing the radial coordinates.
    theta Matrix containing the angular coordinates.
    ===== ==========================================

    .. note:: For rho>1 the returned value is 0

    .. note:: Values for rho<0 are silently returned as rho=0
    """

    Rnm = rnm(n, m, rho)

    NC = sqrt(2 * (n + 1))
    # S = (n - abs(m)) / 2

    if m > 0:
        Zmn = NC * Rnm * cos(m * theta)
    # las funciones cos() y sin() de scipy tienen problemas cuando la grilla
    # tiene dimension cero

    elif m < 0:
        Zmn = NC * Rnm * sin(m * theta)
    else:
        Zmn = sqrt(0.5) * NC * Rnm
    return Zmn


##
# def Ty_Mon(ypow,xpow,dim=-1):
# n=xpow+ypow
# if dim==-1:
# Ty=N.zeros((n+1,n+1))
# Ty[ypow,xpow]=1
# return Ty
# else:
# Ty=N.zeros((dim+1,dim+1))
# Ty[ypow,xpow]=1
# return Ty
##
##
##
# def Ty_Mat2Vec(PolyMat):
# """
# Retorna un polinomio en forma de vector en la base de los Taylor
# """
# n=PolyMat.shape[0]
# PolyVec=[]
# for i in range(n):
# for j in range(i+1):
# x_pow=i-j
# y_pow=j
# PolyVec.append(PolyMat[y_pow,x_pow])
# return N.asarray(PolyVec)
##
##
##
##
def zernike2taylor(n, m):
    """
    Returns the 2D taylor polynomial, that represents the given zernike
    polynomial

    **ARGUMENTS**

        n,m     n and m orders of the Zernike polynomials

    **RETURN VALUE**
        Poly2D instance containing the polynomial
    """
    TC = zeros((n + 1, n + 1))
    mm = (n - m) / 2
    am = abs(m)
    if m > 0:
        B = mm
        if n % 2 == 0:
            p = 1
            q = (m / 2) - 1
        else:
            p = 1
            q = (m - 1) / 2
    else:
        B = n - mm
        if n % 2 == 0:
            p = 0
            q = -(m / 2)
        else:
            p = 0
            q = -(m + 1) / 2
    for i in range(q + 1):
        for j in range(B + 1):
            for k in range(B - j + 1):
                c = (
                    pow(-1, i + j)
                    * binomial(am, 2 * i + p)
                    * binomial(B - j, k)
                    * int(factorial(n - j))
                    / (
                        int(factorial(j))
                        * int(factorial(mm - j))
                        * int(factorial(n - mm - j))
                    )
                )
                x_pow = 2 * (i + k) + p
                y_pow = n - 2 * (i + j + k) - p
                TC[y_pow, x_pow] = TC[y_pow, x_pow] + c

    n = TC.shape[0] - 1
    cohef = [0.0] * ord2i(n)
    for i in range(ord2i(n)):
        px, py = i2pxpy(i)
        cohef[i] = TC[px, py]
    return Poly2D(cohef)


def i2nm(i):
    """
    Return the n and m orders of the i'th zernike polynomial

    ========= == == == == == == == == == == == == == == == ===
    i          0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 ...
    n-order    0  1  1  2  2  2  3  3  3  3  4  4  4  4  4 ...
    m-order    0 -1  1 -2  0  2 -3 -1  1  3 -4 -2  0  2  4 ...
    ========= == == == == == == == == == == == == == == == ===
    """
    # Calculate the polynomial order
    # Order      0   1   2   3   4
    # initindex  0   1   3   6   10
    ia = array(i)
    n = (1 + (sqrt(8 * (ia) + 1) - 3) / 2).astype(int)
    ni = n * (n + 1) / 2
    m = -n + 2 * (i - ni)
    return n, m


class ZernikeXY(object):
    """
    Class used to evaluate the zernike polinomial with coheficients
    given by cohef, in Cartesian coordinates.


    This class uses an internal Taylor representation of the polynomial
    for all the calculations. The internal taylor representation, is
    generator from the cohef argument given in the constructor.
    If cohef is changed, the internal taylor representation is built
    again.

    **ARGUMENTS**

    cohef -- List containing the Zernike polynomial coheficients.

    The coheficients given in cohef are enumerated as follows

    ========= == == == == == == == == == == == == == == == ===
    index      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 ...
    n-order    0  1  1  2  2  2  3  3  3  3  4  4  4  4  4 ...
    m-order    0 -1  1 -2  0  2 -3 -1  1  3 -4 -2  0  2  4 ...
    ========= == == == == == == == == == == == == == == == ===
    """

    def __init__(self, cohef=[0]):

        self.cohef = cohef

    def __set_cohef__(self, cohef):

        # Save the coheficient list
        self.__cohef__ = cohef

        # Generate the taylor representation of the zernike polinomial

        p = Poly2D([0])

        for i, c in enumerate(cohef):
            n, m = i2nm(i)
            p = p + c * zernike2taylor(n, m)

        self.poly = p

    def __get_cohef__(self):
        return self.__cohef__

    cohef = property(
        __get_cohef__, __set_cohef__, None, "Coefficient list of the zernike polynomial"
    )

    def eval(self):
        """Not implemented yet"""
        pass

    def evalm(self, x, y, mask=True):
        """
        Evaluate the zernike polynomial

        Method used to evaluate the zernike polynomial at the coordinates
        given by the 2D arrays x,y

        **ARGUMENTS**

            ==== ==============================================================
            x    2D array containing the X coordinates where the Zernike
                 polynomial is to be evaluated.
            y    2D array containing the Y coordinates where the Zernike
                 polynomial is to be evaluated.
            mask Flag indicating if the evaluated values are to be masked.
                 If mask is True, the polynomial will only be evaluated at
                 the pupil sqrt(x**2+y**2)<1. , and a masked array is returned.
            ==== ==============================================================

        """

        if mask:
            r = x**2 + y**2
            m = where(r < 1, False, True)
            retval = masked_array(self.poly.meval(x, y), m)
        else:
            retval = self.poly.meval(x, y)
        return retval

    def gpu_eval(self, x, y):
        retval = self.poly.gpu_eval(x, y)
        return retval


##
##
# def Eval_Poly(M,x,y):
# x=N.asarray(x)
# y=N.asarray(y)
# Result=0
# Pows=M.nonzero()
# y_pow=Pows[0]
# x_pow=Pows[1]
# for i in range(y_pow.shape[0]):
# Result=Result+M[y_pow[i],x_pow[i]]*pow(x,x_pow[i])*pow(y,y_pow[i])
# return Result
##
# def Eval_Poly_Comb(Ml,C,x,y):
# M=C[0]*Ml[0]
# for j in range(1,C.shape[0]):
# M=M+C[j]*Ml[j]
# return Eval_Poly(M,x,y)
##
##
##
# def test(N,M):
# R,T=Polar_Array(Rmax=1,DS=0.005)
# P.figure(1)
# i=1
# for n in range(N):
# for m in range(-M,M):
# if abs(m)<=n and (n-abs(m))%2==0:
# P.subplot(N,2*M,i)
# Z=ZK_Poly(n,m,R,T)
# P.imshow(Z)
# P.axis("off")
# i=i+1
# P.savefig("zern.png",dpi=300)
# P.show()
# test(3,3)
