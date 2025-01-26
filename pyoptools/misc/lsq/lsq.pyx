from numpy import array, sqrt, power, ones, empty, arange
from numpy.linalg import solve
cimport cython
from numpy.linalg import solve
from pyoptools.misc.function_2d.poly_2d.poly_2d cimport *
# from pyoptools.misc.Poly2D import *

cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
def polyfit2d(x, y, z, int order=2):
    """
    Perform a 2D polynomial fit using least squares.

    This function fits a polynomial of a specified order to the provided
    two-dimensional data (x, y, z) using a least-squares approach. The fitting
    is done using a Vandermonde matrix, and the coefficients of the polynomial
    are solved to minimize the sum of squared residuals.

    Parameters
    ----------
    x : array-like
        1D array representing the x-coordinates of the data points.
    y : array-like
        1D array representing the y-coordinates of the data points.
    z : array-like
        1D array representing the z-values (dependent variable) at the given
        (x, y) points.
    order : int, optional
        The order of the polynomial to fit (default is 2).

    Returns
    -------
    ret_poly : Poly2D
        A `Poly2D` object representing the fitted polynomial function.
    e : float
        The root mean square error (RMSE) of the fit, calculated as the square
        root of the mean of the squared differences between the observed and
        fitted values.

    Examples
    --------
    >>> x = [1.0, 2.0, 3.0]
    >>> y = [1.0, 2.0, 3.0]
    >>> z = [1.0, 4.0, 9.0]
    >>> ret_poly, error = polyfit2d(x, y, z, order=2)
    >>> print(ret_poly)
    Poly2D object with coefficients: [ ... ]
    >>> print(error)
    0.1234

    Notes
    -----
    The polynomial fit is performed using a Vandermonde matrix construction
    for two-dimensional data. The matrix is formed based on the powers of the
    input coordinates (x and y), and the least-squares solution is obtained
    by solving the corresponding normal equations.
    """
    cdef int nc, nd, p, n, ir, ic, id, powx, powy
    cdef double sum
    cdef double[:, :] XP, YP, mat, vec, cohef
    cdef int[:, :] potx, poty
    cdef int[:] px, py

    # Compute the number of coefficients
    nc = (order + 2) * (order + 1) // 2  # Number of polynomial terms

    # Convert input data to NumPy arrays and memoryviews
    xa = array(x, order="C", dtype="double")
    ya = array(y, order="C", dtype="double")
    za = array(z, order="C", dtype="double")

    # Create memoryviews for the input arrays
    cdef double[:] xa_mv = xa
    cdef double[:] ya_mv = ya
    cdef double[:] za_mv = za

    n = xa.shape[0]

    # Create indices array and memoryview
    indices = arange(n, dtype="int32")
    cdef int[:] indices_mv = indices

    # Compute powers for x and y
    px, py = indices_to_powers(indices_mv)

    nd = len(xa)

    # Create Vandermonde matrices with memoryviews
    XP = ones((nd, 2 * (order + 1)), order="C", dtype="double")
    YP = ones((nd, 2 * (order + 1)), order="C", dtype="double")

    # Create memoryviews for Vandermonde matrices
    cdef double[:, :] XP_mv = XP
    cdef double[:, :] YP_mv = YP

    # Compute Vandermonde matrices for x and y
    for p in range(1, 2 * (order + 1)):
        for id in range(nd):
            XP_mv[id, p] = XP_mv[id, p - 1] * xa_mv[id]
            YP_mv[id, p] = YP_mv[id, p - 1] * ya_mv[id]

    # Calculate X and Y powers
    potx = empty((nc, nc), dtype="int32")
    poty = empty((nc, nc), dtype="int32")

    # Create memoryviews for power matrices
    cdef int[:, :] potx_mv = potx
    cdef int[:, :] poty_mv = poty

    # This was replaced by the loop
    # potx_mv[:, :] = array(px, dtype="int32")[:]
    # poty_mv[:, :] = array(py, dtype="int32")[:]

    for i in range(nc):
        for j in range(nc):
            potx_mv[i, j] = px[j]
            poty_mv[i, j] = py[j]

    # Create temporary matrices for the transpose operation
    cdef int[:, :] potx_t = empty((nc, nc), dtype="int32")
    cdef int[:, :] poty_t = empty((nc, nc), dtype="int32")

    # Create transpose matrices
    for i in range(nc):
        for j in range(nc):
            potx_t[i, j] = px[i]  # This creates the transpose effect
            poty_t[i, j] = py[i]  # This creates the transpose effect

    # The loop replace this as it is not allowed in cython memory views
    # potx_mv[:, :] += potx_mv.T
    # poty_mv[:, :] += poty_mv.T

    # Add the original and transpose matrices
    for ir in range(nc):
        for ic in range(nc):
            potx_mv[ir, ic] = potx_mv[ir, ic] + potx_t[ir, ic]
            poty_mv[ir, ic] = poty_mv[ir, ic] + poty_t[ir, ic]

    # Create Vandermonde matrix
    mat = empty((nc, nc), dtype="double")
    cdef double[:, :] mat_mv = mat

    for ir in range(nc):
        for ic in range(nc):
            powx = potx_mv[ir, ic]
            powy = poty_mv[ir, ic]
            sum = 0.0
            for id in range(nd):
                sum += XP_mv[id, powx] * YP_mv[id, powy]
            mat_mv[ir, ic] = sum

    # Right-hand side vector
    vec = empty((nc, 1), dtype="double")
    cdef double[:, :] vec_mv = vec

    for ic in range(nc):
        sum = 0.0
        powx = px[ic]
        powy = py[ic]
        for id in range(nd):
            sum += za_mv[id] * XP_mv[id, powx] * YP_mv[id, powy]
        vec_mv[ic, 0] = sum

    # Solve for the coefficients
    cohef = solve(mat_mv, vec_mv)
    cdef double[:, :] cohef_mv = cohef
    ret_poly = Poly2D(cohef_mv[:, 0])

    # Calculate error
    e = sqrt(power(za_mv - ret_poly.eval_1d(xa_mv, ya_mv), 2).mean())

    return ret_poly, e


@cython.boundscheck(False)
@cython.wraparound(False)
def vander_matrix(x, y, z, int order=2):
    cdef int nc, nd, p
    cdef int[:] px, py

    nc= (order+2)*(order+1)//2  # _ord2i_(order)
    xa=array(x, order="C", dtype="double")
    ya=array(y, order="C", dtype="double")
    za=array(z, order="C", dtype="double")

    # Create indices array and memoryview
    n = xa.shape[0]
    indices = arange(n, dtype="int32")
    cdef int[:] indices_mv = indices

    # Compute powers for x and y
    px, py = indices_to_powers(indices_mv)

    nd=len(xa)
    XP=ones((nd, 2*(order+1)), order="C", dtype="double")
    YP=ones((nd, 2*(order+1)), order="C", dtype="double")

    for p in range(1, 2*(order+1)):
        XP[:, p]=XP[:, p-1]*xa
        YP[:, p]=YP[:, p-1]*ya

    # Calculating X and Y powers
    potx=empty((nc, nc), dtype="int")
    poty=empty((nc, nc), dtype="int")
    potx[:, :]=px
    potx=potx+potx.T

    poty[:, :]=py
    poty=poty+poty.T

    # Creating the Vandermonde matrix
    mat=empty((nc, nc))
    # The python code was changed to a C++ inline
    # for  ir in range(nc):
    #    for ic in range(nc):
    #        powx=potx[ir, ic]
    #        powy=poty[ir, ic]
    #        tv=XP[:, powx]*YP[:, powy]
    #        mat[ir, ic]=tv.sum()
    cdef int powx, powy, ir, ic, id
    cdef double sum
    for ir in range(nc):
        for ic in range(nc):
            powx=potx[ir, ic]
            powy=poty[ir, ic]
            sum=0
            for id in range(nd):
                sum=sum+XP[id, powx]*YP[id, powy]
            mat[ir, ic]=sum

    # imat=pinv(mat)

    vec=empty((nc, 1))
    # The python code was changed to C++ inline
    # for ic in range (nc):
    #    tv=XP[:, px[ic]]*YP[:, py[ic]]*z
    #    vec[ic, 0]=tv.sum()

    for ic in range(nc):
        sum=0
        powx=px[ic]
        powy=py[ic]
        for id in range(nd):
            sum=sum+za[id]*XP[id, powx]*YP[id, powy]
        vec[ic, 0]=sum
    return mat, vec


# Local definitions for the optimized order 2 fitter

_potxo2=array([[0.,  1.,  0.,  2.,  1.,  0.],
               [1.,  2.,  1.,  3.,  2.,  1.],
               [0.,  1.,  0.,  2.,  1.,  0.],
               [2.,  3.,  2.,  4.,  3.,  2.],
               [1.,  2.,  1.,  3.,  2.,  1.],
               [0.,  1.,  0.,  2.,  1.,  0.]])
_potyo2=array([[0.,  0.,  1.,  0.,  1.,  2.],
               [0.,  0.,  1.,  0.,  1.,  2.],
               [1.,  1.,  2.,  1.,  2.,  3.],
               [0.,  0.,  1.,  0.,  1.,  2.],
               [1.,  1.,  2.,  1.,  2.,  3.],
               [2.,  2.,  3.,  2.,  3.,  4.]])
_pxo2, _pyo2=(array([0, 1, 0, 2, 1, 0]), array([0, 0, 1, 0, 1, 2]))


def polyfito2(x, y, z):
    """
    Polyfit function optimized for polynomials of order 2
    """
    nc= 6  # _ord2i_(order)
    xa=array(x)
    ya=array(y)
    px, py= _pxo2, _pyo2  # _i2pxpy_(range(0, nc))
    nd=len(xa)
    XP=ones((nd, 6))  # ones((nd, 2*(order+1)))
    YP=ones((nd, 6))
    # for p in range(1, 6):
    #    XP[:, p]=XP[:, p-1]*xa
    #    YP[:, p]=YP[:, p-1]*ya
    cdef int p, id
    for p in range(1, 6):
        for id in range(nd):
            XP[id, p]=XP[id, p-1]*xa[id]
            YP[id, p]=YP[id, p-1]*ya[id]

    # Calculating X and Y powers
    # potx=empty((nc, nc))
    # poty=empty((nc, nc))
    # potx[:, :]=px
    # potx=potx+potx.T
    # print "*",potx
    # poty[:, :]=py
    # poty=poty+poty.T
    # print "**", poty
    potx=_potxo2
    poty=_potyo2
    # Creating the Vandermonde matrix
    mat=empty((nc, nc))
    # The python code was changed to a C++ inline
    # for  ir in range(nc):
    #    for ic in range(nc):
    #        powx=potx[ir, ic]
    #        powy=poty[ir, ic]
    #        tv=XP[:, powx]*YP[:, powy]
    #        mat[ir, ic]=tv.sum()
    cdef int powx, powy
    cdef int ir, ic
    cdef double sum
    for ir in range[nc]:
        for ic in range(nc):
            powx=potx[ir, ic]
            powy=poty[ir, ic]
            sum=0
            for id in range(nd):
                sum=sum+XP[id, powx]*YP[id, powy]
            mat[ir, ic]=sum

    # imat=pinv(mat)

    vec=empty((nc, 1))

    for ic in range(nc):
        tv=XP[:, px[ic]]*YP[:, py[ic]]*z
        vec[ic, 0]=tv.sum()

    cohef= solve(mat, vec)
    ret_poly=Poly2D(cohef, px=px, py=py)

    # Calculate error. Verify if this is the best way
    ep=cohef[0]+\
        cohef[1]*xa+\
        cohef[2]*ya+\
        cohef[3]*xa*xa+\
        cohef[4]*xa*ya+\
        cohef[5]*ya*ya

    e= sqrt(power(array(z)-ep, 2).mean())

    return ret_poly, e


# Se elimina este código pues no debe ser importante, y esta usando weave
#
# _potxo1= array([[ 0.,   1.,   0.],
#           [ 1.,   2.,   1.],
#           [ 0.,   1.,   0.]])
#     _potyo1= array([[ 0.,   0.,   1.],
#           [ 0.,   0.,   1.],
#           [ 1.,   1.,   2.]])
#
# _pxo1, _pyo1=(array([0, 1, 0]), array([0, 0, 1]))
#
# def polyfito1(x, y, z):
#     xa=array(x)
#     ya=array(y)
#     px, py=_pxo1, _pyo1
#     nd=len(xa)
#     XP=ones((nd, 4))
#     YP=ones((nd, 4))
#     code="""
#             for(int p=1;p<4;p++)
#             {
#                 for(int id=0;id<nd;id++)
#                 {
#                     XP(id,p)=XP(id,p-1)*xa(id);
#                     YP(id,p)=YP(id,p-1)*ya(id);
#                 }
#             }
#     """
#     err = weave.inline(code,['nd', 'XP','YP','xa', 'ya'],
#                         type_converters=converters.blitz,
#                         compiler = 'gcc')
#
#     #Calculating X and Y powers
#     potx=_potxo1
#     poty=_potyo1
#
#     #Creating the Vandermonde matrix
#     mat=empty((3, 3))
#
#     code = """
#             int powx,powy;
#             for(int ir=0;ir<3;ir++)
#                 for(int ic=0;ic<3;ic++)
#                     {
#                         powx=potx(ir,ic);
#                         powy=poty(ir,ic);
#                         double sum=0;
#                         for(int id=0;id<nd;id++)
#                             sum=sum+XP(id,powx)*YP(id,powy);
#                         mat(ir,ic)=sum;
#                     }
#             """
#     err = weave.inline(code,['potx', 'poty', 'XP', 'YP', 'mat', 'nd'],
#                         type_converters=converters.blitz,
#                         compiler = 'gcc')
#     ######################
#
#
#     #imat=pinv(mat)
#
#     vec=empty((3, 1))
#
#     code="""
#             int powx,powy;
#             for(int ic=0;ic<3;ic++)
#             {
#                 double sum=0;
#                 powx=px(ic);
#                 powy=py(ic);
#                 for(int id=0;id<nd;id++)
#                     sum=sum+z(id)*XP(id,powx)*YP(id,powy);
#                 vec(ic,0)=sum;
#             }
#     """
#     err = weave.inline(code,[ 'nd', 'XP','YP','px', 'py', 'vec', 'z'],
#                         type_converters=converters.blitz,
#                         compiler = 'gcc')
#     cohef= solve(mat, vec)
#     ret_poly=poly2d(cohef)
#
#     #Calculate error. Verify if this is the best way
#     e= sqrt(power(array(z)-ret_poly.eval(xa, ya), 2).mean())
#
#     return ret_poly, e
#
#


# ~ def test(niter=10,Omax=4,  ndat=500):
    # ~ from numpy.random import rand
    # ~ for n in range(niter):
    #    # ~ ncohef=ord2i(Omax)
    #    # ~ cohef=100*rand(ncohef)
    #    # ~ X=[]
    #    # ~ Y=[]
    #    # ~ Z=[]
    #    # ~ p=poly2d(cohef)
    #    # ~ for j in range(ndat):
    #        # ~ x=100.*rand()
    #        # ~ y=100.*rand()
    #        # ~ X.append(x)
    #        # ~ Y.append(y)
    #        # ~ Z.append(p.eval(x, y))
    #    # ~ co, er=polyfit(X, Y, Z, Omax)
    #    # ~ print cohef
    #    # ~ print sqrt(power(cohef-co.cohef.T, 2).sum())/ncohef
    #    # ~ print "*"
    # ~  Test using the open opt library it seems not to be needed
    # ~ def polyfitOO(x, y, z, order):
    # ~
    # ~ def f(C, D):
    #    # ~ x=D[0]
    #    # ~ y=D[1]
    #    # ~ a=poly2d(cohef=C)
    #    # ~ return a.eval(x, y)
    # ~ ND=_ord2i_(order)
    # ~ p, e=polyfit(array(x), array(y), array(z), order=order)
    # ~ C=p.cohef
    # ~ print C
    # ~ D=array([x, y]).transpose()
    # ~ Z=array(z).transpose()
    # ~ lb = ND*[-inf]
    # ~ ub = ND*[inf]
    # ~ p = DFP(f, C , D, Z, lb=lb, ub=ub)
    # ~ r = p.solve('nlp:ralg', plot=0, iprint = 10)
    # ~ #print C
    # ~ #print 'solution: '+str(r.xf)+'\n||residuals||^2 = '+str(r.ff)
    # ~ ret_poly=poly2d(r.xf)
    # ~ xa=array(x)
    # ~ ya=array(y)
    # ~ e= sqrt(power(array(z)-ret_poly.eval(xa, ya), 2).mean())
    # ~ return ret_poly, e
