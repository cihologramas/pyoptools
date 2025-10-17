
import numpy as N
from numpy import (
    array,
    sin,
    cos,
    float64,
    dot,
    sqrt,
    floor,
    meshgrid,
    zeros,
    where,
    pi,
    isnan,
    nonzero,
    rint,
    linspace,
    arange,
    argwhere,
    mean,
)
from numpy.ma import is_masked, MaskedArray
from numpy.ma import array as ma_array

from scipy import interpolate
from matplotlib.tri import Triangulation


"""Auxiliary functions and classes
"""


def rot_x(tx):
    """Returns the transformation matrix for a rotation around the X axis"""
    return array(
        [[1.0, 0.0, 0.0], [0.0, cos(tx), -sin(tx)], [0.0, sin(tx), cos(tx)]]
    ).astype(float64)


def rot_y(ty):
    """Returns the transformation matrix for a rotation around the Y axis"""
    return array(
        [[cos(ty), 0.0, sin(ty)], [0.0, 1, 0.0], [-sin(ty), 0.0, cos(ty)]]
    ).astype(float64)


def rot_z(tz):
    """Returns the transformation matrix for a rotation around the Z axis"""
    return array(
        [[cos(tz), -sin(tz), 0.0], [sin(tz), cos(tz), 0.0], [0.0, 0.0, 1.0]]
    ).astype(float64)


# ~ def rot_mat(r):
# ~ '''Returns the transformation matrix for a rotation around the Z,Y,X axes
# ~
# ~ The rotation is made first around the Z axis, then around the Y axis, and
# ~ finally around the X axis.
# ~
# ~ Parameters
# ~
# ~ r= (rx,ry,rz)
# ~ '''
# ~
# ~ c=cos(r)
# ~ s=sin(r)
# ~
# ~ rx=array([[1. , 0., 0.],
# ~ [0. , c[0],-s[0]],
# ~ [0. , s[0], c[0]]])
# ~
# ~ ry=array([[ c[1], 0., s[1]],
# ~ [ 0., 1., 0.],
# ~ [-s[1], 0., c[1]]])
# ~
# ~
# ~ rz=array([[ c[2],-s[2], 0.],
# ~ [ s[2], c[2], 0.],
# ~ [ 0., 0., 1.]])
# ~
# ~
# ~ tm=dot(rz,dot(ry,rx))
# ~
# ~ return tm

#  To improve speed, this routine was moved to cmisc.pyx
# ~ def rot_mat_i(r):
# ~ '''Returns the inverse transformation matrix for a rotation around the Z,Y,X axes
# ~
# ~ Parameters
# ~
# ~ r= (rx,ry,rz)
# ~ '''
# ~
# ~ c=cos(r)
# ~ s=sin(r)
# ~
# ~ rx=array([[ 1., 0., 0.],
# ~ [ 0., c[0], s[0]],
# ~ [ 0.,-s[0], c[0]]])
# ~
# ~ ry=array([[ c[1], 0.,-s[1]],
# ~ [ 0., 1., 0.],
# ~ [ s[1], 0., c[1]]])
# ~
# ~
# ~ rz=array([[ c[2], s[2], 0.],
# ~ [-s[2], c[2], 0.],
# ~ [ 0., 0., 1.]])
# ~
# ~ # Nota: se hizo una prueba para optimizar escribirndo la expresión del producto
# ~ # escalar, y el resultado fue considerablemente mas lento, toca revisar
# ~
# ~
# ~ return dot(rx,dot(ry,rz))


def cross(a, b):
    """3D Vector product producto vectorial"""
    x1, y1, z1 = a
    x2, y2, z2 = b
    return array((y1 * z2 - y2 * z1, x2 * z1 - x1 * z2, x1 * y2 - x2 * y1))


def wavelength2RGB(wl):
    """Function to approximate and RGB tuple from the wavelength value

    Parameter:

    wavelength wavelength in um

    if the wavelength is outside the visible spectrum returns (0,0,0)
    Original code found at:

    http://www.physics.sfasu.edu/astro/color/spectra.html

    """

    R, G, B = 0.0, 0.0, 0.0

    if (wl >= 0.380) & (wl < 0.440):
        R = -1.0 * (wl - 0.440) / (0.440 - 0.380)
        G = 0.0
        B = 1.0

    if (wl >= 0.440) & (wl < 0.490):
        R = 0.0
        G = (wl - 0.440) / (0.490 - 0.440)
        B = 1.0

    if (wl >= 0.490) & (wl < 0.510):
        R = 0.0
        G = 1.0
        B = -1.0 * (wl - 0.510) / (0.510 - 0.490)

    if (wl >= 0.510) & (wl < 0.580):
        R = (wl - 0.510) / (0.580 - 0.510)
        G = 1.0
        B = 0.0
    if (wl >= 0.580) & (wl < 0.645):
        R = 1.0
        G = -1.0 * (wl - 0.645) / (0.645 - 0.580)
        B = 0.0
    if (wl >= 0.645) & (wl < 0.780):
        R = 1.0
        G = 0.0
        B = 0.0
    # LET THE INTENSITY FALL OFF NEAR THE VISION LIMITS

    if wl >= 0.700:
        sss = 0.3 + 0.7 * (0.780 - wl) / (0.780 - 0.700)
    elif wl < 0.420:
        sss = 0.3 + 0.7 * (wl - 0.380) / (0.420 - 0.380)
    else:
        sss = 1

    R = R * sss
    G = G * sss
    B = B * sss
    return (R, G, B)


def matrix_interpolation(M, i, j, type="bilinear"):
    """Returns the interpolated value of a matrix, when the indices i,j are
    floating    point numbers.

    **ARGUMENTS**
        ==== =====================================================
        M    Matrix to interpolate
        i,j  Indices to interpolate
        type Interpolation type. supported types: nearest,bilinear
        ==== =====================================================

    """
    mi, mj = M.shape
    if i < 0 or i > mi - 2 or j < 0 or j > mj - 2:
        raise IndexError("matrix Indexes out of range")
    # Allowed interpolation types
    inter_types = [
        "nearest",
        "bilinear",
    ]
    if type not in inter_types:
        raise ValueError(
            "Interpolation type not allowed. The allowed types"
            " are: {0}".format(inter_types)
        )
    if type == "nearest":
        iri = int(round(i))
        irj = int(round(j))
        return M[iri, irj]
    elif type == "bilinear":
        i_s, j_s = floor((i, j))
        # calc 1
        m = M[i_s:i_s + 2, j_s:j_s + 2]
        iv = array([1 - (i - i_s), i - i_s])
        jv = array(
            [
                [
                    1 - (j - j_s),
                ],
                [
                    j - j_s,
                ],
            ]
        )
        return dot(iv, dot(m, jv))[0]
        # dx=i-i_s
        # dy=j-j_s
        # print i, j, i_s, j_s,  dx, dy
        # p1=dx*dy*M[i_s, j_s]
        # p2=(1.-dx)*dy*M[i_s+1, j_s]
        # p3=dx*(1.-dy)*M[i_s, j_s+1]
        # p4=(1.-dx)*(1.-dy)*M[i_s+1, j_s+1]
        # return p1+ p2+ p3+ p4
    print("error")
    return 1.0


def hitlist2int(x, y, z, xi, yi):
    """Function that estimates an intensity distribution on a plane from a
    ray hitlist
    """
    # if xi.ndim != yi.ndim:
    #    raise TypeError("inputs xi and yi must have same number of dimensions (1 or 2)")
    # if xi.ndim != 1 and xi.ndim != 2:
    #    raise TypeError("inputs xi and yi must be 1D or 2D.")
    # if not len(x)==len(y)==len(z):
    #    raise TypeError("inputs x,y,z must all be 1D arrays of the same length")
    # remove masked points.
    # if hasattr(z,'mask'):
    #    x = x.compress(z.mask == False)
    #    y = y.compress(z.mask == False)
    #    z = z.compressed()

    # if xi.ndim == 1:
    #    xi,yi = meshgrid(xi,yi)

    # triangulate data

    tri = Triangulation(x, y)

    # calculate triangles area
    # ntriangles=tri.circumcenters.shape[0]
    coord = array(zip(tri.x, tri.y))

    # I=zeros((ntriangles, ))
    # xc=zeros((ntriangles, ))
    # yc=zeros((ntriangles, ))
    # for i in range(ntriangles):
    #     i1, i2, i3=tri.triangle_nodes[i]
    #     p1=coord[i1]
    #     p2=coord[i2]
    #     p3=coord[i3]
    #     v1=p1-p2
    #     v2=p3-p2
    #     I[i]=1./(abs(v1[0]*v2[1]-v1[1]*v2[0]))
    #     # the circumcenter data from the triangulation, has some problems so we
    #     # recalculate it
    #     xc[i], yc[i]=(p1+p2+p3)/3.
    # The previous code was replaced by the following code
    ###
    i1 = tri.triangles[:, 0]
    i2 = tri.triangles[:, 1]
    i3 = tri.triangles[:, 2]
    p1 = coord[i1]
    p2 = coord[i2]
    p3 = coord[i3]

    v1 = p1 - p2
    v2 = p3 - p2
    intensity = abs(1.0 / (v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]))

    c = (p1 + p2 + p3) / 3.0
    xc = c[:, 0]
    yc = c[:, 1]
    ###

    # Because of the triangulation algorithm, there are some really high values
    # in the intensity data. To filter these values, remove the 5% points of the
    # higher intensity.
    ni = int(0.1 * len(intensity))
    j = intensity.argsort()[:-ni]
    xc = xc[j]
    yc = yc[j]
    intensity = intensity[j]
    intensity = intensity / intensity.max()

    #    #print tri.circumcenters[:, 0]
    #    #print tri.circumcenters.shape
    #    print ntriangles,  tri.circumcenters[:, 0].shape,  tri.circumcenters[:, 0].flatten().shape

    # itri=delaunay.Triangulation(xc,yc)
    # inti=itri.linear_interpolator(I)
    # xi,yi = meshgrid(xi,yi)
    # d1=itri(xi, yi)

    # Interpolacion con Splines
    # di=interpolate.SmoothBivariateSpline(xc, yc, I)
    # d1=di(xi,yi)

    # Interpolacion nn, y generación de pupila
    xi, yi = meshgrid(xi, yi)
    d1 = interpolate.griddata(xc, yc, intensity, xi, yi)

    return d1


def hitlist2int_list(x, y):
    """Function that estimates an intensity distribution on a plane from a
    ray hitlist. Returns the intensity samples as an x,y,I list
    """

    # if xi.ndim != yi.ndim:
    #    raise TypeError("inputs xi and yi must have same number of dimensions (1 or 2)")
    # if xi.ndim != 1 and xi.ndim != 2:
    #    raise TypeError("inputs xi and yi must be 1D or 2D.")
    # if not len(x)==len(y)==len(z):
    #    raise TypeError("inputs x,y,z must all be 1D arrays of the same length")
    # remove masked points.
    # if hasattr(z,'mask'):
    #    x = x.compress(z.mask == False)
    #    y = y.compress(z.mask == False)
    #    z = z.compressed()

    # if xi.ndim == 1:
    #    xi,yi = meshgrid(xi,yi)

    # triangulate data

    tri = Triangulation(x, y)

    # calculate triangles area
    # ntriangles=tri.circumcenters.shape[0]
    coord = array(zip(tri.x, tri.y))

    # I=zeros((ntriangles, ))
    # xc=zeros((ntriangles, ))
    # yc=zeros((ntriangles, ))
    # for i in range(ntriangles):
    #     i1, i2, i3=tri.triangle_nodes[i]
    #     p1=coord[i1]
    #     p2=coord[i2]
    #     p3=coord[i3]
    #     v1=p1-p2
    #     v2=p3-p2
    #     I[i]=1./(abs(v1[0]*v2[1]-v1[1]*v2[0]))
    #     # the circumcenter data from the triangulation, has some problems so we
    #     # recalculate it
    #     xc[i], yc[i]=(p1+p2+p3)/3.
    # The previous code was replaced by the following code
    ###
    i1 = tri.triangles[:, 0]
    i2 = tri.triangles[:, 1]
    i3 = tri.triangles[:, 2]
    p1 = coord[i1]
    p2 = coord[i2]
    p3 = coord[i3]

    v1 = p1 - p2
    v2 = p3 - p2
    intensity = abs(1.0 / (v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0]))

    c = (p1 + p2 + p3) / 3.0
    xc = c[:, 0]
    yc = c[:, 1]
    ###

    # Because of the triangulation algorithm, there are some really high values
    # in the intensity data. To filter these values, remove the 5% points of the
    # higher intensity.
    ni = int(0.1 * len(intensity))
    j = intensity.argsort()[:-ni]
    xc = xc[j]
    yc = yc[j]
    intensity = intensity[j]
    intensity = intensity / intensity.max()

    #    #print tri.circumcenters[:, 0]
    #    #print tri.circumcenters.shape
    #    print ntriangles,  tri.circumcenters[:, 0].shape,  tri.circumcenters[:, 0].flatten().shape

    # itri=delaunay.Triangulation(xc,yc)
    # inti=itri.linear_interpolator(I)
    # xi,yi = meshgrid(xi,yi)
    # d1=itri(xi, yi)

    # Interpolacion con Splines
    # di=interpolate.SmoothBivariateSpline(xc, yc, I)
    # d1=di(xi,yi)

    return xc, yc, intensity


def unwrapv(inph, in_p=(), uv=2 * pi):
    """Return the input matrix unwrapped the value given in uv

    This is a vectorized routine, but is not as fast as it should
    """

    if not is_masked(inph):
        fasei = MaskedArray(inph, isnan(inph))
    else:
        fasei = inph.copy()

    size = fasei.shape
    nx, ny = size
    # If the initial unwraping point is not given, take the center of the image
    # as initial coordinate
    if in_p == ():
        in_p = (int(size[0] / 2), int(size[1] / 2))

    # Create a temporal space to mark if the points are already unwrapped
    # 0 the point has not been unwrapped
    # 1 the point has not been unwrapped, but it is in the unwrapping list
    # 2 the point was already unwrapped

    fl = N.zeros(size)

    # List containing the points to unwrap
    l_un = [in_p]
    fl[in_p] = 1

    # unwrapped values
    faseo = fasei.copy()
    XI_, YI_ = meshgrid(range(-1, 2), range(-1, 2))
    XI_ = XI_.flatten()
    YI_ = YI_.flatten()
    while len(l_un) > 0:
        # remove the first value from the list
        unp = l_un.pop(0)
        # l_un[0:1]=[]
        XI = XI_ + unp[0]
        YI = YI_ + unp[1]
        # Remove from the list the values where XI is negative
        nxi = XI > -1
        nyi = YI > -1
        nxf = XI < nx
        nyf = YI < ny
        n = nonzero(nxi & nyi & nxf & nyf)
        lco = zip(XI[n], YI[n])

        # Put the coordinates of unwrapped the neighbors in the list

        # And check for wrapping
        nv = 0
        wv = 0

        for co in lco:
            if (fl[co] == 0) & (faseo.mask[co] is False):
                fl[co] = 1
                l_un.append(co)
            elif fl[co] == 2:
                wv = wv + rint((faseo[co] - faseo[unp]) / uv)
                nv = nv + 1

        if nv != 0:
            wv = wv / nv
            # if wv>=0: wv=int(wv+0.5)
            # else: wv=int(wv-0.5)
        fl[unp] = 2
        faseo[unp] = faseo[unp] + wv * uv

    return faseo


def unwrap_py(inph, in_p=(), uv=2 * pi):
    """Return the input matrix unwrapped the value given in uv

    The same as unwrapv, but using for-s, written in python
    """
    if not is_masked(inph):
        fasei = MaskedArray(inph, isnan(inph))
    else:
        fasei = inph

    nx, ny = (fasei.shape[0], fasei.shape[1])

    # If the initial unwraping point is not given, take the center of the image
    # as initial coordinate
    if in_p == ():
        in_p = (int(nx / 2), int(ny / 2))

    # Create a temporal space to mark if the points are already unwrapped
    # 0 the point has not been unwrapped
    # 1 the point has not been unwrapped, but it is in the unwrapping list
    # 2 the point was already unwrapped

    fl = zeros((nx, ny))

    # List containing the points to unwrap
    l_un = [in_p]
    fl[in_p] = 1

    # unwrapped values
    faseo = fasei.copy()

    while len(l_un) > 0:
        # remove the first value from the list
        cx, cy = l_un.pop(0)

        # Put the coordinates of unwrapped the neighbors in the list
        # And check for wrapping
        nv = 0
        wv = 0

        for i in range(cx - 1, cx + 2):
            for j in range(cy - 1, cy + 2):
                if (i > -1) and (i < nx) and (j > -1) and (j < ny):
                    if (fl[i, j] == 0) & (faseo.mask[i, j] is False):
                        fl[i, j] = 1
                        l_un.append((i, j))
                    elif fl[i, j] == 2:
                        wv = wv + rint((faseo[i, j] - faseo[cx, cy]) / uv)
                        nv = nv + 1
        if nv != 0:
            wv = wv / nv

        fl[cx, cy] = 2
        faseo[cx, cy] = faseo[cx, cy] + wv * uv

    return faseo


def interpolate_g(xi, yi, zi, xx, yy, knots=10, error=False, mask=None):
    """Create a grid of zi values interpolating the values from xi,yi,zi

    **ARGUMENTS**

        ========= ==================================================================
        xi,yi,zi  1D Lists  or arrays containing the values to use as base for the
                  interpolation
        xx,yy     1D vectors or lists containing the output coordinates
        samples   tuple containing the shape of the output array.
        knots     number of knots to be used in each direction
        error     if set to true, half of the points (x, y, z) are used to create
                  the interpolation, and half are used to evaluate the interpolation
                  error
        ========= ==================================================================


    """
    xi = array(xi)
    yi = array(yi)
    zi = array(zi)

    # print xi
    # print yi
    # print zi
    assert xi.ndim == 1, "xi must ba a 1D array or list"
    assert yi.ndim == 1, "yi must ba a 1D array or list"
    assert zi.ndim == 1, "zi must ba a 1D array or list"

    assert xx.ndim == 1, "xx must ba a 1D array or list"
    assert yy.ndim == 1, "yy must ba a 1D array or list"

    assert len(xi) == len(yi) and len(xi) == len(
        zi
    ), "xi, yi, zi must have the same number of items"

    if error:
        # Create a list of indexes to be able to select the points that are going
        # to be used as spline generators, and as control points
        idx = where(arange(len(xi)) % 2 == 0, False, True)

    # Use only half of the samples to create the Spline,
    if error is True:
        isp = argwhere(idx)
        ich = argwhere(not idx)

        xsp = xi[isp]
        ysp = yi[isp]
        zsp = zi[isp]

        xch = xi[ich]
        ych = yi[ich]
        zch = zi[ich]
    else:
        xsp = xi
        ysp = yi
        zsp = zi

    # Distribute homogeneously the knots
    xk = linspace(xsp.min(), xsp.max(), knots)
    yk = linspace(ysp.min(), ysp.max(), knots)

    # LSQBivariateSpline using some knots gives smaller error than
    # SmoothBivariateSpline
    di = interpolate.LSQBivariateSpline(xsp, ysp, zsp, xk[1:-1], yk[1:-1])
    # print xsp,ysp,zsp
    # di=interpolate.SmoothBivariateSpline(xsp, ysp, zsp)

    # Evaluate error
    if error:
        zch1 = di.ev(xch, ych)
        er = (zch.flatten() - zch1).std()

    if mask is None:
        # d=griddata(xi,  yi,  zi,  xx, yy) #
        d = di(xx, yy).transpose()
    else:
        d = ma_array(di(xx, yy).transpose(), mask=mask)

    if error:
        return d, er
    else:
        return d


def spot_info(C):
    """
    Function that gives information about the average radius of the rays hitting
    a CCD.

    Args:
        C (:class:`~pyoptools.raytrace.comp_lib.CCD`): Instance of the CCD to
            gather information from.

    Returns:
       tuple: (Rm, (Xm, Ym), (Xmm , Ymm), Rmax )

       Where
           - float: Rm is the mean Radius of the spot diagram
           - tuple: (Xm, Ym) are the mean X and Y values of the spot diagram
           - tuple: (Xmm, Ymm) are the mean radius in X and in Y
           - float: Rmax is the maximum radius

       All the values are measured from the central coordinate (Xm, Ym) of the
       spot diagram.
    """
    X = []
    Y = []
    for hl in C.hit_list:
        x, y, z = hl[0]
        X.append(x)
        Y.append(y)
    xm = mean(X)
    ym = mean(Y)
    X = array(X) - xm
    Y = array(Y) - ym
    R = sqrt(X**2 + Y**2)
    XR = sqrt(X**2)
    YR = sqrt(Y**2)
    return mean(R), (xm, ym), (mean(XR), mean(YR)), R.max()


# Fin Funciones auxiliares
