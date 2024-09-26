# standard imports

# third-party imports
import numpy as np
from numpy.fft import fft, ifft
from numpy import zeros, pi, indices, exp, sqrt

# local imports


# Fractional Fourier Transform


# TODO better docstring
def frft(x, alpha):
    assert x.ndim == 1, "x must be a 1 dimensional array"

    x_expanded_to_2d = np.expand_dims(x, axis=1)
    result_2d = _frft2(x_expanded_to_2d, alpha)
    result_1d = np.squeeze(result_2d)

    return result_1d


def _frft2(x, alpha):
    assert x.ndim == 2, "x must be a 2 dimensional array"
    m, n = x.shape
    # TODO please remove this confusing comment. Is it 'm' or 'm-1' ?
    # TODO If 'p = m', more code cleaning is easy to do.
    p = m  # m-1 # deveria incrementarse el sigiente pow
    y = zeros((2 * p, n), dtype=complex)
    z = zeros((2 * p, n), dtype=complex)

    j = indices(z.shape)[0]
    y[(p - m) // 2:(p + m) // 2, :] = x * exp(
        -1.0j * pi * (j[0:m] ** 2) * float(alpha) / m
    )

    z[0:m, :] = exp(1.0j * pi * (j[0:m] ** 2) * float(alpha) / m)
    z[-m:, :] = exp(1.0j * pi * ((j[-m:] - 2 * p) ** 2) * float(alpha) / m)

    d = exp(-1.0j * pi * j**2 ** float(alpha) / m) * ifft(
        fft(y, axis=0) * fft(z, axis=0), axis=0
    )

    return d[0:m]


# TODO better docstring
def frft2(x, alpha):
    x1 = _frft2(x, alpha)
    return _frft2(x1.transpose(), alpha).transpose()


def rs_kernel(x=0.0, y=0.0, z=0.0, n=1.0):
    """Calculate the Rayleigh Sommerfeld propagation Kernel, for a source point
    at the origin, and a observation point at (x,y,z)
    """
    wavelength = 0.442e-3
    k = 2.0 * pi * n / wavelength
    R2 = x**2 + y**2 + z**2
    R = sqrt(R2)
    ikR = 1.0j * k * R
    R3 = R2 * R
    return (z / (2 * pi * R3)) * exp(ikR) * (1.0 - ikR)


# ~ nd=1024
# ~ data=zeros((nd,nd))
# ~ csize=10
# ~ data[nd/2-csize:nd/2+csize,nd/2-csize:nd/2+csize]=1
# ~
# ~ nx,ny=data.shape
# ~
# ~ dx,dy=0.5e-3,0.5e-3
# ~
# ~ X,Y=indices((nx,ny))
# ~
# ~ ux=(X-nx/2)*dx; uy=(Y-ny/2)*dy #mirar cual es el tamano adecuado para hacer el calculo
# ~
# ~ rs=rs_kernel(x=ux, y=uy, z=250., n=1.)
# ~
# ~
# ~
# ~ RS=fft2(fftshift(rs))
# ~ U=frft2(data,1)
# ~ #U=fft2(data)
# ~
# ~ odata=frft2(RS*U,-1)
# ~ #odata=ifft2(RS*U)
# ~ #figure()
# ~ #imshow(data)
# ~ figure()
# ~ imshow(abs(odata));colorbar()
# ~
# ~
