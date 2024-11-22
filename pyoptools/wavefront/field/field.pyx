# ------------------------------------------------------------------------------
# Copyright (c) 2009, Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amezquita Orozco
# Description:     Field definition module
# Symbols Defined: Field
# ------------------------------------------------------------------------------
#

from warnings import warn

from numpy import dot, zeros, abs, meshgrid, pi, exp, where, angle, sqrt as npsqrt, \
    indices, empty, complex128, complex64, newaxis, column_stack, max, array, \
    linspace, dot, zeros_like, float64, arange, isnan, angle,  mgrid, rint, \
    float32, uint32, exp


from numpy.fft import fft2, ifft2, fftshift, ifftshift
from numpy.ma import array as maarray,  getmask,  getmaskarray
# deprecated use scipy.interpolate.griddata
# from matplotlib.mlab import griddata

from scipy.integrate import simpson as simps
from scipy.ndimage import map_coordinates
# deprecated in scipy 1.2.1; change to imageio
# from scipy.misc import imread
from imageio import imread


# import pp
# import pyopencl as cl


from pyoptools.misc.cmisc.cmisc cimport rot_mat

from pyoptools.raytrace.ray import Ray


# from pyoptools.misc.cmisc.cmisc import unwrap


cdef extern from "math.h":
    double M_PI
    double sqrt(double)

cimport cython


# Nota, los x y los ys parecen estar invertidos al graficar
# Auxiliary functions used to calculate the rotation coheficients


def s(u, n):
    # u=u+1.
    rv=(1.-exp(2.j*pi*u))/(1.-exp(2.j*pi*u/n))

    # insert the result value there u/n = 1
    rv=where(isnan(rv), n, rv)
    return rv


def s2d(kx, ky, nx, ny):
    ux=kx - (arange(nx).astype(float64))  # -int(nx/2))
    uy=ky - (arange(ny).astype(float64))  # -int(ny/2))
    ux= ux.reshape(1, nx)
    uy=uy.reshape(1, ny)
    #    print ux.shape, uy.shapefield.py
    sx=s(ux, nx)
    sy=s(uy, ny)

    # print sx.shape, sy.shape
    # TODO: Verify if the correct expression is
    # dot(sx.transpose(),sy) or dot(sy.transpose(),sx)
    # retval=dot(sx.transpose(),sy)
    retval=dot(sy.transpose(), sx)
    return retval


cdef class Field:
    """
    Class that defines an optical field distribution.

    The sampled field uses an regular Cartesian sampling.

    **ARGUMENTS:**

        ====== =====================================================
        data   2D numpy array containing the complex field information.
               The number of samples are given by data's shape.
        psize  Pixel size of the sample field (resolution)
        amp_im Filename of the image containing the amplitude of the
               field
        ph_im  Filename of the image containing the phase of the field
        amp_n  Floating point number used to normalize the amplitude
               1/255. by default
        ph_n   Floating point number used to normalize the phase.
               2*pi/255 by default.
        ====== =====================================================

    The field data can be given as a numpy complex array, or as a set of
    2 images, one containing the amplitude of the field, and the other
    containing the phase of the field. The amp_n, and the ph_n attributes
    are used to normalize the field. If the images are color images, they
    are flattened to produce a grayscale image. See: scipy.misc.imread
    """
    cdef public object data, psize
    cdef public double l

    def __init__(self, data=None, psize=(10e-3, 10e-3), l=442e-6,
                 amp_im=None, ph_im=None, amp_n=1/255., ph_n=2*pi/255.):

        # Comparing data with None does not work, It seems that cython has a bug
        # when data is a numpy arra and it is compared with none
        assert (not(data is None)) ^ (not(amp_im is None) or not(ph_im is None)), \
            "The initizilation can be a numpy array, or a set of "+\
            "images but not both"

        # Array containing the complex field distribution.

        if data is None:
            if amp_im is not None:
                amp=imread(amp_im, True)
            else:
                amp=1.

            if ph_im is not None:
                ph=imread(ph_im, True)
            else:
                ph=0.
            self.data=amp_n*amp*exp(1.j*ph*ph_n)
        else:
            self.data=data
        # Size of each pixel (sample). It can be given as a tuple psize=(h,w), or as a
        # floating point number psize=h=w.
        self.psize=psize

        # Wave length of the field
        self.l=l

    # Evaluate the possibility to use a field representation different to spatial
    # Maybe an angular spectral representation

    property res:
        '''Returns a tuple (dx,dy) containing the resolution(size of each pixel)
        '''

        def __get__(self):
            try:
                dx, dy=self.psize
            except ValueError:
                dx=self.psize
                dy=dx
            return dx, dy

    property size:
        '''Returns a tuple containing the size of the area where the
        field is sampled
        '''

        def __get__(self):
            dx, dy=self.res
            nx, ny=self.shape
            return array((nx*dx, ny*dy))

    property shape:
        '''Returns the shape of the data contained in the field'''

        def __get__(self):
            return self.data.shape

    property field_sample_coord:
        '''
        Returns the 2D arrays X and Y containing the field sampling
        coordinates
        '''

        def __get__(self):
            dx, dy=self.res
            nx, ny=self.shape
            X, Y=indices((nx, ny))
            ux=(X-nx/2)*dx
            uy=(Y-ny/2)*dy
            return ux, uy

    property xsamples:
        '''Return the sample points in X'''

        def __get__(self):
            dx, _dy=self.res
            nx, _ny=self.shape
            # TODO: Check if this can be calculated using linspaces
            return (arange(nx)-nx/2)*dx

    property ysamples:
        '''Return the sample points in Y'''

        def __get__(self):
            _dx, dy=self.res
            _nx, ny=self.shape
            # TODO: Check if this can be calculated using linspace
            return (arange(ny)-ny/2)*dy

    property phase:
        '''Returns an array containing the unwrapped phase of the optical field'''

        def __get__(self):
            print("Warning: The phase is not being cached. Need to fix this")
            a=self.angle
            return a  # NEED TO FIX UNWRAP AND UNCOMMENT THIS unwrap(a)

    property angle:
        ''' Returns the wrapped phase of the field (mod 2 pi)'''

        def __get__(self):
            a= maarray(angle(self.data),  mask=getmask(self.data))
            return a

    property mask:
        '''Return the mask of the data. True indicates masked values (invalid), False
        indicates valid data.
        '''

        def __get__(self):
            return getmaskarray(self.data)

    def abs(self):
        '''Returns an array containing the absolute value of the field
        '''
        return abs(self.data)

    def intensity(self):
        '''Returns an array containing the intensity of the field
        '''
        # TODO: Check if this is correct
        return abs(self.data)**2

    def conj(self):
        '''Returns the conjugated field
        '''
        return Field(data=self.data.conj(), psize=self.psize, l=self.l)

    # Definition of arithmetic functions for fields
    def __add__(self, a):
        assert isinstance(a, Field), "Error: can not add a Field with a %s" %[
                          type(a)]
        assert self.data.shape==a.data.shape, \
               "Error: both fields must have the same shape"
        assert self.psize==a.psize , "Error: both fields must have the samepixel size"
        assert self.l==a.l , "Error: Both fields must have the same wavelength"
        return Field(data=self.data+a.data, psize=self.psize, l=self.l)

    def __sub__(self, a):
        assert isinstance(a, Field), "Error: can not add a Field with a %s" %[
                          type(a)]
        assert self.data.shape==a.data.shape, \
               "Error: both fields must have the same shape"
        assert self.psize==a.psize , "Error: both fields must have the samepixel size"
        assert self.l==a.l , "Error: Both fields must have the same wavelength"

        return Field(data=self.data-a.data, psize=self.psize, l=self.l)

    def __mul__(self, other):
        # Note, this is different than standard python, here mul==rmul.
        # We need to see what is going on
        if isinstance(self, Field):
            f=self
            a=other
        else:
            f=other
            a=self

        # If you are multipliying 2 fields (1 field* 1 phase mask)
        if isinstance(a, Field):
            # TODO: Find a good way to do the assertions
            # assert(a.l==f.l),"The field and the mask must be for the same wavelength"
            # assert(a.psize==f.psize),"The pixel size in the field and the mask must "
            #                          "be the same"
            return Field(data=a.data*f.data, psize=f.psize, l=f.l)

        # If you are multipliying a field by a numpy array or a number
        return Field(data=a*f.data, psize=f.psize, l=f.l)

    # End of arithmetic functions

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def _rs_kernel(self, double[:, :]x,
                   double[:, :] y, double z, double n):
        """Calculate the Rayleigh Sommerfeld propagation Kernel, for a source point
        at the origin, and a observation point at (x,y,z)
        """
        cdef double k=2.*M_PI*n/self.l
        cdef int nx, ny, i, j
        cdef double xd, yd, R2, R, R3
        # cdef double complex ikR
        cdef complex ikR
        nx=x.shape[0]
        ny=x.shape[1]
        cdef double[:, :] result=zeros((nx, ny), dtype=complex64)
        for i in range(nx):
            for j in range(ny):
                xd=x[i, j]
                yd=y[i, j]
                R2=xd*xd+yd*yd+z*z
                R=sqrt(R2)
                R3=R2*R
                ikR=1j*k*R
                result[i, j]=((z/(2.*M_PI*R3))*exp(ikR)*(1.-ikR))

        return result

    def propagate_rsc_sc(self, z, scale=(1, 1)):
        """Propagate the field a distance z , scaling the output shape
        by the values given in scale. The resolution remains the same.
        """
        cdef int kx, ky, nx, ny, i, j
        cdef double dsx, dsy, dx, dy
        kx=scale[0]
        ky=scale[1]
        dx=self.res[0]
        dy=self.res[1]
        nx=self.shape[0]
        ny=self.shape[1]

        pd=zeros((kx*nx, ky*ny), dtype=complex128)

        for i in range(kx):
            for j in range(ky):
                dsx=(-(<double>kx)/2*nx+nx/2+i*nx)*dx
                dsy=(-(<double>ky)/2*ny+ny/2+j*ny)*dy
                rfd=self.propagate_rsc_d((dsx, dsy, z), n=1.)
                pd[i*nx:(i+1)*nx, j*ny:(j+1)*ny]=rfd.data

        return Field(data=pd, psize=self.psize, l=self.l)

    def propagate_rsc(self, z, n=1., shape=None):
        """Propagate the field a distance z making the convolution with the
        Rayleigh Sommerfeld kernel.

        To avoid noise, the kernel is calculated increasing the size of the matrix.
        """
        if shape is None:
            nx, ny=self.data.shape
        else:
            nx, ny=shape

        dx, dy=self.res

        X, Y=indices((2*nx, 2*ny))
        ux=(X-nx)*dx
        uy=(Y-ny)*dy
        # del X
        # del Y
        krn=self._rs_kernel(ux, uy, z, n)
        # del ux
        # del uy
        krn=fftshift(krn)  # TODO: Check if it is fftshift or ifftshift
        fkrn=fft2(krn)
        # del(krn)
        nf=self.resize((2*nx, 2*ny))
        # ux,uy=self.field_sample_coord
        retval=ifft2(fft2(nf.data)*fkrn)
        # fftconvolve(self.data, self._rs_kernel(x=ux,y=uy,z=z,n=n), mode='same')
        retval=retval[nx/2:nx/2+nx, ny/2:ny/2+ny]*dx*dy
        return Field(data=retval, psize=self.psize, l=self.l)

    def propagate_rsc_d(self, d=(0., 0., 0.), n=1, shape=None):
        """Propagate the field a distance z making the convolution with the
        Rayleigh Sommerfeld kernel.

        To avoid noise, the kernel is calculated increasing the size of the matrix.
        """
        if shape is None:
            nx, ny=self.data.shape
        else:
            nx, ny=shape

        dx, dy=self.res

        x, y, z=d

        X, Y=indices((2*nx, 2*ny))
        ux=(X-nx)*dx+x
        uy=(Y-ny)*dy+y
        # del X
        # del Y
        krn=self._rs_kernel(ux, uy, z, n)
        # del ux
        # del uy
        krn=fftshift(krn)  # TODO: Check if it is fftshift or ifftshift
        fkrn=fft2(krn)
        # del(krn)
        nf=self.resize((2*nx, 2*ny))
        # ux,uy=self.field_sample_coord
        retval=ifft2(fft2(nf.data)*fkrn)
        # fftconvolve(self.data, self._rs_kernel(x=ux,y=uy,z=z,n=n), mode='same')
        retval=retval[nx/2:nx/2+nx, ny/2:ny/2+ny]*dx*dy
        return Field(data=retval, psize=self.psize, l=self.l)

    def _ursi(self, x=0., y=0., z=0., n=1.):
        """Calculate the field at x,y,z , using the Rayleigh Sommerfeld integral
        """
        krs=self._rs_kernel(x, y, z, n)
        # krs=self._rs_kernel(ux+x,uy+y,z,n=n)
        dx, dy=self.res

        # TODO: Check if the samplings dx,dy are in the right place

        # TODO: Check if it is possible to do the integral using the method in
        # Calculation of the Rayleigh–Sommerfeld diffraction integral by exact
        # integration of the fast oscillating factor
        return simps(simps(krs*self.data, dx=dx), dx=dy)

    def propagate_ae(self, z, n=1, shape=None):
        """
        Propagate the field a distance z using angular spectrum method.

        .. warning::
            Check if n is implemented and working
        """

        if shape is None:
            nx, ny=self.data.shape
        else:
            nx, ny=shape

        dx, dy=self.res
        sx=dx*nx
        sy=dy*ny
        # Phase function calculation
        X, Y=indices((nx, ny))
        X=X-nx/2
        Y=Y-ny/2

        X=ifftshift(X)
        Y=ifftshift(Y)

        # Check if the fftshift was OK
        assert(X[0, 0]==0)
        assert(Y[0, 0]==0)

        sq=(1./self.l)**2 -(X/sx)**2-(Y/sy)**2
        # TODO: Check what does this condition mean
        ph=exp(2.*pi*1.j*z*npsqrt(sq))
        ph=where(sq<0, 0., ph)

        pf=ifft2(fft2(self.data, (nx, ny))*ph)
        return Field(data=pf, psize=self.psize, l=self.l)

    def propagate_ae_d(self, d, n, shape=None):
        """Propagate the field a distance z using angular spectrum method.
        The observation plane is shifted by x,y units
        """

        if shape is None:
            nx, ny=self.data.shape
        else:
            nx, ny=shape

        (x, y, z)=d

        dx, dy=self.res
        sx=dx*nx
        sy=dy*ny
        # Phase function calculation
        X, Y=indices((nx, ny))
        X=X-nx/2
        Y=Y-ny/2

        X=ifftshift(X)
        Y=ifftshift(Y)

        # Check if the fftshift was OK
        assert(X[0, 0]==0)
        assert(Y[0, 0]==0)
        zeta=X/sx
        eta=Y/sy
        sq=(1./self.l)**2 -(zeta)**2-(eta)**2
        # TODO: Check what does this condition mean
        ph=exp(2.*pi*1.j*(x*zeta+y*eta+z*sqrt(sq)))
        ph=where(sq<0, 0., ph)
        pf=ifft2(fft2(self.data, (nx, ny))*ph)
        return Field(data=pf, psize=self.psize, l=self.l)

    def propagate_rsi_gpu(self, z, n, dfield=None, offset=(0., 0.), gpu=0):
        print("propagate_rsi_gpu: This method is not working fine needs to be checked")
        import pyopencl as cl
        pl=cl.get_platforms()[0]
        dev=pl.get_devices(device_type=cl.device_type.GPU)
        print(dev)

        # Inicio del programa

        prg_src = """
        __kernel void rsk(
                          __global float *x,
                          __global float *y,
                          __global float *ur,
                          __global float *ui,
                          unsigned h,
                          unsigned w,
                          float Z,
                          float K,
                          float dx,
                          float dy,
                          __global float *resr,
                          __global float *resi)
        {
            int nWidth = get_global_size(0);
            int nHeight = get_global_size(1);

            int ox=get_global_id(0); // Toma los indices en X
            int oy=get_global_id(1); // Toma los indices en Y
            int oid= oy*nWidth+ox;

            float X=(float)(ox-nWidth/2)*dx;
            float Y=(float)(oy-nHeight/2)*dy;
            float pi=3.14159265358979323846264338327950288;
            resr[oid]=0;
            resi[oid]=0;
            float Z2=Z/2.;
            for(unsigned i=0;i<h*w;i++)
            {
                float wi=2.;
                float wj=2.;

                int ci=i%w;
                int cj=i/w;

                //Hay que verificar si los i y los j estan correctos y si los
                //coheficientes estan al derecho
                if ((ci==0)||(ci==(w-1)))wi=1;
                else if(ci%2!=0) wi=4;

                if ((cj==0)||(cj==(h-1)))wj=1;
                else if(cj%2!=0) wj=4;

                float ww=wi*wj/9.;



                float Xx=X-x[i];
                float Yy=Y-y[i];
                float R2=pow(Xx,2.f)+pow(Yy,2.f)+pow(Z,2.f);
                float R= pow(R2,0.5f);
                float R3=R2*R;
                float sinkr=sin(R*K);
                float coskr=cos(R*K);
                float piR2=pi*R2;
                float piR3=pi*R3;

                //real 1/2*k*z*sin(R*k)/(piR2) + 1/2*z*cos(R*k)/(piR3)
                //imag -1/2*k*z*cos(R*k)/(piR2) + 1/2*z*sin(R*k)/(piR3)
                float re= K*Z2*sinkr/(piR2) + Z2*coskr/(piR3);
                float im= -K*Z2*coskr/(piR2) + Z2*sinkr/(piR3);
                resr[oid] +=  ww*(re*ur[i]-im*ui[i])*dx*dy;
                resi[oid] +=  ww*(re*ui[i]+im*ur[i])*dx*dy;
            }
        }
        """

        ctx0 = cl.Context(devices=[dev[gpu]])
        queue0 = cl.CommandQueue(ctx0)
        prg0 = cl.Program(ctx0, prg_src).build()
        mf = cl.mem_flags

        if dfield is None:
            ushape=self.data.shape
            dx, dy=self.res
        else:
            ushape=dfield.data.shape
            dx, dy=dfield.res

        ox=offset[0]
        oy=offset[1]
        # Get the coordinates from the object field and translate them
        # according to the offset
        Xo, Yo=self.field_sample_coord
        Xo=(Xo-ox).astype(float32)
        Yo=(Yo-oy).astype(float32)

        # Get the info from the object field
        Uor=self.data.real.astype(float32)
        Uoi=self.data.imag.astype(float32)

        # Create the object field buffers for the GPU
        Xob0 = cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Xo)
        Yob0 = cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Yo)
        Uobr0= cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Uor)
        Uobi0= cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=Uoi)

        # Create the buffers for the image field
        ifr=empty(ushape, dtype=float32)
        ifi=empty(ushape, dtype=float32)

        # Create the image buffers for the GPU

        ifbr=cl.Buffer(ctx0, mf.WRITE_ONLY, ifr.nbytes)
        ifbi=cl.Buffer(ctx0, mf.WRITE_ONLY, ifi.nbytes)

        prg0.rsk(queue0, ushape,
                 Xob0,
                 Yob0,
                 Uobr0,
                 Uobi0,
                 uint32(Xo.shape[0]),
                 uint32(Xo.shape[1]),
                 float32(z),
                 float32(2.*pi/self.l),
                 float32(dx),
                 float32(dy),
                 ifbr,
                 ifbi,)
        #        local_size=(16,16)
        #        )
        # print "nota: Al usar el metodo de integración de simpson, hay que tener"
        #       "un buen"
        # print "muestreo en el campo objeto. Si se tiene un solo pixel en 1 y lo demas"
        #       " en 0 "
        # print "Puede haber problemas, por que el muestreo no es adecuado."
        # print "Para generar frentes de onda creados por un filtro espacial, tratar de"
        #        "hallar "
        # print "la solución teorica"

        cl.enqueue_read_buffer(queue0, ifbr, ifr).wait()
        cl.enqueue_read_buffer(queue0, ifbi, ifi).wait()
        U=ifr+1.j*ifi
        if dfield is None:
            dfield=Field(data=U, psize=(dx, dy), l=self.l)
        dfield.data=U
        return dfield

    def propagate_rsi1(self, z, n, dfield=None):
        """Propagate the field a distance z, using the Rayleigh Sommerfeld integral.
        If dest is not given, the returned field has the same size and discretization
        of the propagated filed, else the Field class instance (dfield) is filled
        with propagated values, keeping its size and discretization. The value of
        the propagated field is returned.
        """
        from cfield import _propagate_rsi
        return _propagate_rsi(self, z, n, dfield)

        # """    from cfield import _ursi
        # #cpus=detectCPUs()

        # #job_server=pp.Server()

        # #TODO: Need to find a way to not recalculate the _RS_KERNEL
        # if dfield==None:
        #     X,Y=self.field_sample_coord
        # else:
        #     X,Y=dfield.field_sample_coord
        # U=empty((X.shape[0]*X.shape[1],), complex)
        # xf=X.flatten()
        # yf=Y.flatten()

        # cxy=column_stack((xf[:,newaxis],yf[:,newaxis]))
        # for i,(xi,yi) in enumerate(cxy):
        #    #U[i]=self._ursi(xi,yi, z, n=n)
        #     U[i]=_ursi(self,xi,yi, z, n=n)

        # if dfield==None:
        #     return Field(data=U.reshape(X.shape),psize=self.psize,l=self.l)
        # else:
        #     dfield.data=U.reshape(X.shape)
        #    return dfield
        # """

    def propagate_rsi(self, z, n, dfield=None):
        """Propagate the field a distance z, using the Rayleigh Sommerfeld integral.
        If dest is not given, the returned field has the same size and discretization
        of the propagated filed, else the Field class instance (dfield) is filled
        with propagated values, keeping its size and discretization. The value of
        the propagated field is returned.
        """

        # from cfield import _ursi
        # cpus=detectCPUs()

        # job_server=pp.Server()

        # TODO: Need to find a way to not recalculate the _RS_KERNEL
        if dfield is None:
            X, Y=self.field_sample_coord
        else:
            X, Y=dfield.field_sample_coord
        U=empty((X.shape[0]*X.shape[1],), complex)
        xf=X.flatten()
        yf=Y.flatten()

        cxy=column_stack((xf[:, newaxis], yf[:, newaxis]))
        for i, (xi, yi) in enumerate(cxy):
            U[i]=self._ursi(xi, yi, z, n=n)
            # U[i]=_ursi(self,xi,yi, z, n=n)

        if dfield is None:
            return Field(data=U.reshape(X.shape), psize=self.psize, l=self.l)
        else:
            dfield.data=U.reshape(X.shape)
            return dfield

    def check_z(self, n=1., shape=None):
        """Given the field characteristics, return the tipical z's that should
        be used in the different propagation algorithms.
        """
        # TODO: Deduce equations
        dx, dy=self.res
        if shape is None:
            nx, ny=self.data.shape
        else:
            nx, ny=shape
        sx, sy=dx*nx, dy*ny

        if dx!=dy:
            warn("The field has different X and Y resolutions. Using average")

        if nx!=ny:
            warn("The field has different X and Y samplings. Using average")

        if sx!=sy:
            warn("The field has different X and Y dimensions. Using average")

        dx=(dx+dy)/2.
        sx=(sx+sy)/2.
        nx=(nx+ny)/2.

        U0=1./(2.*dx)
        dU=1./(sx)

        # Maximum distance  that should be used for angular espectral
        # propagation
        ae_max= 1/(2.*(sqrt((n/self.l)**2 -(U0-dU)**2)-
                   sqrt((n/self.l)**2 -(U0)**2)))

        # Minimum distance to be used by Rayleigh Sommerfeld
        rs_min=(n/self.l)*sqrt(((self.l/(2*n))**2- sx**
                                2 -(sx-dx)**2)**2 -4*sx**2*(sx-dx)**2)

        # Minimum distance to be used by the fresnel propagation method (FFT)

        fr_min= dx**2*n*nx/self.l

        return {"ae_max": float(ae_max), "rs_min": float(
            rs_min), "fr_min": float(fr_min)}

    # method(This, Trait(None,Float,Array(shape=(3,))))

    def propagate_rs(self, z, n=1, shape=None):
        """Method that calculates the free space propagation of an optical field.

        This method calculates the optical propagation of a field, using
        the Rayleigh Sommerfeld propagation kernel.
        Depending on the propagation distance, this method  selects between
        the angular spectrum method, and the Rayleigh Sommerfeld convolution method

        **ARGUMENTS:**

            ===== =======================================================
            z     Propagation distance.
            n     Media refraction index (Not implemented yet assume n=1)
            shape Tuple indicating the shape (number of pixels) of the
                  returned field. If shape = None, the field size is
                  preserved.
            ===== =======================================================

            The other parameters are kept.

        **RETURN VALUE:**

            Field instance containing the propagated optical field
        """
        if shape is None:
            nx, ny=self.shape
        else:
            nx, ny=shape

        # Get the approximations validity limits
        d=self.check_z(n=n, shape=shape)
        dae=d["ae_max"]
        drs=d["rs_min"]
        print(dae, drs, z)
        dx, _dy=self.res
        if z< dae:
            print("ae")
            return self.propagate_ae(z, n, shape)
        elif z>drs:
            print("rs")
            return self.propagate_rsc(z, n, shape)
        else:
            print("Error: It is not possible to calculate using this distance")
            # Corregir esto, haciendo un remuestreo del campo, y usando espectro angular
            # TODO: Fix this when dx and dy are different
            a=1./(2*z)+sqrt((n/self.l)**2-(1/(2*dx))**2)
            nxmin=int(2/(1-2*dx*sqrt((n/self.l)**2-a**2)))
            nymin=nxmin
            print(nx, nxmin)
            if (nxmin-nx) %2!=0:
                nxmin=nxmin+1
            if (nymin-ny) %2!=0:
                nymin=nymin+1

            f=self.resize((nxmin, nxmin))
            pf=f.propagate_ae(z, n)
            # TODO: Use the next 2 power so the FFT is optimized
            warn("Calculating propagation using a resized field")

            return pf.resize((nx, ny))

    def propagate_fraunhofer(self, z):
        """
        Propagate the field using the Fraunhofer approximation

        **ARGUMENT:**

                =  ====================
                z  Propagation distance
                =  ====================
        """

        sx, sy=self.size
        l=self.l
        if z>0:

            data=fftshift(fft2(ifftshift(self.data)))
            # data=fftshift(fft2(self.data))
            dxp=l*z/sx
            dyp=l*z/sy

        if z<0:
            # ~ data=fftshift(ifft2(ifftshift(self.data)))
            data=ifftshift(ifft2(fftshift(self.data)))
            dxp=-l*z/sx
            dyp=-l*z/sy

        print("warning: The phase factor is not present. Must be corrected")
        return Field(data=data, psize=(dxp, dyp), l=self.l)

    def propagate_fresnel(self, z):
        '''
        Calculate the fresnel transform, using the FFT algorithm.

        **ARGUMENT:**

                =  ====================
                z  Propagation distance
                =  ====================


        '''
        dx, dy=self.res

        sx, sy=self.size
        l=self.l

        dxp=abs(l*z/sx)
        dyp=abs(l*z/sy)

        X, Y=meshgrid(self.xsamples, self.ysamples)

        Xp, Yp=dxp*X/dx, dyp*Y/dy

        ph=exp(1.j*pi*(X**2+Y**2)/(z*l))

        ph1=exp((2.j*pi/l)*z)/(1.j*l*z)
        ph2=exp(1.j*pi*(Xp**2+Yp**2)/(z*l))

        if z>0:
            data=ph1*ph2*fftshift(fft2(ifftshift(self.data*ph)))
        elif z<0:
            data=ph1*ph2*ifftshift(ifft2(fftshift(self.data*ph)))

        return Field(data=data, psize=(dxp, dyp), l=self.l)

    def propagate(self, r, n=1., method="auto", fix=False):
        """Method that calculates the free space propagation of an optical field.


        **ARGUMENTS**

        ======= ========================================================
        r       Vector that goes from the source plane origin, to the
                destination plane origin. If r is not a vector but a
                floating point number, it assumes that the vector is
                (0,0,r).
        n       Media refraction index (Not implemented yet assume n=1)
        method  Propagation method:

                "ae"   Angular spectrum propagation method

                "rsc"  Convolution with the Rayleigh Sommerfeld Kernel

                "rsi"  Propagate the field using the Rayleigh Sommerfeld
                       integral

                "???"

                "auto" Automatically select the best propagation method
                       given the field conditions.


        fix     If set to true, the field is re-sampled or resized, so
                the sampling conditions for the selected propagation
                method are met.
        ======= ========================================================

        **RETURN VALUE:**

        Field instance containing the propagated optical field
        """

        try:
            x, y, z=r
        except ValueError:
            x, y, z= 0., 0., r

        dx, dy=self.res
        nx, ny=self.data.shape
        sx, _sy=dx*nx, dy*ny

        m=None

        if method=="auto":
            if z<=self.check_z(n)["ae_max"]:
                m="ae"
            else:
                m="crs"
        else:
            m=method

        # TODO: , Warning n not used in calculations must be fixed
        print("Warning n not used in calculations must be fixed")
        if m=="ae":
            tf=self
            if z>self.check_z(n)["ae_max"]:
                wm="The propagation distance (z=%f) is greater than the "\
                   "recommended (z_max=%f) for the angular espectrum "\
                   "propagation method" %(z, self.check_z(n)["ae_max"])
                warn(wm)
                # TODO: Deduce equations
                if fix:
                    # TODO: Fix this when dx and dy are different
                    a=1./(2*z)+sqrt((n/self.l)**2-(1/(2*dx))**2)
                    nxmin=int(2/(1-2*dx*sqrt((n/self.l)**2-a**2)))
                    # TODO: Use the next 2 power so the FFT is optimized
                    warn("Calculating propagation using a resized field")

                    # This is done because the number of rows and columns to add
                    # must be even.
                    # TODO: This is not OK, if the number of columns is even and rows
                    # is odd (or backwards), this will not work. FIXME
                    print(nxmin)
                    try:
                        tf=self.resize((nxmin, nxmin))

                    except ValueError:
                        tf=self.resize((nxmin+1, nxmin+1))
            retval=tf.propagate_ae(z, n)
        elif m=="rsc" or m=="rsi":
            if z<self.check_z(n)["rs_min"]:
                wm="The propagation distance (z=%f) is smaller than the "\
                   "recommended (z_min=%f) for the Rayleigh Sommerfeld integral "\
                   "propagation method" %(z, self.check_z(n)["ae_max"])
                warn(wm)
                # TODO: Deduce equations
                if fix:
                    # TODO: Fix this when sx and sy are different
                    k=2*pi/self.l
                    dxmax=pi/(k*sx)*sqrt(sx**2+z**2)

                    rx=1+int(sx/dxmax)
                    warn("Calculating using a resampled field")
                    tf=self.resample((rx, rx))
                else:
                    tf=self

            if m=="rsc":
                retval=tf.propagate_rsc(z, n)
            else:
                retval=tf.propagate_rsi(z, n)
        else:
            raise ValueError("No suitable propagation method found")
            retval=None

        return retval

    def resize(self, samples):
        """Returns a resized optical field. The resolution is not modified.

        * samples *
            Tuple (nx,ny) indicating the new number of samples of the
            optical field.

        Note: In all resizes, the origin of the optical field is
        preserved (the origin is always at the center)
        """

        nx, ny=self.data.shape
        nnx, nny=samples

        mnx=max((nx, nnx))
        mny=max((ny, nny))

        # The number of rows and columns to add must be even
        if (((mnx-nx) %2 !=0)or((mny-ny) %2 !=0)):
            raise ValueError(
                "The number of rows and columns to add  or remove must be even")

        ndata=zeros((mnx, mny), complex)
        ix=(mnx-nx)/2
        iy=(mny-ny)/2

        ndata[ix:ix+nx, iy:iy+ny]=self.data

        if nx>nnx:
            ndata=ndata[(nx-nnx)/2:(nx-nnx)/2+nnx, :]

        if ny>nny:
            ndata=ndata[:, (ny-nny)/2:(ny-nny)/2+nny]

        retval=Field(data=ndata, psize=self.psize, l=self.l)
        return retval

    def resample(self, samples):
        """Returns a resampled optical field The size is not modified.

        * res *
            Tuple (rx,ry) indicating the new width and height of the optical
            field in pixels.

        Note: In all resamples, the origin of the optical field is
        preserved (the origin is always at the center)
        """
        dx, dy=self.res
        nx, ny=self.data.shape
        sx=dx*nx
        sy=dy*ny
        Zr=self.data.real
        Zi=self.data.imag
        nnx, nny=samples
        xvals=linspace(0, nx, nnx, endpoint=True)
        yvals=linspace(0, ny, nny, endpoint=True)
        X, Y=meshgrid(xvals, yvals)
        coords = array([X, Y])
        outgr = map_coordinates(Zr, coords, order=1)
        outgi = map_coordinates(Zi, coords, order=1)

        data=outgr+1.j*outgi

        retval=Field(data=data, psize=(sx/nnx, sy/nny), l=self.l)
        return retval

    def tilt(self, r=(0., 0., 0.)):
        """Rotate around the origin the field observation plane

        Idea Taken from:
            Free-space beam propagation between arbitrarily oriented planes based
            on full diffraction theory: a fast Fourier transform approach.
            The interpolation algorithm suggested in the paper, does not work. This
            routine is based on a novel interpolation algorithm.

        r
            Tuple (rx,ry,rz) where rx is the rotation around the x axis, ry is
            the rotation around the y axis and rz is the rotation around the z axis
            that must be issued to the object plane, to obtain the image plane.

        The rotations are applied first to the z axis, and then to the y axis and
        finally to the x axis.
        (need to verify this to check if this is consistent with wxRayTrace).

        """

        # TODO: Put the complete reference

        dx, dy=self.res
        nx, ny=self.data.shape
        sx=dx*nx
        sy=dy*ny

        RM=rot_mat(r)

        # check if there is a need to resample the wavefront
        z_o=0
        e_o=0
        k_o=0
        sq=(1./self.l)**2 -(z_o)**2-(e_o)**2
        k_o=sqrt(sq)
        z_r, e_r, k_r=dot(RM, [z_o, e_o, k_o])
        kr=z_r*sx
        lr=e_r*sy

        fk=2.*abs(kr)/nx
        fl=2.*abs(lr)/ny
        fm=[fk, fl]
        fm.sort()
        nf=self
        if fm[-1]>1:
            fm=int(fm[-1])+1
            print("Resampling field ", fm)
            nf=self.resample((nx*fm, ny*fm))
            nx, ny=nf.data.shape
            dx, dy=nf.res
            sx=dx*nx
            sy=dy*ny

        ix, iy=indices((nx, ny)).astype(float)
        nx2=int(nx/2)
        ny2=int(ny/2)
        # Calculate the displacement matrix for the original planewave spectrum
        c1=exp(2.j*pi/nx*(nx2*(ix-nx2)))*exp(2.j*pi/ny*(ny2*(iy-ny2)))

        ix=(ix-int(nx/2)).flatten()
        iy=(iy-int(ny/2)).flatten()

        # Calculate the K vector components in the original space.
        zeta_o = ix/sx
        eta_o = iy/sy

        sq=(1./nf.l)**2 -(zeta_o)**2-(eta_o)**2
        kappa_o=npsqrt(sq)

        # Rotate to obtain the K vectors in the rotated space

        # TODO: Find a better way to do this

        # Obtain the k vectors in the original space
        A=fftshift(fft2(nf.data))*c1
        Af=A.flatten()

        # Matrix to put the rotated coefficients vectors
        Ar=zeros_like(A)

        for i in range(nx*ny):
            if sq[i]>0:
                # Calculate the rotated k vector and the corresponding indexes
                z_r, e_r, k_r=dot(RM, [zeta_o[i], eta_o[i], kappa_o[i]])

                kr=z_r*sx+int(nx/2)
                lr=e_r*sy+int(ny/2)

                # Calculate the displacement matrix for the rotated planewave
                # spectrum
                c2=exp(-2.j*pi/nx*(nx2*(kr-nx2)))* exp(-2.j*pi/ny*(ny2*(lr-ny2)))
                Ar=Ar+Af[i]*c2*s2d(kr, lr, nx, ny)

        data=ifft2(ifftshift(Ar/(nx*ny)))
        return Field(data=data, psize=nf.psize, l=nf.l)

    def rayrep(self, nx, ny, eps=1e-2):
        """
        Method to calculate the ray representation of the wavefront.

        nx  Number of samples to use in the x direction
        ny  Number of samples to use in the Y direction
        eps Allowed maximum value for lap a/ a (eikonal condition)

        Note: This ray representation has only information about the phase.
            the intensity, is not yet represented in the ray representation.
            this needs to be solved.
        """
        # Get field amplitude and unwrapped phase
        a=abs(self.data)
        ph=self.phase

        # Get field mask
        fmask=self.mask

        # TODO: Check if it is better to calculate the gradient using a polynomial fit
        # might be more accurate

        # Calculate using finite differences: f' (x)~(f(x+h)-f(x-h))/(2h)
        #                                   f''(x)~(f(x+h)+f(x-h)-2f(x))/(h^2)

        # Create displaced matrices and repeat the data at the borders
        fx1=zeros_like(ph)
        fx2=zeros_like(ph)
        fy1=zeros_like(ph)
        fy2=zeros_like(ph)

        fx1[1: -1, :]=ph[:-2, :]
        fx1[0, :]=fx1[1, :]
        fx1[-1, :]=fx1[-2, :]

        fx2[1: -1, :]=ph[2:, :]
        fx2[0, :]=fx2[1, :]
        fx2[-1, :]=fx2[-2, :]

        fy1[:, 1: -1]=ph[:, :-2]
        fy1[:, 0]=fy1[:, 1]
        fy1[:, -1]=fy1[:, -2]

        fy2[:, 1: -1]=ph[:, 2:]
        fy2[:, 0]=fy2[:, 1]
        fy2[:, -1]=fy2[:, -2]

        rx, ry=self.res

        # Calculate the partial derivatives respect X and Y and normalize them
        dx=((fx2-fx1)/(2*rx))/(2*pi/self.l)
        dy=((fy2-fy1)/(2*ry))/(2*pi/self.l)

        # Check for the eikonal condition

        # Calculate the second partial derivatives of a
        # Create displaced matrices and mark invalid data
        mask=zeros_like(a)
        mask[0, :]=True
        mask[-1, :]=True
        mask[:, 0]=True
        mask[:, -1]=True

        ax1=maarray(zeros_like(a), mask=mask)
        ax2=maarray(zeros_like(a), mask=mask)
        ay1=maarray(zeros_like(a), mask=mask)
        ay2=maarray(zeros_like(a), mask=mask)

        ax1[1: -1, :]=a[:-2, :]
        ax2[1: -1, :]=a[2:, :]
        ay1[:, 1: -1]=a[:, :-2]
        ay2[:, 1: -1]=a[:, 2:]

        a2x=(ax1+ax2-2*a)/(rx**2)
        a2y=(ay1+ay2-2*a)/(ry**2)

        # TODO: Find a better criteria to calculate rz
        rz=(rx+ry)/2.
        a0=self.propagate(-rz/10000., n=1., method="ae").abs()
        a1=a
        a2=self.propagate(rz/10000. , n=1., method="ae").abs()
        a2z=(a0+a2-2.*a1)/(rz**2)

        lapa=(a2x+a2y+a2z)/a

        # TODO: Idea. We need to study a better way to justify the criteria.
        #      Also we need to study how big the laplacian can be and what
        #      error it creates.

        if abs(lapa).max()>eps:
            # TODO: establish a better warning system, and apply it to all warnings
            # in the program
            print("Warning: the conversion from wavefront to rays, does not fulfil")
            print("the eikonal equation max(lap(a)/a)", abs(lapa).max())

        # TODO: Check if we need a condition for the continuity of the phase, or
        # if this condition is enough

        # Resample to the number of rays
        onx, ony=self.shape
        coord=mgrid[1:onx-2:1.j*nx, 1:ony-2:1.j*ny]

        coord=tuple(rint(coord).astype("i"))

        dx=dx[coord]
        dy=dy[coord]

        # Calculate partial derivative respect to Z
        dz2= 1.-dx**2-dy**2

        # Check that there are no inconsistencies in the data
        # TODO: Fix the assert
        # assert alltrue(dz2>=0), "There is aninconsistency in the gradient calculation"
        dz=sqrt(dz2)

        X, Y=self.field_sample_coord
        X=X[coord]
        Y=Y[coord]
        fmask=fmask[coord]
        fmask=fmask.flatten()

        ph=ph[coord]
        ph=ph.flatten()

        dx=dx.flatten()
        dy=dy.flatten()
        dz=dz.flatten()
        X=X.flatten()
        Y=Y.flatten()

        print("Warning: The assignment of the optical path assumes n=1 \n"
              "This must be corrected.")

        rl=[]
        dirs=zip(dx, dy, dz, X, Y,  fmask, ph)
        for dx, dy, dz, x, y,  mask,  path in dirs:
            if dz<0:
                print(dz)
            if mask is False:
                rl.append(Ray(origin=(x, y, 0), direction=(dx, dy, dz),
                          wavelength=self.l, pop=self.l*path/(2.*pi)))
        return rl  # ,lapa,a2x,a2y,a2z
