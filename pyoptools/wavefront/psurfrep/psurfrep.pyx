# cython: profile=True

import numpy as np
cimport numpy as np


#from numpy import dot,arctan2, polyval, pi, indices, sort, cumsum, where
from numpy.fft import fft2, fftshift, ifftshift
import numpy.ma as ma
from pyoptools.misc import *
from time import time
from pyoptools.misc.Poly2D import ord2i,poly2d
cdef extern from "math.h":
    double sqrt(double) nogil
    double atan2(double,double) nogil


from pyoptools.wavefront.field import Field

#@cython.boundscheck(False) # turn of bounds-checking for entire function
#@cython.wraparound(False)

class PSurf:
    """Class used to describe the wave transfer function of an optical
    surface using a polynomial representation algorithm.
    """
    
    def __init__(self, s, ni=1,nr=1,ilimit=0, slimit=1.3, l=.633e-3 , step=0.05, order=10, rsamples=(500,500),zb=None):
        """
        
        **ARGUMENTS**

        s  -- Optical surface to model
        ni,nr  -- Refraction index from the incident and refracted sides
        ilimit -- Inferior limit for incidence angle of the plane wave in radians
        slimit -- Superior limit for incidence angle of the plane wave in radians
        l      -- Wavelength that will be used in the simulation. We need to find a solution for any wavelength
        step   -- Step to be used to generate the interpolation data
        order  -- Order of the Taylor interpolation used
        rsamples -- Tuple containing the number of ray samples to be used en 
                    each direction
        zb    -- Z position of the plane where the measurementas are made. The origin 
                 is the vertex of the surface. If None, an estimate position is taken.
        
        Notes: 
               - The pupil is normalized to the radius of the lens
               - This method asumes a circular shaped pupil
               - The phase is returned as an optical path
               - The intensity is normalized to 1
            
        """
        self.pf,self.pi,self.zm=s.pw_cohef(ni,nr,ilimit, slimit, step, order, rsamples,zb)
        
        #get the number of polynomial coheficients
        self.nc=ord2i(order)
        
        self.l=l
    def wf_polys(self,np.ndarray[np.double_t, ndim=1] k):
        '''Return the polynomials that describe the wavefront
        
        (fase_poly, i_poly)
        '''
        cdef double r, iang 
        cdef int i
        r=sqrt(k[0]**2+k[1]**2)
        iang =atan2(r,k[2])
        
        cdef np.ndarray[np.double_t, ndim=1] cf=np.zeros((self.nc,), dtype=np.double)
        cdef np.ndarray[np.double_t, ndim=1] ci=np.zeros((self.nc,), dtype=np.double)
        for i in range(self.nc):
            pof=self.pf[i]
            cf[i]=np.polyval(pof,iang)
            poi=self.pi[i]
            ci[i]=np.polyval(poi,iang)
        
        df=poly2d(cf)
        di=poly2d(ci)
        return df,di
        
    
    def pw_evaluate(self,np.ndarray[np.double_t, ndim=1] k, samples=(512,512), gpu=True):
        """Plane wave evaluate
        
        **Arguments:**
        
        k -- Propagation direction of the incident plane wave 
        
        """
        cdef double r, iang 
        cdef int i
        r=sqrt(k[0]**2+k[1]**2)
        iang =atan2(r,k[2])
        
        cdef np.ndarray[np.double_t, ndim=1] cf=np.zeros((self.nc,), dtype=np.double)
        ci=np.zeros((self.nc,))
        for i in range(self.nc):
            pof=self.pf[i]
            cf[i]=np.polyval(pof,iang)
            poi=self.pi[i]
            ci[i]=np.polyval(poi,iang)
        
        df=poly2d(cf)
        di=poly2d(ci)
        
        xxl,dx=np.linspace(-1., 1., samples[0],retstep=True)
        yyl,dy=np.linspace(-1., 1., samples[1],retstep=True)
        xx,yy=np.meshgrid(xxl,yyl)
        
        #Create the circular mask
        
        # There is a problem. Some of the pixels of the border of the aperture get an error
        # Too big. I don't think this will be an issue, because it is just in the border.
        
        rm=np.where(xx**2+yy**2>1,True,False)
        
        rotang=-atan2(k[1],k[0])
        
        if not gpu:
            fr=df.mevalr(xx,yy,rotang)
            ir=di.mevalr(xx,yy,rotang)
        else:
            fr=df.gpu_evalr(xxl,yyl,rotang)
            ir=di.gpu_evalr(xxl,yyl,rotang)
        
        f=ma.masked_array(fr, mask=rm)
        
        
        inte=ma.masked_array(ir, mask=rm)
        fdata=inte*np.exp(1.j*2*np.pi/self.l*f)
        

        #fdata=fdata*fdata[samples[0]/2,samples[1]/2]/abs(fdata[samples[0]/2,samples[1]/2])
        print "**-->",np.angle(fdata[samples[0]/2,samples[1]/2])
        return Field(fdata,psize=(dx,dy), l=self.l)

        
        
    #@cython.boundscheck(False) # turn of bounds-checking for entire function
    #@cython.wraparound(False)        
    def propagate(self,wf , pt=0.001, samples=(512,512)):
        """
        Propagate a Wavefront through the modelled surface.
        
        wf -- wave front to propagate
        pt -- Percentage of the power spectrum that must be truncated. 
              If no truncation is intended, set to 1
        samples -- Samples of the resulting wave
        """
        
        #Calculate the plane wave spectrum
        data=wf.data
        cdef double dx,dy,sx,sy
        cdef int nx,ny
        nx,ny=wf.shape
        dx,dy=wf.res
        
        sx=dx*nx; sy=dy*ny
        print wf.shape,wf.res
        # Phase function calculation
        #cdef np.ndarray[np.complex128_t, ndim=3] a = \
        #np.zeros((3,3,3), dtype=np.complex128)
        
        cdef np.ndarray[np.double_t, ndim=2] X #= np.empty([nx, ny], dtype=np.double) 
        cdef np.ndarray[np.double_t, ndim=2] Y #= np.zeros([nx, ny], dtype=np.double) 
        
        
        X,Y=np.indices((nx,ny), dtype=np.double)
        X=X-nx/2; Y=Y-ny/2
        #cdef int i,j
        #for i in range(nx):
        #    for j in range(ny):
        #        X[i,j]=i-nx/2
        #        Y[i,j]=j-ny/2
        
        
        

        X=np.fft.ifftshift(X)
        Y=np.fft.ifftshift(Y)
        
        # Check if the fftshift was OK
        assert(X[0,0]==0)
        assert(Y[0,0]==0)
        
        #Calculate the k vector corresponding to each sample of the fft2
        kx=(X/sx)
        ky=(Y/sy)
        print sx,sy
        kz=np.sqrt((1./wf.l)**2 -kx**2-ky**2)
        
        #Calculate the amplitude and phase of each plane wave
        a=np.fft.fft2(data)
        
        #Eliminate the least important components
        
        sd=np.sort(abs(a.flatten()))
        #Obtain something similar to the cumulative power, and normalize it
        cs=np.cumsum(sd)
        cs=cs/cs[-1]
        #Get the value where the truncation must occur
        i=where(cs>pt)[0][0]-1
        #~ print cs
        #~ print i

        #eliminate the unwanted components
        #ks=where(abs(a)>sd[i],a,0)
        ki,kj=np.where(abs(a)>sd[i])
        #####
        k=np.array((kx[ki[0],kj[0]],ky[ki[0],kj[0]],kz[ki[0],kj[0]]))
        # got_cl is importen from misc. Misc imports it from Poly2D
        r=a[ki[0],kj[0]]*self.pw_evaluate(k,samples,got_cl)
        print "*****",len(ki),k/k[2],abs(a[ki[0],kj[0]]),k[0],k[1],k[2]
        
        
        for i in range(len(ki)-1):
            print "*",i
            k=np.array((kx[ki[i+1],kj[i+1]],ky[ki[i+1],kj[i+1]],kz[ki[i+1],kj[i+1]]))
            print k/k[2],abs(a[ki[i+1],kj[i+1]])
            r+=a[ki[i+1],kj[i+1]]*self.pw_evaluate(k,samples,got_cl)
            
            ####TODO: Important------------------------------- ###### Need to take into account the phase............

            ### Need to add the polynomials before evaluating the resulting field. Have to check if this is possible, or if 
            ### it is better to add the evaluated fields
            ### Or if it is possible to do an interpolation in the fourier plane
            
            #####resf=resf+resf*exp(1.j*)
            #~ if i%1 ==0:
                #~ print "**", time()-ti
                #~ print "***",i, len(ki)
                #~ ti=time()
        return r
