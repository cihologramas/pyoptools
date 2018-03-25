#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#cython: profile=True

cdef extern from "math.h":
    double pow(double,double)
    double sin(double)
    double cos(double)

cdef bint got_cl=True


        
try: 
    import pyopencl as cl
except ImportError: 
    got_cl=False

import numpy as np
cimport numpy as np
cimport cython
np.import_array()

#Import C stdlib
from libc.stdlib cimport *
from pyoptools.misc.cmisc.cmisc cimport zero_vec, zero_mat, empty_mat

cimport cython

## Module initialization
if got_cl:
    init_cl()  


cdef init_cl():
        """Routine to initialize the OpenGL kernels"""
        global ctx0
        global queue0
        global prg0
        global prg1
        global prg2
        
        #Look for  the first available GPU
        #for GPU, create the contexts and the programs
        #Its better not tu use CPU for calculatibg the polynomial. 
        #It is as slow as an optimized cython program.
        #TODO: Needs to check if this is true
        pls=cl.get_platforms()
        foundgpu=False
        
        for pl in pls:
            devs=pl.get_devices(device_type=cl.device_type.ALL)
            for dev in devs:
                if dev.type==cl.device_type.GPU:
                    foundgpu=True
                    break
                
            if foundgpu: break
        
        ## Inicio del programa
	##Need to see how can cl_amd_fp64 can be automatically changed
	##to cl_fp64
        prg_src = """
        #pragma OPENCL EXTENSION cl_khr_fp64 : enable
        __kernel void rsk(
                            __global double *Xb,
                            __global double *Yb,
                            __global double *PXb,
                            __global double *PYb,
                            __global double *COb,
                            unsigned np,
                            unsigned ord,
                            __global double *RESb)
        {
            int nWidth = get_global_size(0);
            int nHeight = get_global_size(1);
            
            int ox=get_global_id(0); // Toma los indices en X
            int oy=get_global_id(1); // Toma los indices en Y
            int oid= oy*nWidth+ox;
           
            double Xp;
            double Yp;
            double X=Xb[ox];
            double Y=Yb[oy];
            double RES;
            RES=0;
            
            for(unsigned i=0;i<np;i++)
            {
                // Because of memory limitations, the powers can not be precomputed
                Xp=1;
                for(unsigned j=0;j<PXb[i];j++)
                    Xp=Xp*X;
                
                Yp=1;
                for(unsigned j=0;j<PYb[i];j++)
                    Yp=Yp*Y;
                    
                RES +=  COb[i]*Xp*Yp;
            }
            RESb[oid]=RES;
        }
        """
    
        ctx0 = cl.Context(devices=[dev])
        queue0 = cl.CommandQueue(ctx0)
        prg0 = cl.Program(ctx0,prg_src ).build()
         

        prg_src = """
        #pragma OPENCL EXTENSION cl_khr_fp64 : enable
        __kernel void rsk(
                          __global double *PXb,
                          __global double *PYb,
                          __global double *COb,
                          unsigned np,
                          unsigned ord,
                          __global double *RESb)
        {
            int nWidth = get_global_size(0);
            int nHeight = get_global_size(1);
            
            int ox=get_global_id(0); // Toma los indices en X
            int oy=get_global_id(1); // Toma los indices en Y
            int oid= oy*nWidth+ox;
            double Xp;
            double Yp;
            double X=2.*ox/(nWidth-1)-1.;
            double Y=2.*oy/(nHeight-1)-1.;
            double RES;
            RES=0;
            
            for(unsigned i=0;i<np;i++)
            {
            // Because of memory limitations, the powers can not be precomputed
                Xp=1;
                for(unsigned j=0;j<PXb[i];j++)
                    Xp=Xp*X;
                Yp=1;
                for(unsigned j=0;j<PYb[i];j++)
                    Yp=Yp*Y;
                RES +=  COb[i]*Xp*Yp;

            }
            RESb[oid]=RES;
        }
        """
        prg1 = cl.Program(ctx0,prg_src ).build()
        
        prg_src = """
        #pragma OPENCL EXTENSION cl_khr_fp64 : enable
        __kernel void rsk(
                          __global double *Xb,
                          __global double *Yb,
                          __global double *PXb,
                          __global double *PYb,
                          __global double *COb,
                          double sinr,
                          double cosr,
                          unsigned np,
                          unsigned ord,
                          __global double *RESb)
        {
            int nWidth = get_global_size(0);
            int nHeight = get_global_size(1);
            
            int ox=get_global_id(0); // Toma los indices en X
            int oy=get_global_id(1); // Toma los indices en Y
            int oid= oy*nWidth+ox;
            //rx=  x[i,j]*cosr-y[i,j]*sinr
            //ry= x[i,j]*sinr+y[i,j]*cosr
            double Xp;
            double Yp;
            double X=Xb[ox]*cosr-Yb[oy]*sinr;
            double Y=Xb[ox]*sinr+Yb[oy]*cosr;
            double RES;
            RES=0;
            
            for(unsigned i=0;i<np;i++)
            {
                // Because of memory limitations, the powers can not be precomputed
                Xp=1;
                for(unsigned j=0;j<PXb[i];j++)
                    Xp=Xp*X;
                Yp=1;
                for(unsigned j=0;j<PYb[i];j++)
                    Yp=Yp*Y;
                    
                RES +=  COb[i]*Xp*Yp;
            }
            RESb[oid]=RES;
        }
        """
        
        prg2 = cl.Program(ctx0,prg_src ).build()    


cdef class poly2d:  
    '''Class to define a 2D polynomial
        
        .. math::
            z=c0+
            c1*x+c2*y+
            c3*x^2+c4*x*y+c5*y^2+
            c6*x^3+c7*x^2*y+c8*x*y^2+c9*y^3+...
    '''
    def __cinit__(self, coh):
        """
        """
        
        cohef=np.array(coh, dtype=np.float64)
        if not cohef.flags['C_CONTIGUOUS']:
            # Array is not contiguous, need to make contiguous copy
            c = cohef.copy('C')
        else:
            c = cohef
        
        
        ##### Save the coheficients so can be used in a python way or a C fast way
        self.cohef=c    
        self.cohef_c= <np.float64_t*>np.PyArray_DATA(self.cohef)        
        #####
        
        self.clen=np.uint32(len(self.cohef))
        #~ if px==None or py==None:
        
        #### Save the powers so thy can be used in a python way or a C fast way
        # Remember if they are not saved as a class attribute, they will get garbage
        # collected
        
        self.px,self.py=i2pxpy(range(len(self.cohef)))
        
        self.px=np.array(self.px, dtype=np.float64)
        self.py=np.array(self.py, dtype=np.float64)
        
        if not self.px.flags['C_CONTIGUOUS']:
            self.px=self.px.copy('C')
        
        if not self.py.flags['C_CONTIGUOUS']:
            self.py=self.py.copy('C')
        
        self.px_c=<np.float64_t*>np.PyArray_DATA(self.px)
        self.py_c=<np.float64_t*>np.PyArray_DATA(self.py)
        #####
        
        
        
        self.px64=self.px.astype(np.float64)
        self.py64=self.py.astype(np.float64)
        mx=self.px.max()
        my=self.py.max()
        self.order=np.uint32((np.array((mx,my))).max())
        self.dx=None
        self.dy=None
        
    def __reduce__(self):
        args=(self.cohef,) #self.intensity,self.wavelength,self.n ,self.label,self.parent,self.pop,self.orig_surf)
        return(type(self),args)
        
        
    def __add__(self,other):
        if isinstance(self,poly2d): #Note, this is different than standard python, here mul==rmul. We need to see what is going on
            p=self
            a=other
        else:
            p=other
            a=self
        if isinstance(a,poly2d):
            o1=p.order
            o2=a.order
            if o1>o2:o=o1
            else: o=o2
            ncohef=zero_vec( pxpy2i(0,o)+1)
            #qncohef=np.zeros((pxpy2i(0,o)+1,), dtype=np.float64)
            ncohef[:len(a.cohef)]=a.cohef
            ncohef[:len(p.cohef)]=ncohef[:len(p.cohef)]+p.cohef
            return poly2d(ncohef)
        
        
        return NotImplemented 
    
    def __sub__(self,other):
        if isinstance(self,poly2d): #Note, this is different than standard python, here mul==rmul. We need to see what is going on
            p=self
            a=other
        else:
            p=other
            a=self
        if isinstance(a,poly2d):
            o1=p.order
            o2=a.order
            if o1>o2:o=o1
            else: o=o2
            ncohef=zero_vec( pxpy2i(0,o)+1)
            #ncohef=np.zeros((pxpy2i(0,o)+1,), dtype=np.float64)
            ncohef[:len(a.cohef)]=-a.cohef
            ncohef[:len(p.cohef)]=ncohef[:len(p.cohef)]+p.cohef
            return poly2d(ncohef)
        
        
        return NotImplemented 
    
    def __neg__(self):
        return poly2d(-self.cohef)
        
        
    @cython.boundscheck(False) 
    def __mul__(self,other):
		
        cdef int rxp,ryp,axp,ayp,pxp,pyp
        cdef int o1,o2
     
        if isinstance(self,poly2d): #Note, this is different than standard python, here mul==rmul. We need to see what is going on
            p=self
            a=other
        else:
            p=other
            a=self
        cdef np.ndarray[np.float64_t, ndim=1,  mode="c"] ncohef
        cdef unsigned int i,j,ir
        if isinstance(a,(float,int)):
            cohef=p.cohef*a
            return poly2d(cohef)
        
        elif isinstance(a,poly2d):
            o1=p.order
            o2=a.order
            ncohef=zero_vec( pxpy2i(0,o1+o2)+1)#np.zeros((pxpy2i(0,o1+o2)+1,), dtype=np.float64)
            
            
            for i in range((<poly2d>p).clen):#len(p.cohef)):
                for j in range((<poly2d>a).clen):#len(a.cohef)):
                    axp=<int>(<poly2d>a).px_c[j]
                    ayp=<int>(<poly2d>a).py_c[j] #i2pxpy(j)
                    pxp=<int>(<poly2d>p).px_c[i]
                    pyp=<int>(<poly2d>p).py_c[i]#i2pxpy(i)
                    
                    rxp=axp+pxp
                    ryp=ayp+pyp
                    
                    ir=pxpy2i(rxp,ryp)
                    ncohef[ir]=ncohef[ir]+(<poly2d>p).cohef_c[i]*(<poly2d>a).cohef_c[j]
            return poly2d(ncohef)
        
         
        return NotImplemented 
        
    def __str__(self):
        retval=""
        for i in range(len(self.cohef)):
            px,py=i2pxpy(i)
            c=self.cohef[i]
            
            if c!=0:
                if c!=1 or (c==1 and px==0 and py==0): retval=retval+str(c)
                if px!=0:
                    retval=retval+"x"
                    if px!=1: retval=retval+"^"+str(px)
                if py!=0:
                    retval=retval+"y"
                    if py!=1: retval=retval+"^"+str(py)
                
                retval=retval+"+"
        return retval[:-1]
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def dxdy(self):
        """
        Function that calculates the derivative with respect to ``X`` and ``Y`` of a
        given polynomial represented as a matrix.

        **Return Value**
        

            (dx,dy) 
                Tuple containing the Poly2D representing de derivative with 
                respect to ``X``, and the Poly2D containing the derivative with 
                respect to ``Y``
        
        Note: this method gets cached the first time it is called. If you modify
        the cohef attribute, the cache does not get updated
        """
        cdef np.ndarray[np.float64_t, ndim=1] Dx
        cdef np.ndarray[np.float64_t, ndim=1] Dy
        cdef int i,px,py,dxi,dyi
        
        if (self.dx==None) or (self.dy==None):
            Dx=zero_vec( self.clen) #np.zeros(self.cohef.shape,dtype=np.float64)
            Dy=zero_vec( self.clen) #Dy=np.zeros(self.cohef.shape,dtype=np.float64)
            
            #cdef np.ndarray[np.int_t, ndim=1] PX
            #cdef np.ndarray[np.int_t, ndim=1] PY
            
            
            #PX,PY=i2pxpy(range(self.clen))
            
            
            
            
            for i in range(1,self.clen):
                #px=PX[i]
                #py=PY[i]
                px=<int>self.px_c[i]
                py=<int>self.py_c[i]
                dxi=pxpy2i(px-1,py)
                dyi=pxpy2i(px,py-1)
        
                Dx[dxi]=Dx[dxi]+px*self.cohef_c[i]
                Dy[dyi]=Dy[dyi]+py*self.cohef_c[i]
                self.dx=poly2d(Dx)
                self.dy=poly2d(Dy)
        return self.dx,self.dy
        

    cpdef eval(self, x, y):
        '''
        Evaluate the polynomial at x,y
        
        This is the general method written in numpy. In this method x and y
        can have any shape (number, vectors, matrices, etc).
        Because of its generality, this method is not the fastest.
        
        Arguments:
            = =========================================================
            x Number, vector or matrix containing the values of x where
              the polynomial is goint to be evaluated.
            y Number, vector or matrix containing the values of y where
              the polynomial is goint to be evaluated.
            = =========================================================
       
        Return Value:
            
            This method returns the polynomial evaluated Z=P(X,Y). Z will 
            have the same shape as X and Y 
            
        In order to increase speed, use peval to evaluate for a single point,
        or meval to evaluate for 2D matrices
        '''
        Result=np.zeros(np.array(x).shape, dtype=np.double)
        cdef int i
        for i in self.cohef.nonzero()[0]:        
            Result=Result+self.cohef[i]*np.power(x, self.px_c[i])*np.power(y, self.py_c[i])
        return Result


    ###ojr
    cpdef eval_2(self, x, y, key):
        Result=np.zeros(np.array(x).shape, dtype=np.double)
        cdef int i
        for i in key:
            Result=Result+self.cohef[i]*np.power(x, self.px_c[i])*np.power(y, self.py_c[i])
        return Result


    @cython.boundscheck(False) # turn of bounds-checking for entire function
    @cython.wraparound(False)
    cpdef double peval(self, double x, double y):
        '''
        Evaluate the polynomial at the point x,y
              
        Arguments:
            = =========================================================
            x Floating point number containing the x value used to 
              evaluate the polynomial.
            y Floating point number containing the y value used to 
              evaluate the polynomial.
            = =========================================================
       
        Return Value:
            
            This method returns the polynomial evaluated z=P(X,Y).  
            
        '''
        cdef int i,j,nx,ny,ord,co,k
       
        
        cdef double Result=0.
        
        cdef double *X = <double *>malloc((self.order+1)*sizeof(double))
        cdef double *Y = <double *>malloc((self.order+1)*sizeof(double))
        
        X[0]=1
        Y[0]=1
        #Calculate the powers (it is faster this way than using pow from math.h) 
        for ord in range(self.order):
            X[ord+1]=X[ord]*x
            Y[ord+1]=Y[ord]*y
        for k in range(self.clen):
            Result +=  self.cohef_c[k]*X[<int>self.px_c[k]]*Y[<int>self.py_c[k]];
                #   Result[i,j] +=  self.cohef_c[k]*pow(x[i,j],self.px[k])*pow(y[i,j],self.py[k]);
                # Calculating the exponentials using products is a lot faster
        free(X)
        free(Y)
        return Result

    @cython.boundscheck(False) # turn of bounds-checking for entire function
    @cython.wraparound(False)
    def meval(self, np.ndarray[np.float64_t, ndim=2]x, np.ndarray[np.float64_t, ndim=2]y):
        '''
        Evaluate the polynomial at for the values given in the 2D matrices
        x, y.
        
        Arguments:
            = =========================================================
            x 2D matrix the values of x where the polynomial is going to 
              be evaluated.
            y 2D matrix the values of y where the polynomial is going to 
              be evaluated.
            = =========================================================
       
        Return Value:
            
            This method returns the polynomial evaluated Z=P(X,Y). Z will 
            have the same shape as x and y 
        '''
        cdef int i,j,nx,ny,ord,co,k
       
        nx=x.shape[0]
        ny=x.shape[1]
        
        #cdef np.ndarray[np.float64_t, ndim=2] Result=np.zeros(np.array(x).shape, dtype=np.double)
        cdef np.ndarray[np.float64_t, ndim=2] Result=zero_mat(nx, ny )
       
        cdef double *X = <double *>malloc((self.order+1)*sizeof(double))
        cdef double *Y = <double *>malloc((self.order+1)*sizeof(double))
        
        for i in range(nx):
            for j in range(ny):
                X[0]=1
                Y[0]=1
                #Calculate the powers (it is faster this way than using pow from math.h) 
                for ord in range(self.order):
                    X[ord+1]=X[ord]*x[i,j]
                    Y[ord+1]=Y[ord]*y[i,j]
                for k in range(self.clen):
                    Result[i,j] +=  self.cohef_c[k]*X[<int>self.px_c[k]]*Y[<int>self.py_c[k]];
                #   Result[i,j] +=  self.cohef_c[k]*pow(x[i,j],self.px[k])*pow(y[i,j],self.py[k]);
                # Calculating the exponentials using products is a lot faster
        free(X)
        free(Y)        
        
        return Result
        
    @cython.boundscheck(False) # turn of bounds-checking for entire function
    @cython.wraparound(False)
    def mevalr(self, np.ndarray[np.float64_t, ndim=2]x, np.ndarray[np.float64_t, ndim=2]y, double rot=0):
    
        '''
        Evaluate the polynomial at for the values given in the 2D matrices
        x, y rotated.
        
        Arguments:
            === =========================================================
            x   2D matrix the values of x where the polynomial is going to 
                be evaluated.
            y   2D matrix the values of y where the polynomial is going to 
                be evaluated.
            rot Angle to rotate the coordinate points x,y before evaluating
                the polynomial. Given in radians
            === =========================================================
       
        Return Value:
            
            This method returns the polynomial evaluated Z=P(x', y'). Z will 
            have the same shape as x and y. x' and y' are calculated from
            rotating x and y
        '''
        cdef int i,j,nx,ny,ord,co,k
       
        #Calculate the rotation matrix
        #Cos  -sin
        #sin cos
        
        cdef double cosr= cos(rot)
        cdef double sinr= sin(rot)
        
        cdef double rx,ry
        
        nx=x.shape[0]
        ny=x.shape[1]
        
        #cdef np.ndarray[np.float64_t, ndim=2] Result=np.zeros(np.array(x).shape, dtype=np.double)
        cdef np.ndarray[np.float64_t, ndim=2] Result=zero_mat(nx, ny )
        
        cdef double *X = <double *>malloc((self.order+1)*sizeof(double))
        cdef double *Y = <double *>malloc((self.order+1)*sizeof(double))
        
        for i in range(nx):
            for j in range(ny):
                X[0]=1
                Y[0]=1
                #Calculate the powers (it is faster this way than using pow from math.h) 
                for ord in range(self.order):
                    rx= x[i,j]*cosr-y[i,j]*sinr
                    ry= x[i,j]*sinr+y[i,j]*cosr
                    
                    X[ord+1]=X[ord]*rx
                    
                    Y[ord+1]=Y[ord]*ry
                for k in range(self.clen):
                    Result[i,j] +=  self.cohef_c[k]*X[<int>self.px_c[k]]*Y[<int>self.py_c[k]];
                #   Result[i,j] +=  self.cohef_c[k]*pow(x[i,j],self.px[k])*pow(y[i,j],self.py[k]);
                # Calculating the exponentials using products is a lot faster
        free(X)
        free(Y)        
        return Result
    @cython.boundscheck(False) 
    @cython.wraparound(False)   
    def vveval(self, np.ndarray[np.float64_t, ndim=1]x, np.ndarray[np.float64_t, ndim=1]y):
        '''
        Evaluate the polynomials in a 2D mesh defined by the vectors x and y.
        
        Arguments:
            === =========================================================
            x   1D vector with the values of x where the polynomial is going to 
                be evaluated.
            y   1D vector with the values of y where the polynomial is going to 
                be evaluated.
            === =========================================================
       
        Return Value:
            
            This method returns the polynomial evaluated Z=P(x,y). Z will 
            have the shape (nx, ny), where nx is the length of the x vector, 
            and ny is the length of the y vector.
            
        This method is not really faster than meval, but it might be 
        useful when doing parallel programming parallelizing, because needs 
        to transfer a smaller amount of information to the child processes 
        (vectors instead matrices).
        '''
                
        cdef int i,j,nx,ny,ord,co,k
       
        nx=x.shape[0]
        ny=y.shape[0]
        
        #cdef np.ndarray[np.float64_t, ndim=2] Result=np.zeros((nx,ny), dtype=np.double)
        cdef np.ndarray[np.float64_t, ndim=2] Result=zero_mat(nx, ny )
        
        #####Calculate X**n and Y**n
        cdef double **X = <double **>malloc((self.order+1)*sizeof(double *))
       
        cdef double **Y = <double **>malloc((self.order+1)*sizeof(double *))
        
        #X^0 Y^0
        cdef double *Xp = <double *>malloc(nx*sizeof(double))
        cdef double *Yp = <double *>malloc(ny*sizeof(double))
        
       
        for i in range(nx):
            Xp[i]=1.
            
        for i in range(ny):
            Yp[i]=1.
        X[0]=Xp
        Y[0]=Yp
        
       
        for ord in range(self.order):
            Xp = <double *>malloc(nx*sizeof(double))
            Yp = <double *>malloc(ny*sizeof(double))
           
            X[ord+1]=Xp
            Y[ord+1]=Yp
            
            for i in range(nx):
                #X[ord+1][i]=X[ord][i]*x[i]
                Xp[i]=X[ord][i]*x[i]
                
            for i in range(ny):
                #Y[ord+1][i]=Y[ord][i]*y[i]
                Yp[i]=Y[ord][i]*y[i]
        
        #######             
        for j in range(nx):
            for i in range(ny):
                for k in range(self.clen):
                    Result[j,i] +=  self.cohef_c[k]*X[<int>self.px_c[k]][i]*Y[<int>self.py_c[k]][j];
                    #Result[j,i]=X[1][i]
        for i in range(self.order):
            free(X[i])
            free(Y[i])
        free(X)
        free(Y)        
        return Result
        

    def gpu_eval(self,np.ndarray[np.float64_t, ndim=1]x, np.ndarray[np.float64_t, ndim=1]y):
        '''
        Evaluate the polynomial using the GPU, x and y are vectors
        
        
        This method is similar to vveval, but using the GPU
        '''
        
        mf = cl.mem_flags


        nx=x.shape[0]
        ny=y.shape[0]
        
       
        #X=x.astype(np.float64)
        #Y=y.astype(np.float64)
        #PX=self.px64 # The conversion to float64 is made in __init__
       # PY=self.py64
        
        
        CO=self.cohef.astype(np.float64)
        
        #Create the object field buffers for the GPU
        Xb = cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=x)
        Yb = cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=y)
        
        PXb=cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=self.px64)
        PYb=cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=self.py64)
        
        COb=cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=CO)
        #Create the buffer for the result
        #res=np.empty((nx,ny),dtype=np.float64)
        res= empty_mat( ny, nx )
    
        # Create the result buffer for the GPU
        
        RESb=cl.Buffer(ctx0, mf.WRITE_ONLY,res.nbytes)

        prg0.rsk(queue0, (nx,ny),
                    Xb,
                    Yb,
                    PXb,
                    PYb,
                    COb,
                    np.uint32(len(CO)),
                    np.uint32(self.order),
                    RESb,
                    local_size=(16,16)
                    )
        
        
        cl.enqueue_read_buffer(queue0, RESb, res).wait()
        return res
            
    def gpu_eval1(self,samples):
        "Evaluate the polynomial between -1 and 1 using the gpu"
        
        mf = cl.mem_flags


        nx,ny=samples
        
        # Translate the info to floats. Latter I'll check it using doubles
       
        PX=self.px64 # The conversion to float32 is made in __init__
        PY=self.py64
        
        
        CO=self.cohef.astype(np.float64)
        
        #Create the object field buffers for the GPU
        
        PXb=cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=PX)
        PYb=cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=PY)
        
        COb=cl.Buffer(ctx0, mf.READ_ONLY | mf.USE_HOST_PTR, hostbuf=CO)
        #Create the buffer for the result
        res=np.empty((nx,ny),dtype=np.float64)
 
        # Create the result buffer for the GPU
        
        RESb=cl.Buffer(ctx0, mf.WRITE_ONLY,res.nbytes)

        prg1.rsk(queue0, (nx,ny),
                    PXb,
                    PYb,
                    COb,
                    np.uint32(len(CO)),
                    np.uint32(self.order),
                    RESb,                  
                    local_size=(16,16)
                    )
        
        
        cl.enqueue_read_buffer(queue0, RESb, res).wait()
        return res

cpdef i2pxpy(i):
    """Method that returns the x and y powers for a given poly index
    
    ========= = == == == == == == == == == == == == == == ===
    index     0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 ...
    xpower    0  1  0  2  1  0  3  2  1  0  4  3  2  1  0 ...
    ypower    0  0  1  0  1  2  0  1  2  3  0  1  2  3  4 ...
    pol_order 0  1  1  2  2  2  3  3  3  3  4  4  4  4  4 ...
    ========= = == == == == == == == == == == == == == == ===
    """
    # Calculate the polynomial order
    
    cdef np.ndarray ia=np.array(i)
    poly_order=(1+(np.sqrt(8*(ia)+1)-3)/2).astype(int)
    ret_y=np.where(poly_order==0, 0,\
                ia-(((poly_order-1)**2+3*(poly_order-1)+2)/2 - 1)-1).astype(int)
    ret_x=poly_order-ret_y
    return ret_x, ret_y
cpdef int pxpy2i(int px,int py):
    """Method that returns the index given power in x and the power in y
    """
    cdef int po=px+py
    cdef int i0=(((po-1)*2+3)**2-1)/8
    return i0+py
    
cpdef ord2i(o):
    """
    Method that returns the number of coeficients of a polynomial of order o
    
    ===== == == == == == ===
    order  0  1  2  3  4 ...
    i      1  3  6 10 15 ...
    ===== == == == == == ===
    """
    return (o+2)*(o+1)/2
