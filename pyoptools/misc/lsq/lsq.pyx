from numpy import array, sqrt, zeros, where, power, dot,  ones,  empty, inf
from numpy.linalg import inv, pinv, cond,  solve
from pyoptools.misc.Poly2D.Poly2D cimport *
#from pyoptools.misc.Poly2D import *

cimport numpy as np
import numpy as np
cimport cython

#from openopt import DFP


@cython.boundscheck(False)
@cython.wraparound(False)
def polyfit2d(x, y,z, int order=2):
    cdef int nc,nd,p
    nc= (order+2)*(order+1)/2 #_ord2i_(order)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] xa=array(x)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ya=array(y)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] za=array(z)
    cdef np.ndarray[np.int_t, ndim=1, mode="c"] px,py
    
    px, py=i2pxpy(range(0, nc))
    nd=len(xa)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] XP=ones((nd, 2*(order+1)))
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] YP=ones((nd, 2*(order+1)))
    for p in range(1, 2*(order+1)):
        XP[:, p]=XP[:, p-1]*xa
        YP[:, p]=YP[:, p-1]*ya

    #Calculating X and Y powers
    cdef np.ndarray[np.int_t, ndim=2, mode="c"] potx=empty((nc, nc),dtype=np.int)
    cdef np.ndarray[np.int_t, ndim=2, mode="c"] poty=empty((nc, nc),dtype=np.int)
    potx[:, :]=px
    potx=potx+potx.T
    
    poty[:, :]=py
    poty=poty+poty.T

    #Creating the Vandermonde matrix
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] mat=empty((nc, nc))
    #The python code was changed to a C++ inline
    #for  ir in range(nc):
    #    for ic in range(nc):
    #        powx=potx[ir, ic]
    #        powy=poty[ir, ic]
    #        tv=XP[:, powx]*YP[:, powy]
    #        mat[ir, ic]=tv.sum()
    cdef int powx,powy,ir,ic, id
    cdef double sum
    for ir in range(nc):
        for ic in range(nc):
            powx=potx[ir,ic];
            powy=poty[ir,ic];
            sum=0;
            for id in range(nd): 
                sum=sum+XP[id,powx]*YP[id,powy]
            mat[ir,ic]=sum;
                    
    
    
    #imat=pinv(mat)

    cdef np.ndarray[np.double_t, ndim=2, mode="c"] vec=empty((nc, 1))
    #The python code was changed to C++ inline
    #for ic in range (nc):
    #    tv=XP[:, px[ic]]*YP[:, py[ic]]*z
    #    vec[ic, 0]=tv.sum()

    for ic in range(nc):
        sum=0;
        powx=px[ic]
        powy=py[ic]
        for id in range(nd):
            sum=sum+za[id]*XP[id,powx]*YP[id,powy]
        vec[ic,0]=sum;
    
    cohef= solve(mat, vec)
    ret_poly=poly2d(cohef[:,0])

    #Calculate error. Verify if this is the best way
    e= sqrt(power(array(z)-ret_poly.eval(xa, ya), 2).mean())
        
    return ret_poly, e
    
    
@cython.boundscheck(False)
@cython.wraparound(False)
def vander_matrix(x, y,z, int order=2):
    cdef int nc,nd,p
    nc= (order+2)*(order+1)/2 #_ord2i_(order)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] xa=array(x)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] ya=array(y)
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] za=array(z)
    cdef np.ndarray[np.int_t, ndim=1, mode="c"] px,py
    
    px, py=i2pxpy(range(0, nc))
    nd=len(xa)
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] XP=ones((nd, 2*(order+1)))
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] YP=ones((nd, 2*(order+1)))
    for p in range(1, 2*(order+1)):
        XP[:, p]=XP[:, p-1]*xa
        YP[:, p]=YP[:, p-1]*ya

    #Calculating X and Y powers
    cdef np.ndarray[np.int_t, ndim=2, mode="c"] potx=empty((nc, nc),dtype=np.int)
    cdef np.ndarray[np.int_t, ndim=2, mode="c"] poty=empty((nc, nc),dtype=np.int)
    potx[:, :]=px
    potx=potx+potx.T
    
    poty[:, :]=py
    poty=poty+poty.T

    #Creating the Vandermonde matrix
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] mat=empty((nc, nc))
    #The python code was changed to a C++ inline
    #for  ir in range(nc):
    #    for ic in range(nc):
    #        powx=potx[ir, ic]
    #        powy=poty[ir, ic]
    #        tv=XP[:, powx]*YP[:, powy]
    #        mat[ir, ic]=tv.sum()
    cdef int powx,powy,ir,ic, id
    cdef double sum
    for ir in range(nc):
        for ic in range(nc):
            powx=potx[ir,ic];
            powy=poty[ir,ic];
            sum=0;
            for id in range(nd): 
                sum=sum+XP[id,powx]*YP[id,powy]
            mat[ir,ic]=sum;
                    
    
    
    #imat=pinv(mat)

    cdef np.ndarray[np.double_t, ndim=2, mode="c"] vec=empty((nc, 1))
    #The python code was changed to C++ inline
    #for ic in range (nc):
    #    tv=XP[:, px[ic]]*YP[:, py[ic]]*z
    #    vec[ic, 0]=tv.sum()

    for ic in range(nc):
        sum=0;
        powx=px[ic]
        powy=py[ic]
        for id in range(nd):
            sum=sum+za[id]*XP[id,powx]*YP[id,powy]
        vec[ic,0]=sum;
    return mat,vec
    

#Local definitions for the optimized order 2 fitter

_potxo2=array([[ 0.,  1.,  0.,  2.,  1.,  0.], 
                [ 1.,  2.,  1.,  3.,  2.,  1.], 
                [ 0.,  1.,  0.,  2.,  1.,  0.], 
                [ 2.,  3.,  2.,  4.,  3.,  2.], 
                [ 1.,  2.,  1.,  3.,  2.,  1.], 
                [ 0.,  1.,  0.,  2.,  1.,  0.]])
_potyo2=array([[ 0.,  0.,  1.,  0.,  1.,  2.], 
                [ 0.,  0.,  1.,  0.,  1.,  2.], 
                [ 1.,  1.,  2.,  1.,  2.,  3.], 
                [ 0.,  0.,  1.,  0.,  1.,  2.], 
                [ 1.,  1.,  2.,  1.,  2.,  3.], 
                [ 2.,  2.,  3.,  2.,  3.,  4.]])
_pxo2, _pyo2=(array([0, 1, 0, 2, 1, 0]), array([0, 0, 1, 0, 1, 2]))

def polyfito2(x, y, z):
    """
    Polyfit function optimized for polynomials of order 2
    """
    nc= 6  #_ord2i_(order)
    xa=array(x)
    ya=array(y)
    px, py= _pxo2,_pyo2# _i2pxpy_(range(0, nc))
    nd=len(xa)
    XP=ones((nd, 6)) #ones((nd, 2*(order+1)))
    YP=ones((nd, 6))
    #for p in range(1, 6):
    #    XP[:, p]=XP[:, p-1]*xa
    #    YP[:, p]=YP[:, p-1]*ya
    cdef int p,id
    for p in range (1,6):
        for id in range(nd):
            XP[id,p]=XP[id,p-1]*xa[id];
            YP[id,p]=YP[id,p-1]*ya[id];

    #Calculating X and Y powers
    #potx=empty((nc, nc))
    #poty=empty((nc, nc))
    #potx[:, :]=px
    #potx=potx+potx.T
    #print "*",potx
    #poty[:, :]=py
    #poty=poty+poty.T
    #print "**", poty
    potx=_potxo2
    poty=_potyo2
    #Creating the Vandermonde matrix
    mat=empty((nc, nc))
    #The python code was changed to a C++ inline
    #for  ir in range(nc):
    #    for ic in range(nc):
    #        powx=potx[ir, ic]
    #        powy=poty[ir, ic]
    #        tv=XP[:, powx]*YP[:, powy]
    #        mat[ir, ic]=tv.sum()
    cdef int powx,powy
    cdef int ir,ic
    cdef double sum
    for ir in range[nc]:
        for ic in range (nc):
            powx=potx[ir,ic];
            powy=poty[ir,ic];
            sum=0;
            for id in range(nd):
                sum=sum+XP[id,powx]*YP[id,powy];
            mat[ir,ic]=sum;
        
    
    
    #imat=pinv(mat)

    vec=empty((nc, 1))

    for ic in range (nc):
        tv=XP[:, px[ic]]*YP[:, py[ic]]*z
        vec[ic, 0]=tv.sum()

    cohef= solve(mat, vec)
    ret_poly=poly2d(cohef, px=px, py=py)

    #Calculate error. Verify if this is the best way
    ep=cohef[0]+\
       cohef[1]*xa+\
       cohef[2]*ya+\
       cohef[3]*xa*xa+\
       cohef[4]*xa*ya+\
       cohef[5]*ya*ya
       
    e= sqrt(power(array(z)-ep, 2).mean())
        
    return ret_poly, e



##Se elimina este código pues no debe ser importante, y esta usando weave
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


#~ def test(niter=10,Omax=4,  ndat=500):
    #~ from numpy.random import rand
    #~ for n in range(niter):
        #~ ncohef=ord2i(Omax)
        #~ cohef=100*rand(ncohef)
        #~ X=[] 
        #~ Y=[]
        #~ Z=[]
        #~ p=poly2d(cohef)
        #~ for j in range(ndat):
            #~ x=100.*rand()
            #~ y=100.*rand()
            #~ X.append(x)
            #~ Y.append(y)
            #~ Z.append(p.eval(x, y))
        #~ co, er=polyfit(X, Y, Z, Omax)
        #~ print cohef
        #~ print sqrt(power(cohef-co.cohef.T, 2).sum())/ncohef
        #~ print "*"
#~  Test using the open opt library it seems not to be needed
#~ def polyfitOO(x, y, z, order):
    #~ 
    #~ def f(C, D):
        #~ x=D[0]
        #~ y=D[1]
        #~ a=poly2d(cohef=C)
        #~ return a.eval(x, y)
    #~ ND=_ord2i_(order)
    #~ p, e=polyfit(array(x), array(y), array(z), order=order)
    #~ C=p.cohef
    #~ print C
    #~ D=array([x, y]).transpose()
    #~ Z=array(z).transpose()
    #~ lb = ND*[-inf]
    #~ ub = ND*[inf]
    #~ p = DFP(f, C , D, Z, lb=lb, ub=ub)
    #~ r = p.solve('nlp:ralg', plot=0, iprint = 10)
    #~ #print C
    #~ #print 'solution: '+str(r.xf)+'\n||residuals||^2 = '+str(r.ff)
    #~ ret_poly=poly2d(r.xf)
    #~ xa=array(x)
    #~ ya=array(y)
    #~ e= sqrt(power(array(z)-ret_poly.eval(xa, ya), 2).mean())
    #~ return ret_poly, e
