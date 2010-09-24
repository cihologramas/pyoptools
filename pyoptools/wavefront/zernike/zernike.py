"""
Module defining the Zernike polynomials
"""

import types
from scipy import factorial, arange
from numpy import mgrid, sqrt, arccos, where, zeros,transpose, pi, cos, sin, ones
from numpy.ma import masked_array

def polar_array(Rmax=1.,DS=0.1, pr=1.):
    """
    Function that greates 2 square matrices one with rho, and the other with
    theta, to be able to calculate functions using polar coordinates.
    It is similar to the mgrid function.

    Arguments:
    
    Rmax
        Limit the pupil area -Rmax<=X<=Rmax -Rmax<=Y<=Rmax
    DS
        Step between pixels
    pr
        Pupil radius. Used to normalize the pupil.
        
    TODO: This function should be moved to a auxiliary functions module
    """
    
    X,Y= mgrid[-Rmax:Rmax+DS:DS,-Rmax:Rmax+DS:DS]/pr
    r = sqrt(X**2+Y**2)
    th= arccos(transpose(X*1./r))
    th= where(th<2.*pi,th,0)
    th= where(X<0,2.*pi-th,th)
    return r,th

def rnm(n,m,rho):
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
    
    if(type(n) is not int):
        raise Exception, "n must be integer"
    if(type(m) is not int):
        raise Exception, "m must be integer"
    if (n-m)%2!=0:
        raise Exception, "n-m must be even"
    if abs(m)>n:
        raise Exception, "The following must be true |m|<=n"
    mask=where(rho<=1,False,True)
    
    if(n==0 and m==0):
        return  masked_array(data=ones(rho.shape), mask=mask)
    rho=where(rho<0,0,rho)
    Rnm=zeros(rho.shape)
    S=(n-abs(m))/2
    for s in range (0,S+1):
        CR=pow(-1,s)*factorial(n-s)/ \
            (factorial(s)*factorial(-s+(n+abs(m))/2)* \
            factorial(-s+(n-abs(m))/2))
        p=CR*pow(rho,n-2*s)
        Rnm=Rnm+p
    return masked_array(data=Rnm, mask=mask)
    
def zernike(n,m,rho,theta):
    """
    Returns the an array with the Zernike polynomial evaluated in the rho and 
    theta
    
    Arguments:
    
    
    *n*
        n order of the Zernike polynomial
    
    *m*
        m order of the Zernike polynomial
        
    *rho*
        Matrix containing the radial coordinates. 
       
    *theta*
        Matrix containing the angular coordinates.
    
    Note: For rho>1 the returned value is 0
    
    Note: Values for rho<0 are silently returned as rho=0
    """
    
    
    Rnm=rnm(n,m,rho)
    
    NC=sqrt(2*(n+1))
    S=(n-abs(m))/2
    
    if m>0:
        Zmn=NC*Rnm*cos(m*theta)
    #las funciones cos() y sin() de scipy tienen problemas cuando la grilla
    # tiene dimension cero
    
    elif m<0:
        Zmn=NC*Rnm*sin(m*theta)
    else:
        Zmn=sqrt(0.5)*NC*Rnm
    return Zmn
    
##    
##def Ty_Mon(ypow,xpow,dim=-1):
##    n=xpow+ypow
##    if dim==-1:
##        Ty=N.zeros((n+1,n+1))
##        Ty[ypow,xpow]=1
##        return Ty
##    else:
##        Ty=N.zeros((dim+1,dim+1))
##        Ty[ypow,xpow]=1
##        return Ty
##
##
##
##def Ty_Mat2Vec(PolyMat):
##    """
##    Retorna un polinomio en forma de vector en la base de los Taylor
##    """
##    n=PolyMat.shape[0]
##    PolyVec=[]
##    for i in range(n):
##        for j in range(i+1):
##            x_pow=i-j
##            y_pow=j
##            PolyVec.append(PolyMat[y_pow,x_pow])
##    return N.asarray(PolyVec)
##    
##    
##
##
##def ZK2TY(n,m):
##    """
##    Retorna los coeficientes de Taylor del P. de Zernike dado 
##    """
##    TC=N.zeros((n+1,n+1))
##    mm=(n-m)/2
##    am=abs(m)
##    if m>0:
##        B=mm
##        if n%2==0:
##            p=1
##            q=(m/2)-1
##        else:
##            p=1
##            q=(m-1)/2
##    else:
##        B=n-mm
##        if n%2==0:
##            p=0
##            q=-(m/2)
##        else:
##            p=0
##            q=-(m+1)/2
##    for i in range(q+1):
##        for j in range(B+1):
##            for k in range(B-j+1):
##                c=pow(-1,i+j)*binomial(am,2*i+p)*binomial(B-j,k)*int(scipy.factorial(n-j))/(int(scipy.factorial(j))*int(scipy.factorial(mm-j))*int(scipy.factorial(n-mm-j)))
##                x_pow=2*(i+k)+p
##                y_pow=n-2*(i+j+k)-p
##                TC[y_pow,x_pow]=TC[y_pow,x_pow]+c
##    return TC
##
##
##def Eval_Poly(M,x,y):
##    x=N.asarray(x)
##    y=N.asarray(y)
##    Result=0
##    Pows=M.nonzero()
##    y_pow=Pows[0]
##    x_pow=Pows[1]
##    for i in range(y_pow.shape[0]):
##        Result=Result+M[y_pow[i],x_pow[i]]*pow(x,x_pow[i])*pow(y,y_pow[i])
##    return Result
##        
##def Eval_Poly_Comb(Ml,C,x,y):
##    M=C[0]*Ml[0]
##    for j in range(1,C.shape[0]):
##        M=M+C[j]*Ml[j]
##    return Eval_Poly(M,x,y)  
##    
##    
##
##def test(N,M):
##    R,T=Polar_Array(Rmax=1,DS=0.005)
##    P.figure(1)
##    i=1
##    for n in range(N):
##        for m in range(-M,M):
##            if abs(m)<=n and (n-abs(m))%2==0:
##                P.subplot(N,2*M,i)
##                Z=ZK_Poly(n,m,R,T)
##                P.imshow(Z)
##                P.axis("off")
##            i=i+1
##    P.savefig("zern.png",dpi=300)
##    P.show()
###test(3,3)



    
