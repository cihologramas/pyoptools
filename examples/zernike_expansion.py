# -*- coding: utf8 -*-


### Primera prueba de conversion de una onda plana a un polinomio de zernike.


from sys import exit, path
path.append("../")

from pyoptools import *
from time import time
from pylab import *
import scipy.integrate as integrate
import numpy.ma as ma
"""
Tarea de expandir un frente de onda en los polinomios de zernike
"""

    
data=load("./sombrero")
dmax=1.74
dx=dmax/240.
r,t=polar_array(dmax,dmax/240.)
r=r[:-1,:-1]
t=t[:-1,:-1]
import scipy.integrate as integrate
cohef=[]
for n in range(0,9,1):
    for m in range(-n,n+1,2):
        zpol=zernike(n,m,r,t)
        #Note the integral on masked arrays have problems. It does not assume that
        #the mask implies 0 value
        id=data*zpol
        mask=(~ma.getmaskarray(id)).astype(int)
        id=id*mask
        idata=integrate.trapz(integrate.trapz((1./pi)*id,dx=dx),dx=dx)
        #idata=(1./pi)*(data*zpol).sum()*dx*dx
        cohef.append((n,m,idata))

# Calcular el polinomio a partir de los cohef
exp_pol=zeros(r.shape).astype(float)
for i in cohef:
    n,m,c=i
    exp_pol=exp_pol+c*zernike(n,m,r,t)
    print "C%2d,%2d = % 5e error= % 5e"%(n,m,c,(data-exp_pol).std())
figure()
imshow(data-exp_pol)
colorbar()

figure()
imshow(data)
colorbar()

figure()
imshow(exp_pol)
colorbar()

#plot3d(data-exp_pol)
show()
