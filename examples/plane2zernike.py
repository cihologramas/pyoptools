#!/usr/bin/python
# -*- coding: utf8 -*-


### Primera prueba de conversion de una onda plana a un polinomio de zernike.


from sys import exit, path
path.append("../")

from pyoptools import *
from time import time
from pylab import *
import scipy.integrate as integrate

dmax=1.05
dx=dmax/500.
r,t=polar_array(dmax,dx)
r=r[:-1,:-1]
t=t[:-1,:-1]

#ps=Plane(shape=Rectangular(size=(25,15)))
#cs=Plane(shape=Circular(radius=25))

cs=Spherical(shape=Circular(radius=15),curvature=-1/50.)

#Encontrar el z de el borde de la lente superficie.
X,Y,Z=cs.shape.mesh(ndat=(100,100),  topo=cs.topo)
im=X[50].argsort()[0]
print im
zm=Z[50][im]
print zm
I1, d1,er=cs.pw_propagate(Ray(dir=(0.5, 0, 1)),1,1.5, rsamples=(800, 800) , shape=(1000, 1000), knots=10,z=zm)
print "**",er
#figure();imshow(I);colorbar()
#figure();imshow(d);colorbar()
cohef=[]
for n in range(0,10,1):
    for m in range(-n,n+1,2):
        zpol=zernike(n,m,r,t)
        #Note the integral on masked arrays have problems. It does not assume that
        #the mask implies 0 value
        id=d1*zpol
        mask=(~ma.getmaskarray(id)).astype(int)
        id=id*mask

        
        idata=integrate.trapz(integrate.trapz((1./pi)*id,dx=dx),dx=dx)
        #figure();
        #imshow(d1*zpol);colorbar()
        #idata=(1./pi)*(data*zpol).sum()*dx*dx
        cohef.append((n,m,idata))
        print n,"\t",m,"\t",idata
#I2, d2=cs.pw_propagate(Ray(dir=(.3, 0, 1)),1.5,1,  rsamples=(200, 200) ,  isamples=(400, 400))
#print d
exp_pol=zeros(r.shape).astype(float)
for i in cohef:
    n,m,c=i
    exp_pol=exp_pol+c*zernike(n,m,r,t)
    print "C%2d,%2d = % 5e error= % 5e"%(n,m,c,(d1-exp_pol).std())

figure();
title("fase reconstruida")
imshow(exp_pol);colorbar()
figure();
title("fase de los rayos")
imshow(d1);colorbar()
figure()
imshow(d1-exp_pol);colorbar()
show()
