#!/usr/bin/python
# -*- coding: utf8 -*-


### Primera prueba de conversion de una onda plana a un polinomio de taylor.


from sys import exit, path
path.append("../")

from pyoptools import *
from time import time
from pylab import *
import scipy.integrate as integrate
from numpy import meshgrid, arange, power, sqrt, polyfit, polyval
from math import atan
import numpy.ma as ma

#from scipy.misc import imrotate


#cs=Spherical(shape=Circular(radius=15),curvature=-1/50.)

cs=Plane(shape=Circular(radius=15),curvature=-1/50.)

#Encontrar el z de el borde de la lente superficie.

#Calcular los coheficientes.
ti=time()
psurf=PSurf(cs,ni=1,nr=1.5,ilimit=0, slimit=1.3, step=0.05, rsamples=(800,800))# order=10, rsamples=(500,500),zb=None):
print "**",time()-ti
zm=psurf.zm
#Crear un nuevo frente de onda para comparar con el interpolado
iang=.1

#x,y,d =cs.pw_propagate_list(Ray(dir=(0, tan(iang), 1)),1,1.5, rsamples=(800, 800),z=zm)
x,y,d =cs.pw_propagate_list(Ray(dir=(1, 1, 1)),1,1.5, rsamples=(800, 800),z=zm)

xm=15.

#Normalizar la pupila de entrada
x=x/xm; y=y/xm

pi,ei=polyfit2d(x, y, d,order=10)


xx=linspace(-1., 1., 500)
yy=linspace(-1., 1., 500)
I=hitlist2int(x, y, x,  xx, yy)

xx,yy=meshgrid(xx,yy)

rm=where(xx**2+yy**2>1,True,False)

d=ma.masked_array(pi.eval(xx,yy), mask=rm)

di= psurf.pw_evaluate((1,1,1), samples=(500,500))





print di[0]
print di[1]
figure()
imshow(d)
colorbar()

figure()
imshow(di[0])
colorbar()

figure()
imshow(d-di[0])
colorbar()

show()


#I, xe, ye=histogram2d(yi, xi, (xx, yy))
#        I=hitlist2int(xi, yi, xi,  xx, yy)

