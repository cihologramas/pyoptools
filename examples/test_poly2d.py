#!/usr/bin/python
# -*- coding: utf8 -*-

from pylab import *
from scipy.ndimage.interpolation import rotate 


'''Example of a spherical lens simulation

'''
from sys import exit, path
path.append("../")

from pyoptools import *
from time import time
from pylab import *




a=poly2d([0,1,2,4,5,6])
n=4096
xl=linspace(-1,1,n).astype(float64)
yl=linspace(-1,1,n).astype(float64)

x,y=meshgrid(xl,yl)

rotang=pi/2
#a1=a.eval(x,y)

a1=a.eval1(x,y)
a2=rotate(a1,rotang*180/pi,reshape=False)

a3=a.eval1r(x,y,rotang)


a4= a.gpu_evalr(xl,yl,rotang)

#a5= a.gpu_eval1((n,n))

figure();imshow(a1);colorbar()
figure();imshow(a2);colorbar()
figure();imshow(a3);colorbar()
figure();imshow(a4);colorbar()
#figure();imshow(a5); colorbar()

show()
