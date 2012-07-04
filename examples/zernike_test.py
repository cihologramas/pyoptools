# -*- coding: utf8 -*-


### Example to use the taylor form of the zernike polynomials



from pyoptools.all import *
from pylab import *
from numpy.ma import array as marray


#Define the X and Y coordinates to evaluate the polynomial


X,Y=meshgrid(linspace(-1,1,512),linspace(-1,1,512))
mask=where(X**2+Y**2<1,False,True)
pz=zernike2taylor(7,-7)
z=marray(pz.meval(X,Y),mask=mask)
imshow(z)
 
