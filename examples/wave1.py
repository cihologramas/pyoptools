from sys import exit, path
path.append("../")

from pyoptools import *
from time import time
from pylab import *
#???7.5um pinhole
#5um pinhole
phd=zeros((128,128))
ph=Field(data=phd, psize=2.0e-3, l=.442e-3)
X,Y= ph.field_sample_coord
ph.data=exp(-(X**2+Y**2)/.0025**2) #+ exp(-((X+0.01)**2+Y**2)/.0025**2) 

t=time()
ff=ph.propagate_rsc_sc(250,scale=(2,2))
#ff=ph.propagate_rsc(250)
#ff=ph._rs_kernel(x=X+100e-3, y=Y, z=250., n=1.)
print time() -t
#rfd=ph.propagate_ae(0.5,n=1.)
figure()
imshow(ff.angle);colorbar()
#imshow(angle(ff))
#figure()
#imshow(ff.intensity());colorbar()

#figure()
#imshow(ph.intensity());colorbar()

show()

""" Mirar si usando esto se puede hacer un algoritmo para calcular la 
 propagacion a partir de muchas peque#as. Esto se puede ademas usar para
 hacer cambio de escala (tamano de matriz de salida), y mirar si se puede 
 crear un algoritmo donde el campo optico se pueda guardar como algo 
 multiresolucion"""
