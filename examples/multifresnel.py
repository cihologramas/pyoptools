from sys import exit, path
path.append("../")

from pyoptools import *
from time import time
from pylab import *
#???7.5um pinhole
#5um pinhole
s=2048
F=100.
phd=zeros((s,s))
phd[s/2-100:s/2+100,s/2-100:s/2+100]=1
#phd[s/2-100:s/2-80,s/2-100:s/2+100]=1
#phd[s/2-100:s/2+100,s/2-100:s/2-80]=1
#phd[s/2-100:s/2+100,s/2-10:s/2+10]=1
ph=Field(data=phd, psize=2.0e-3, l=.6328e-3)
X,Y= ph.field_sample_coord
ph1=ph.propagate_ae(F,1)
figure()
imshow(ph1.angle[::5,::5]);colorbar()

del(ph)


lens0=where(X**2+Y**2 < 4., exp(-2.j*pi/.6328e-3*(sqrt(X**2+Y**2+F**2)+(X+Y)*0.22/F)),0) #+ exp(-((X+0.01)**2+Y**2)/.0025**2) 
lens1=where(X**2+Y**2 < 4., exp(.5j*pi-2.j*pi/.6328e-3*(sqrt(X**2+Y**2+F**2)+(-X+Y)*0.22/F)),0) #+ exp(-((X+0.01)**2+Y**2)/.0025**2) 
lens2=where(X**2+Y**2 < 4., exp(1.j*pi+-2.j*pi/.6328e-3*(sqrt(X**2+Y**2+F**2)+(X-Y)*0.22/F)),0) #+ exp(-((X+0.01)**2+Y**2)/.0025**2) 
lens3=where(X**2+Y**2 < 4., exp(1.5j*pi-2.j*pi/.6328e-3*(sqrt(X**2+Y**2+F**2)+(-X-Y)*0.22/F)),0) #+ exp(-((X+0.01)**2+Y**2)/.0025**2) 


lens4=exp(-2.j*pi/.6328e-3*(sqrt(X**2+Y**2+F**2)))

#lens=exp(1.j *angle((lens0+lens1+lens3+lens4)/4))
lens=(lens0+lens1+lens2+lens3)/4
del(lens0)
del(lens1)
figure()
imshow(abs(lens));colorbar()

figure()
imshow(angle(lens));colorbar()

ph2=ph1*lens
del(lens)
del(ph1)

ph3=ph2.propagate_ae(2*F,1)
ph3=ph3*lens4
ph3=ph3.propagate_ae(F,1)

del(ph2)

figure()
imshow(ph3.angle[::5,::5]);colorbar()

figure()
imshow(ph3.intensity()[::5,::5]);colorbar()

#~ #imshow(angle(ff))

#~ t=time()
#~ ff=ph.propagate_rsc_sc(250,scale=(2,2))
#~ #ff=ph.propagate_rsc(250)
#~ #ff=ph._rs_kernel(x=X+100e-3, y=Y, z=250., n=1.)
#~ print time() -t
#~ #rfd=ph.propagate_ae(0.5,n=1.)
#~ figure()
#~ imshow(ff.angle);colorbar()
#~ #imshow(angle(ff))
#~ #figure()
#~ #imshow(ff.intensity());colorbar()
#~ 
#~ #figure()
#~ #imshow(ph.intensity());colorbar()
#~ 
#~ show()
#~ 
#~ """ Mirar si usando esto se puede hacer un algoritmo para calcular la 
 #~ propagacion a partir de muchas peque#as. Esto se puede ademas usar para
 #~ hacer cambio de escala (tamano de matriz de salida), y mirar si se puede 
 #~ crear un algoritmo donde el campo optico se pueda guardar como algo 
 #~ multiresolucion"""
show()
