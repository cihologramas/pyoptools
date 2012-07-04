from sys import exit, path
path.append("../")

from pyoptools import *
from time import time
from pylab import *


d=zeros((512,512),dtype=complex)
psize=5./512
a=1.
b=0.5
A=1/psize
B=.5/psize
Z=70.

print (a)**2/(.632e-3*Z)

d[256-A/2.:256+A/2.,256-B/2.:256+B/2.]=1


U=Field(data=d, psize=psize, l=.632e-3)

sx,sy=U.size
figure();imshow(U.intensity(),extent=(-sx/2,sx/2,-sy/2,sy/2));colorbar()
print sx,sy

FU=U.propagate_fourier(Z) #Ojo el nombre en la libreria esta mal
sx,sy=FU.size
print sx,sy
figure();imshow(FU.intensity(),extent=(-sx/2,sx/2,-sy/2,sy/2));colorbar()
title("Fraunhoffer")

FU=U.propagate_ae(Z,1.) #Ojo el nombre en la libreria esta mal
sx,sy=FU.size
print sx,sy
figure();imshow(FU.intensity(),extent=(-sx/2,sx/2,-sy/2,sy/2));colorbar()
title("Espectro Angular")

FU=U.propagate_fresnel(Z) #Ojo el nombre en la libreria esta mal
sx,sy=FU.size
print sx,sy
figure();imshow(FU.intensity(),extent=(-sx/2,sx/2,-sy/2,sy/2));colorbar()
title("Fresnel")


#imshow(angle(ff))
#figure()
#imshow(ff.intensity());colorbar()

#figure()
#imshow(ph.intensity());colorbar()

show()
