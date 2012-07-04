from time import time
from numpy import  hanning
ps=100
phd=zeros((ps, ps))
phd[ps/2-10:ps/2+10, ps/2-10:ps/2+10]=1
ph_o=Field(data=phd, psize=5e-3, l=.442e-3)


for z in [5,10,15,20,100,250]:
    ph=ph_o.propagate_rs(z,1)
    figure()
    imshow(ph.intensity())
    colorbar()
    figure()
    imshow(ph.angle)
    colorbar()





