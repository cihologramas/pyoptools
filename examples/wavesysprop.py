#!/usr/bin/python
# -*- coding: utf8 -*-


'''Example of a spherical lens simulation

'''
from sys import exit, path
path.append("../")

from pyoptools import *
from time import time
from pylab import *
from wavefront.field_old import Field as Fieldo



#### Get the surface hit secuence by propagating the axial ray thru the optical system


BK7=schott['BK7']
SF5=schott['SF5']

l1=SphericalLens(radius=3, thickness=4.93, curvature_as=1./257.57, 
                 curvature_ps=-1./257.57, material= BK7 )
                 

ccd=CCD(size=(3,3))

os=System(complist=[(l1,(0,0,250),(0,0,0)),
                     #(ccd,(0,0,300),(0,0,0)),
                    ],n=1)
                    
#Definicion del rayo axial
iray=Ray()

os.ray_add(iray)
os.propagate()
                 
raylist=[]
nray=iray
while len(nray.childs)>0:
    raylist.append(nray)
    nray=nray.childs[0]                 

# Get the surfaces and the positions 
prop_pos=[]
s_list=[]
n_list=[]

for r in raylist:
    prop_pos.append(r.pos[2])
    s_list.append(r.orig_surf)
    n_list.append(r.n)
prop_pos.append(r.childs[0].pos[2])
s_list.append(r.childs[0].orig_surf)
n_list.append(r.childs[0].n)
# Print the propagation information
for i in range(len(prop_pos)):
    print "%3.5f\t%40s\t%f"%(prop_pos[i],s_list[i],n_list[i])
    



#Simulación del pinhole

phd=zeros((128,128))
ph=Field(data=phd, psize=2.0e-3, l=.442e-3)

#ph=Field(phd, psize=2.0e-3, l=.442e-3)


X,Y= ph.field_sample_coord
ph.data=exp(-((X**2+Y**2)/.0025**2))

print "iniciando propagacion escalada"
t=time()

# para cubrir un diametro de 25.6 mm, se necesita escalar por (100,100)
#Por ahora vamos a suponer que la lente es de 6 mm y no de 25.
# Cuando todo este listo hay que cambiar esto en la definición del sistema optico

# Para probar el programa, se va a escalar por algo mas pequeño. Una vez este listo
# El programa, se hará la simulación completa


ff=ph.propagate_rsc_sc(250,scale=(25,25)) 
print "T=",time() -t
print "Finalizando propagacion escalada"

imshow(ff.angle)

print "Crear la primera superficie"
t=time()
psurf=PSurf(s_list[1],ni=1,nr=1.5,ilimit=0, slimit=0.1, step=0.03, rsamples=(512,512),l=.633e-3)# order=10, rsamples=(500,500),zb=None):
print "T=",time() -t
t=time()

pe=psurf.propagate(ff, pt=0.70, samples=(2048,2048))
print "T=",time() -t

figure()
imshow(pe.angle);colorbar()
show()
