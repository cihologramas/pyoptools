#!/usr/bin/env python
# -*- coding: UTF-8 -*-


'''simulacion del litografo

'''
from pyoptools import *


# Definition of the ray sources at the origin

wavelength=0.442

def error(d,graphics=False):
    r_b=point_source_p(origin=(0.,0.,0.),direction=(0.,0.,0),span=pi/5
                       ,num_rays=(15,15),wavelength=wavelength, label="blue")





    r2=poly2d((0,0,0,1,0,1))
    r4=r2*r2
    r6=r4*r2
    r8=r4*r4
    r10=r8*r2


    #edmund 47146

    #surface 1
    K=-0.533348;    A4=4.202096*10**-4;    A6=-7.493867*10**-5;      A8=2.403049*10**-5;     A10=-3.180361*10**-6
    poly1=A4*r4+ A6*r6 +A8*r8 +A10*r10
    asf1=Aspherical(Kx=K, Ky=K, Ax=1/2.774,Ay=1/2.774,shape=Circular(radius=2.7), 
                     poly=poly1
                    )  
    #Surface 2

    K=-191.5207;     A4=-1.443363*10**-3;     A6=1.4648185*10**-3;       A8=-4.7042041*10**-4;    A10=4.8133394*10**-5

    poly2=A4*r4+ A6*r6 +A8*r8 +A10*r10
    asf2=Aspherical(Kx=K, Ky=K, Ax=1/-14.582,Ay=1/-14.582,shape=Circular(radius=2.7), 
                     poly=poly2
                    )  

    ASF=Component(surflist=[(asf1, (0, 0, 0), (0, 0, 0)), 
                           (asf2, (0, 0, 3.04), (0, 0, 0))
                           ],
                            material=1.62068) #Need to find ECO-550
                            
                            
    #ECO550     w           n
    #        407.7      1.62636
    #        435.8      1.62068
    #       480.0       1.61464
    ccd=CCD(size=(10,10))
    os=System(complist=[(ASF,(0,0,d[0]),(0,pi,0)),
                        (ccd,(0,0,20),(0,0,0)),#Distancia encontrada buscando la interseccion entre 2 rayos cualquiera
                                                                    # del ccd.hitlist
                        ],n=1)





    #Add the ray sources

    os.ray_add(r_b)
    os.propagate()

    if graphics:
        spot_diagram(ccd)
        glPlotFrame(os)

    X,Y,Z=ccd.get_optical_path_data()   
    print d,array(Z).std()

    return array(Z).std()

#from scipy.optimize import fmin
#x0=fmin(error,5)    
x0=[5.1]
print error(x0,True)
