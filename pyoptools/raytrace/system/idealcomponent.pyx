#!/usr/bin/env python
# -*- coding: UTF-8 -*-

# cython: profile=True

from numpy import inf, asarray, pi, alltrue, sometrue, isinf, isnan,array, dot, float64
cimport numpy as np
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.system import System
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.plane cimport Plane

from pyoptools.misc.picklable.picklable cimport Picklable

#from misc import rot_mat, rot_mat_i

from pyoptools.misc.cmisc.cmisc cimport * # This might be changed to a cimport. need to check

from pyoptools.misc.plist.plist cimport plist

class IdealThickLens(System):
    """
    Class defined to create an Ideal ThickLens.
    
    Is is created as a subsystem because the propagation inside has multiple rays
    entrance-> principal plane1
    principal plane1 -> principal plane2
    principal plane 2 -> exit
    fl -> focal lenght
    
    **ARGUMENTS**
    
    ===== ======================================================================
    shape Shape of the lens (Entrance and exit surface's shape).
    th    Thinckness of the lens (Distance between the entrance and exit
          surfaces)
    ppp   Tuple with the principal planes position. The position of each 
          principal plane is measured from its corresponding entrance surface. 
    ===== ======================================================================
          
    The origin of the component is located in the midle point between the 
    entrance and exit surface
    """

    def __init__(self,shape,th, pplanes =(0.,0.),f=100):

        System.__init__(self,n=1)

        self.__E1__ = Plane(shape=shape)
        self.__E2__ = Plane(shape=shape)
        self.__H1__ = Plane(shape=shape)
        self.__H2__ = Plane(shape=shape)

        self.__C1__ = Component(surflist = [(self.__E1__,(0,0,0),(0,0,0))])
        self.__C2__ = Component(surflist = [(self.__E2__,(0,0,0),(0,0,0))])
        self.__C3__ = Component(surflist = [(self.__H1__,(0,0,0),(0,0,0))])
        self.__C4__ = Component(surflist = [(self.__H2__,(0,0,0),(0,0,0))])

        self.complist["C1"]=(self.__C1__,(0,0,-th/2.),(0,0,0))
        self.complist["C2"]=(self.__C2__,(0,0,th/2.),(0,0,0))
        self.complist["H1"]=(self.__C3__,(0,0,-th/2.+pplanes[0]),(0,0,0))
        self.complist["H2"]=(self.__C2__,(0,0,th/2.+pplanes[1]),(0,0,0))

        self.pplanes = pplanes
        self.f=f


    def propagate_ray(self,Ray ri):

        pos = ri.pos
        dir = ri.dir
        wav = ri.wavelength
        label = ri.label
        parent = ri
        n=ri.n

        D,C,S = self.distance(ri)

        #Verificar si entra por E1 o E2
        if S==self.__E1__:
            ##Codigo a usar si el rayo entra por E1
            #Propagar hasta  H1
            C,P,D =self.complist["H1"]
            #Cambiando el rayo al systema de coordenadas de la superficie E1, que
            #en este caso en particular es el mismo de la componente C1
            R=ri.ch_coord_sys(P,D)
            #Calcular el punto de corte con H1 sin verificar apertura
            PI=self.__H1__._intersection(R)
            #Crear el rayo que va entre H1 y H2
            if self.complist["H1"][1][2] < self.complist["H2"][1][2]:
                R1=Ray(pos=PI, dir=(0,0,1), wavelength=wav, n=n)
            else:
                R1=Ray(pos=PI, dir=(0,0,-1), wavelength=wav, n=n)
            R2=R1.ch_coord_sys_inv(P,D)
            ri.add_child(R2)

            #Propagar hasta H2
            C,P,D =self.complist["H2"]
            R=R2.ch_coord_sys(P,D)
            PI=self.__H2__._intersection(R) #No se verifica apertura

            #### Calculate the reffraction in H2
            rx,ry,rz = ri.dir
            FP=ri.dir*self.f/abs(rz)
            d=FP-PI
            if self.f<0:
                    d=-d

            R3=Ray(pos=PI, dir=d, wavelength=wav, n=n)
            R4=R3.ch_coord_sys_inv(P,D)
            R2.add_child(R4)

            ##Propagar hasta C2
            C,P,D = self.complist["C2"]
            R=R4.ch_coord_sys(P,D)
            PI = self.__E2__.intersection(R)
            if self.__E2__.shape.hit(PI):
                R5=Ray(pos=PI, dir = R.dir, wavelength=wav, n=n)
                R6=R5.ch_coord_sys_inv(P,D)
                R4.add_child(R6)
            else:
                R4.intensity = 0.
        else:
            #codigo a usar si el rayo entra por E2
            #Propagar hasta  H2
            C,P,D =self.complist["H2"]
            #Cambiando el rayo al systema de coordenadas de la superficie E1, que
            #en este caso en particular es el mismo de la componente C1
            R=ri.ch_coord_sys(P,D)
            #Calcular el punto de corte con H1 sin verificar apertura
            PI=self.__H2__._intersection(R)
            #Crear el rayo que va entre H1 y H2
            if self.complist["H2"][1][2] < self.complist["H1"][1][2]:
                R1=Ray(pos=PI, dir=(0,0,1), wavelength=wav, n=n)
            else:
                R1=Ray(pos=PI, dir=(0,0,-1), wavelength=wav, n=n)
            R2=R1.ch_coord_sys_inv(P,D)
            ri.add_child(R2)

            #Propagar hasta H1
            C,P,D =self.complist["H1"]
            R=R2.ch_coord_sys(P,D)
            PI=self.__H2__._intersection(R) #No se verifica apertura

            #### Calculate the reffraction in H2
            rx,ry,rz = ri.dir
            FP=ri.dir*self.f/abs(rz)
            d=FP-PI

            if self.f<0:
                    d=-d

            R3=Ray(pos=PI, dir=d, wavelength=wav, n=n)
            R4=R3.ch_coord_sys_inv(P,D)
            R2.add_child(R4)

            ##Propagar hasta C1

            C,P,D = self.complist["C1"]
            R=R4.ch_coord_sys(P,D)
            PI = self.__E2__.intersection(R)
            if self.__E2__.shape.hit(PI):
                R5=Ray(pos=PI, dir = R.dir, wavelength=wav, n=n)
                R6=R5.ch_coord_sys_inv(P,D)
                R4.add_child(R6)
            else:
                R4.intensity = 0.

        return ri


    def distance(self,Ray ri):
        """Solo se tienen en cuenta C1 y C2, que son los limites de entrada y 
        salida del sistema.
        """
        #cdef np.ndarray P,D
        cdef list dist_list=[]
        cdef list pi_list=[]
        cdef list surf_list=[]
        #print self.complist


        for comp  in [self.complist["C1"],self.complist["C2"]]:
            C,P,D =comp
            #C,P,D = i
            #Reorientar el rayo, al sistema de coordenadas del elemento
            #y calcular el recorrido del rayo hasta chocar con la
            #el elemento
            R=ri.ch_coord_sys(P,D)

            Dist=C.distance(R)

            dist_list.append(Dist[0])

            pi_list.append(Dist[1])
            surf_list.append(Dist[2])

        mini=asarray(dist_list).argmin()

        return  dist_list[mini],pi_list[mini],surf_list[mini]