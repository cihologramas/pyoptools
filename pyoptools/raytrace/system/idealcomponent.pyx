from pyoptools.raytrace.component.component import Component
from pyoptools.raytrace.system.system import System
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.plane cimport Plane
from math import isinf
from numpy import asarray
cimport numpy as np


class IdealThickLens(System):
    """
    Class defined to create an Ideal ThickLens.
    
    Is is created as a subsystem because the propagation inside has multiple
    rays
    entrance-> principal plane1
    principal plane1 -> principal plane2
    principal plane 2 -> exit
    fl -> focal length
    
    **ARGUMENTS**    
    ===========================================================================
    shape          Shape of the lens (Entrance and exit surface's shape).
    thickness      Thinckness of the lens (Distance between the entrance and
                   exit surfaces)
    princ_planes   Tuple with the principal planes position. The position of
                   each principal plane is measured from its corresponding
                   entrance surface.
    pupils         (pupil_pos, pupil_diam,pupil_rs)
                    pupil_pos: pupil position measured from the active side
                    pupil_shape: pupil shape
                    pupil_rs: pupil reference surface. If true the reference 
                    surface will be E1. If false it will be E2 (see the source 
                    code).
                    If True it will be an entrance pupil for the rays entering
                    through surface E1, and an exit pupil for rays exiting
                    through E1, or an exit pupil for rays entering through E2,
                    and and an entrance pupil exiting through E2.
                    
                    If false, E1 and E2 are switched.  
                    None if no pupil is defined
                   
    complete_trace If set to false the trace between the principal rays will
                   not be shown. Still the trace from and to the entrance and
                   exit surfaces will be exact
    ==========================================================================
          
    The origin of the component is located in the middle point between the
    entrance and exit surface
    """

    def __init__(self, shape,thickness, princ_planes =(0.,0.), pupils=None,
                 f=100, complete_trace=False):

        System.__init__(self,n=1)

        self.__E1__ = Plane(shape=shape)
        self.__E2__ = Plane(shape=shape)
        self.__H1__ = Plane(shape=shape)
        self.__H2__ = Plane(shape=shape)

        if pupils is None:
            self.__PUP1__ = False
            self.__PUP2__ = False
        else:
            pup_pos = pupils[0]
            pup_shape = pupils[1]
            pup_rs = pupils[2]

            if pup_rs:
                self.__PUP1__ = True
                self.__PUP2__ = False
            else:
                self.__PUP1__ = False
                self.__PUP2__ = True

            if self.__PUP1__:
                self.__P1__ = Plane(shape=pup_shape)
                self.__C5__ = Component(surflist = [(self.__P1__,(0,0,0),(0,0,0))])
                self.complist["P1"] = (self.__C5__,(0,0,-thickness/2+pup_pos),(0,0,0))
                self.__PUP1__ = True

            if self.__PUP2__:
                self.__P2__ = Plane(shape=pup_shape)
                self.__C6__ = Component(surflist = [(self.__P2__,(0,0,0),(0,0,0))])
                self.complist["P2"] = (self.__C6__,(0,0,thickness/2+pup_pos),(0,0,0))

        self.__C1__ = Component(surflist = [(self.__E1__,(0,0,0),(0,0,0))])
        self.__C2__ = Component(surflist = [(self.__E2__,(0,0,0),(0,0,0))])
        self.__C3__ = Component(surflist = [(self.__H1__,(0,0,0),(0,0,0))])
        self.__C4__ = Component(surflist = [(self.__H2__,(0,0,0),(0,0,0))])

        self.complist["E1"]=(self.__C1__,(0,0,-thickness/2.),(0,0,0))
        self.complist["E2"]=(self.__C2__,(0,0,thickness/2.),(0,0,0))
        self.complist["H1"]=(self.__C3__,(0,0,-thickness/2.+princ_planes[0]),(0,0,0))
        self.complist["H2"]=(self.__C4__,(0,0,thickness/2.+princ_planes[1]),(0,0,0))

        self.f=f
        self.complete_trace = complete_trace

    def propagate_ray(self, Ray ri):

        pos = ri.pos
        dir = ri.dir
        wav = ri.wavelength
        label = ri.label
        parent = ri
        n=ri.n

        D,C,S = self.distance(ri)


        #Verificar si entra por E1 o E2
        if S==self.__E1__:
            P1 = self.complist["P1"] if self.__PUP1__ else False
            E1 = self.complist["E1"]
            H1 = self.complist["H1"]
            H2 = self.complist["H2"]
            E2 = self.complist["E2"]
            P2 = self.complist["P2"] if self.__PUP2__ else False
        else:
            P1 = self.complist["P2"] if self.__PUP2__ else False
            E1 = self.complist["E2"]
            H1 = self.complist["H2"]
            H2 = self.complist["H1"]
            E2 = self.complist["E1"]
            P2 = self.complist["P1"] if self.__PUP1__ else False

        ## Verificar si el rayo pasa por la pupila de entrada
        if P1:
            C,P,D = P1
            R=ri.ch_coord_sys(P,D)
            PI=C["S0"][0].intersection(R)
            if isinf(PI[0]):
                ST=True
                #Aunque el rayo le pega a la superficie de entrada, el rayo no
                # pasa por la pupila de entrada
            else:
                ST=False
        else: #No se verifican la pupila de entrada
            ST=False

        #Encontrar el punto de corte en la superficie de entrada E1
        C,P,D = E1
        R = ri.ch_coord_sys(P,D)
        PI = C["S0"][0].intersection(R)


        #Generar el rayo que va de E1 a H1
        if E1[1][2]<H1[1][2]:
            R = Ray(pos = PI, dir = R.dir,wavelength = wav, n=n,intensity = ri.intensity)
        else:
            R = Ray(pos = PI, dir = -R.dir,wavelength = wav, n=n,intensity = ri.intensity)

        R_E1 = R.ch_coord_sys_inv(P,D)

        if ST:
            #El rayo no pasa la pupila de entrada. Pintarlo hasta la superficie
            #de entrada
            R_E1.intensity=0
            ri.add_child(R_E1)
            return ri

        #Propagar hasta  H1
        C,P,D =H1
        R=R_E1.ch_coord_sys(P,D)
        #Calcular el punto de corte con H1 sin verificar apertura
        PI=C["S0"][0]._intersection(R)

        #Crear el rayo que va entre H1 y H2
        if H1[1][2] < H2[1][2]:
            R=Ray(pos=PI, dir=(0,0,1), wavelength=wav, n=n,intensity = ri.intensity)
        else:
            R=Ray(pos=PI, dir=(0,0,-1), wavelength=wav, n=n,intensity = ri.intensity)
        R_H1=R.ch_coord_sys_inv(P,D)

        #Propagar hasta H2
        C,P,D =H2
        R=R_H1.ch_coord_sys(P,D)
        PI=C["S0"][0]._intersection(R) #No se verifica apertura

        #### Calculate the refraction in H2
        rx,ry,rz = ri.dir
        FP=ri.dir*self.f/abs(rz)
        d=FP-PI
        if self.f<0:
            d=-d

        ##Crear el rayo que va entre H2 y E2
        if H2[1][2] < E2[1][2]:
            R=Ray(pos=PI, dir=d, wavelength=wav, n=n,intensity = ri.intensity)
        else:
            R=Ray(pos=PI, dir=-d, wavelength=wav, n=n,intensity = ri.intensity)

        R_H2=R.ch_coord_sys_inv(P,D)

        ##Propagar hasta E2
        C,P,D = E2
        R=R_H2.ch_coord_sys(P,D)
        PI = C["S0"][0].intersection(R)

        #Verificar la apertura de salida
        if isinf(PI[0]):
            PE2 = False
        else:
            PE2 = True

        ##Verificar la pupila de salida
        if P2:
            C0,P0,D0 = P2
            R=R_H2.ch_coord_sys(P,D)
            PII=C0["S0"][0].intersection(R)
            if isinf(PII[0])  or isinf(PII[0])  or isinf(PII[0]):
                ST2=False
            else:
                ST2=True
        else: # No se verifica la pupila de salida
            ST2=True

        if ST2 and PE2:
            ie2 = ri.intensity
        else:
            ie2=0

        if H2[1][2] < E2[1][2]:
            R=Ray(pos=PI, dir = R.dir, wavelength=wav, intensity = ie2, n=n)
        else:
            R=Ray(pos=PI, dir = -R.dir, wavelength=wav, intensity = ie2, n=n)
        R_E2=R.ch_coord_sys_inv(P,D)



        if self.complete_trace:
            ri.add_child(R_E1)
            R_E1.add_child(R_H1)
            R_H1.add_child(R_H2)
            R_H2.add_child(R_E2)
        else:
            R_E1_E2=Ray(pos=R_E1.pos, dir = R_E2.pos-R_E1.pos, wavelength=wav, n=n,intensity = ri.intensity)
            ri.add_child(R_E1_E2)
            R_E1_E2.add_child(R_E2)

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


        for comp  in [self.complist["E1"],self.complist["E2"]]:
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