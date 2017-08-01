from numpy import asarray
cimport numpy as np
from pyoptools.raytrace.component.component import Component
from pyoptools.raytrace.system.system import System
from pyoptools.raytrace.ray.ray cimport Ray
from pyoptools.raytrace.surface.plane cimport Plane

class IdealThickLens(System):
    """
    Class defined to create an Ideal ThickLens.
    
    Is is created as a subsystem because the propagation inside has multiple rays
    entrance-> principal plane1
    principal plane1 -> principal plane2
    principal plane 2 -> exit
    fl -> focal lenght
    
    **ARGUMENTS**    
    ============== ===============================================================
    shape          Shape of the lens (Entrance and exit surface's shape).
    thickness      Thinckness of the lens (Distance between the entrance and exit
                   surfaces)
    princ_planes   Tuple with the principal planes position. The position of each 
                   principal plane is measured from its corresponding entrance 
                   surface. 
    complete_trace If set to false the trace between the principal rays will not
                   be shown. Still the trace from and to the entrace and exit 
                   surfaces will be exact 
    ============== ===============================================================
          
    The origin of the component is located in the midle point between the 
    entrance and exit surface
    """

    def __init__(self,shape,thickness, princ_planes =(0.,0.),f=100,complete_trace=False):

        System.__init__(self,n=1)

        self.__E1__ = Plane(shape=shape)
        self.__E2__ = Plane(shape=shape)
        self.__H1__ = Plane(shape=shape)
        self.__H2__ = Plane(shape=shape)

        self.__C1__ = Component(surflist = [(self.__E1__,(0,0,0),(0,0,0))])
        self.__C2__ = Component(surflist = [(self.__E2__,(0,0,0),(0,0,0))])
        self.__C3__ = Component(surflist = [(self.__H1__,(0,0,0),(0,0,0))])
        self.__C4__ = Component(surflist = [(self.__H2__,(0,0,0),(0,0,0))])

        self.complist["E1"]=(self.__C1__,(0,0,-thickness/2.),(0,0,0))
        self.complist["E2"]=(self.__C2__,(0,0,thickness/2.),(0,0,0))
        self.complist["H1"]=(self.__C3__,(0,0,-thickness/2.+princ_planes[0]),(0,0,0))
        self.complist["H2"]=(self.__C2__,(0,0,thickness/2.+princ_planes[1]),(0,0,0))

        self.f=f
        self.complete_trace = complete_trace

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
            #Encontrar el punto de corte en E1
            C,P,D = self.complist["E1"]
            R=ri.ch_coord_sys(P,D)
            PI=self.__E1__._intersection(R)

            #Generar el rayo que va de E1 a H1
            if self.complist["E1"][1][2]<self.complist["H1"][1][2]:
                R = Ray(pos = PI, dir = R.dir,wavelength = wav, n=n)
            else:
                R = Ray(pos = PI, dir = -R.dir,wavelength = wav, n=n)
            R_E1 = R.ch_coord_sys_inv(P,D)

            #Propagar hasta  H1
            C,P,D =self.complist["H1"]
            R=R_E1.ch_coord_sys(P,D)
            #Calcular el punto de corte con H1 sin verificar apertura
            PI=self.__H1__._intersection(R)
            #Crear el rayo que va entre H1 y H2
            if self.complist["H1"][1][2] < self.complist["H2"][1][2]:
                R=Ray(pos=PI, dir=(0,0,1), wavelength=wav, n=n)
            else:
                R=Ray(pos=PI, dir=(0,0,-1), wavelength=wav, n=n)
            R_H1=R.ch_coord_sys_inv(P,D)

            #Propagar hasta H2
            C,P,D =self.complist["H2"]
            R=R_H1.ch_coord_sys(P,D)
            PI=self.__H2__._intersection(R) #No se verifica apertura

            #### Calculate the reffraction in H2
            rx,ry,rz = R_H1.dir
            FP=R_H1.dir*self.f/abs(rz)
            d=FP-PI
            if self.f<0:
                    d=-d
            R=Ray(pos=PI, dir=d, wavelength=wav, n=n)
            R_H2=R.ch_coord_sys_inv(P,D)

            ##Propagar hasta E2
            C,P,D = self.complist["E2"]
            R=R_H2.ch_coord_sys(P,D)
            PI = self.__E2__.intersection(R)
            if self.__E2__.shape.hit(PI):
                R=Ray(pos=PI, dir = R.dir, wavelength=wav, n=n)
                R_E2=R.ch_coord_sys_inv(P,D)
            else:
                R_H2.intensity = 0.

            if self.complete_trace:
                ri.add_child(R_E1)
                R_E1.add_child(R_H1)
                R_H1.add_child(R_H2)
                if R_H2.intensity !=0:
                    R_H2.add_child(R_E2)
            else:
                R_E1_E2=Ray(pos=R_E1.pos, dir = R_E2.pos-R_E1.pos, wavelength=wav, n=n)
                ri.add_child(R_E1_E2)
                if R_H2.intensity !=0:
                    R_E1_E2.add_child(R_E2)
                else:
                    R_E1_E2.intensity = 0

        else:
            #codigo a usar si el rayo entra por E2
            # Es una copia del cvodigo en el if, cambiando el orden.
            # Hay que buscar una forma mas elegante de hacer esto

            #Encontrar el punto de corte en E2

            C,P,D = self.complist["E2"]
            R=ri.ch_coord_sys(P,D)
            PI=self.__E2__._intersection(R)

            #Generar el rayo que va de E2 a H2
            if self.complist["E2"][1][2]<self.complist["H2"][1][2]:
                R = Ray(pos = PI, dir = R.dir,wavelength = wav, n=n)
            else:
                R = Ray(pos = PI, dir = -R.dir,wavelength = wav, n=n)
            R_E2 = R.ch_coord_sys_inv(P,D)

            #Propagar hasta  H2
            C,P,D =self.complist["H2"]
            R=R_E2.ch_coord_sys(P,D)
            #Calcular el punto de corte con H2 sin verificar apertura
            PI=self.__H2__._intersection(R)
            #Crear el rayo que va entre H2 y H1
            if self.complist["H2"][1][2] < self.complist["H1"][1][2]:
                R=Ray(pos=PI, dir=(0,0,1), wavelength=wav, n=n)
            else:
                R=Ray(pos=PI, dir=(0,0,-1), wavelength=wav, n=n)
            R_H2=R.ch_coord_sys_inv(P,D)

            #Propagar hasta H1
            C,P,D =self.complist["H1"]
            R=R_H2.ch_coord_sys(P,D)
            PI=self.__H1__._intersection(R) #No se verifica apertura

            #### Calculate the reffraction in H1
            rx,ry,rz = R_H2.dir
            FP=R_H2.dir*self.f/abs(rz)
            d=FP-PI
            if self.f<0:
                    d=-d
            R=Ray(pos=PI, dir=d, wavelength=wav, n=n)
            R_H1=R.ch_coord_sys_inv(P,D)

            ##Propagar hasta E2
            C,P,D = self.complist["E1"]
            R=R_H1.ch_coord_sys(P,D)
            PI = self.__E1__.intersection(R)
            if self.__E1__.shape.hit(PI):
                R=Ray(pos=PI, dir = R.dir, wavelength=wav, n=n)
                R_E1=R.ch_coord_sys_inv(P,D)
            else:
                R_H2.intensity = 0.

            if self.complete_trace:
                ri.add_child(R_E2)
                R_E2.add_child(R_H2)
                R_H2.add_child(R_H1)
                if R_H1.intensity !=0:
                    R_H1.add_child(R_E1)
            else:
                R_E2_E1=Ray(pos=R_E2.pos, dir = R_E1.pos-R_E2.pos, wavelength=wav, n=n)
                ri.add_child(R_E2_E1)
                if R_H1.intensity !=0:
                    R_E2_E1.add_child(R_E1)
                else:
                    R_E2_E1.intensity = 0

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