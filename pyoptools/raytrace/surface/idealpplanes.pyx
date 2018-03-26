#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Modulo con clases y funciones auxiliares.
"""


from numpy import array, dot, inf, float64, zeros, asarray,isnan, mean, sqrt
from pyoptools.misc.definitions import inf_vect

from pyoptools.raytrace.surface.surface import Surface
from pyoptools.raytrace.ray.ray import Ray


class IdealPPlanes(Surface):
    """Clase que representa un par de superficies principales ideales. Se utiliza para crear
    lentes ideales gruesas
    """
    def __init__(self,f=100,d=20,*args, **kwargs):
        """
        f: Focal length
        d: Distance between planes
        """
        Surface.__init__(self,*args, **kwargs)
        #self.phm=phm
        #dxdy=phm.dxdy()
        #self.phx=<poly2d>dxdy[0]
        #self.phy=<poly2d>dxdy[1]
        #self.M=M

        self.f=f
        self.d=d

        #Add attributes to the state list
        #self.addkey("f")


    def topo(self, x, y):
        return zeros(asarray(x).shape)

    def _intersection(self,A):
        """Returns the intersection point between a ray and an the XY plane

        """
        #N_=array([0.,0.,1.])

        P1=A.pos     # Punto que pertenece al rayo "Origen" del rayo
        L1= A.dir #Vector paralelo a la linea


        #if dot(N_,L1) ==0 : return inf_vect
        if L1[2] ==0 : return inf_vect

        #print N_,P1,L1
        #print dot(N_,-P1),dot(N_,L1)
        #u=dot(N_,-P1)/dot(N_,L1)
        u1=-(P1[2]-self.d/2.)/L1[2]
        u2=-(P1[2]+self.d/2.)/L1[2]
        if abs(u1)<abs(u2):
            if u1!=0:
                u=u1
            else:
                u=u2
        else:
            if u2!=0:
                u=u2
            else:
                u=u1
        # Si u es muy grande, no hay intersección

        retval=P1+u*L1
        return retval


    def normal(self, ri):
        """Method that returns the normal to the surface
        """
        N_=array((0.,0.,1.)).astype(float64)
        return (N_)

    def propagate(self,ri,ni,nr):

        PI,P=self.int_nor(ri)
        #l=ri.wavelength*1e-3 # Express wavelength in millimeters so all the method works in millimeters
        rx,ry,rz=ri.dir
        #Get the focussing point as the point where the principal ray hits the focal plane

        FP=ri.dir*self.f/abs(rz)
        #Las ecuaciones de refraccion usadas funcionan para el caso donde el plano
        #esta en Z=0. Como aca el plano no esta en z=0, hay que moverlo para el
        #calculo, y luego moverlo nuevamente a las posiciones adecuadas
        if PI[2]<0:
            PI[2]=0
            d=FP-PI
            PI[2]=self.d/2.
        else:
            PI[2]=0
            d=FP-PI
            PI[2]=-self.d/2.

        ret=[]
        if self.reflectivity!=1:
                ret.append(Ray(pos=PI,dir=d,
                                intensity=ri.intensity,
                                wavelength=ri.wavelength,n=ni,label=ri.label, orig_surf=self.id))
        if self.reflectivity!=0:
                #print "not 0"
                ret.append(Ray(pos=PI,dir=PI-FP,
                                intensity=ri.intensity,
                                wavelength=ri.wavelength,n=ni,label=ri.label, orig_surf=self.id))

        #print self.reflectivity
        return ret