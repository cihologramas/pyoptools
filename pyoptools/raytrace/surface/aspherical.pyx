#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# cython: profile=True
#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco <AUTHOR>
# All rights reserved.
# 
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  
# 
# 
# Author:         Ricardo Amézquita Orozco
# Description:    Aspherical surface definition module
#------------------------------------------------------------------------------
#
'''Module that defines support for Aspherical optical surfaces
'''


cdef extern from "math.h":
    double sqrt(double)
    
from numpy import  array, asarray, arange, polyadd, polymul, polysub, polyval,\
     dot, inf, roots, zeros, meshgrid, where, abs,sqrt as npsqrt

cimport numpy as np

from pyoptools.raytrace.surface.surface cimport Surface
from pyoptools.raytrace.ray.ray cimport Ray
from scipy.optimize import fsolve,ridder,newton,brentq,brenth,fminbound

#from ray_trace.surface.taylor_poly import eval_poly,  Poly_DyDx
from pyoptools.misc.Poly2D cimport *
from pyoptools.misc.definitions import inf_vect

cdef class Aspherical(Surface):
    """**Class that defines a high order aspherical surface**
    
    An aspherical surface is defined as::
    
      Z=(Ax*x**2+Ay*y**2)/(1+sqrt(1-(1+Kx)*Ax**2*x**2-(1+Ky)*Ay**2*y**2))+ poly2d()
    
    The poly2d is defined by a array in the same way as it is defined in the
    TaylorPoly Class
    
    Example 
        >>> cs=Aspherical(shape=Rectangle(size=(5,5)),Ax=.5,Ay=.3,Kx=.1, Ky=.1\ 
                            poly =poly2d((0,1,1)))
    """

    
    cdef public double Ax,Ay,Kx,Ky
    cdef public object poly
    cdef public double zmax,zmin
    
    def __init__(self,Ax=0.,Ay=0.,Kx=0.,Ky=0.,poly=None, *args, **kwargs):
        Surface.__init__(self, *args, **kwargs)
        self.Ax=Ax
        self.Ay=Ay
        self.Kx=Kx
        self.Ky=Ky
        self.poly=poly
        
        z=self.shape.mesh(ndat=(200,200),topo=self.topo)[2]
        cdef double zmax,zmin,dz
        zmax=z.max()
        zmin=z.min()
        dz=zmax-zmin
        
        # to make sure the limits of t are found correctly
        self.zmax=zmax+0.01*dz
        self.zmin=zmin-0.01*dz

        #if abs(self.zmin)<abs(self.zmax):
        #    self.zmin=0
        #else:
        #    self.zmax=0

        #Add attributes to the state list        
        self.addkey("Ax")
        self.addkey("Ay")
        self.addkey("Kx")
        self.addkey("Ky")
        self.addkey("poly")
        self.addkey("zmax")
        self.addkey("zmin")
        
    #~ def __reduce__(self):
        #~ 
        #~ args=(self.Ax, self.Ay, self.Kx, self.Ky, self.poly, self.reflectivity, self.shape)
        #~ return(type(self),args,self.__getstate__())
    
    
            
    cpdef topo(self, x, y):
        """**Returns the Z value for a given X and Y**
        
        This method returns the topography of the aspherical surface to be 
        used to plot the surface.
        """
        Ax=self.Ax
        Ay=self.Ay
        Kx=self.Kx
        Ky=self.Ky
        
        Z0=(Ax*x**2+Ay*y**2)/(1+npsqrt(1-(1+Kx)*Ax**2*x**2+(-(1+Ky))*Ay**2*y**2))
        if self.poly!=None:
            try:
                Z1=self.poly.meval(x,y)
            except (TypeError,ValueError):
                Z1=self.poly.eval(x,y)
        else:
            Z1=0.
        return Z0+Z1
        
    cpdef np.ndarray normal(self, int_p):
        """**Return the vector normal to the surface**
        
        This method returns the vector normal to the asphere at a point 
        ``int_p=(x,y,z)``.
        
        Note: It uses ``x`` and ``y`` to calculate the ``z`` value and the normal. 
        """
        cdef double Ax,Ay,Kx,Ky,x,y,z,dxA,dyA,dxP,dyP
        Ax=self.Ax
        Ay=self.Ay
        Kx=self.Kx
        Ky=self.Ky
        
        x, y, z= int_p
        
        dxA=(2*Ax*x)/(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1) +(Ax**2*(Kx+1)*x*(Ay*y**2+Ax*x**2))/(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)*(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)**2)
        dyA=(2*Ay*y)/(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)-(Ay**2*(-Ky-1)*y*(Ay*y**2+Ax*x**2))/(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)*(sqrt(Ay**2*(-Ky-1)*y**2-Ax**2*(Kx+1)*x**2+1)+1)**2)

        if self.poly!=None:
            Dx,Dy=self.poly.dxdy()
            dxP=Dx.peval(x,y)
            dyP=Dy.peval(x,y)
        else:
            dxP=0.
            dyP=0.
        
        N_=array((dxA+dxP, dyA+dyP, -1))
        
        return N_/sqrt(dot(N_, N_))
    
    cpdef double __f1(self,double t, Ray iray):
        cdef double Ax,Ay,Kx,Ky,Ox,Oy,Oz,Dx,Dy,Dz
        Ax=self.Ax
        Ay=self.Ay
        Kx=self.Kx
        Ky=self.Ky
        
        #Ox, Oy, Oz = iray.pos
        
               
        Ox=iray.cpos[0]
        Oy=iray.cpos[1]
        Oz=iray.cpos[2]
        
        #Dx, Dy, Dz = iray.dir
        
        Dx=iray._dir[0]
        Dy=iray._dir[1]
        Dz=iray._dir[2]
        
        cdef double X,Y,Z
        X=Dx*t+Ox
        Y=Dy*t+Oy
        Z=Dz*t+Oz
        return (Ay*Y**2+Ax*X**2)/(sqrt(Ay**2*(-Ky-1)*Y**2-Ax**2*(Kx+1)*X**2+1)+1)+self.poly.peval(X,Y) -Z
    
    cpdef double __f2(self,double t, iray):
        cdef double Ax,Ay,Kx,Ky,Ox,Oy,Oz,Dx,Dy,Dz
        
        Ax=self.Ax
        Ay=self.Ay
        Kx=self.Kx
        Ky=self.Ky
        
        #Ox, Oy, Oz = iray.pos
        Ox=iray.cpos[0]
        Oy=iray.cpos[1]
        Oz=iray.cpos[2]
        
        #Dx, Dy, Dz = iray.dir
        Dx=iray._dir[0]
        Dy=iray._dir[1]
        Dz=iray._dir[2]
        
        return (Ay*(Dy*t+Oy)**2+Ax*(Dx*t+Ox)**2)/(sqrt(Ay**2*(-Ky-1)*(Dy*t+Oy)**2-Ax**2*(Kx+1)*(Dx*t+Ox)**2+1)+1)-Dz*t-Oz
        
    
        
    cpdef _intersection(self,Ray iray):
        '''**Point of intersection between a ray and the asphere**

        This method returns the point of intersection  between the surface
        and the ray. This intersection point is calculated in the coordinate
        system of the surface.
        
        It uses an iterative process to calculate the intersection point.
        
        
           iray -- incident ray

        iray must be in the coordinate system of the surface
        '''
        
        #z=pz+t*dz t=(z-pz)/dz
        #Find the limits for t
        cdef double ta,tb,t,fa,fb,tm,tta,ttb,dt
        
        ta=(self.zmax-iray.pos[2])/iray.dir[2]
        tb=(self.zmin-iray.pos[2])/iray.dir[2]
        
        if self.poly!=None:
            fa=self.__f1(ta,iray)
            fb=self.__f1(tb,iray)
            if (fa<0 and fb>0) or (fa>0 and fb<0):
                t=brentq(self.__f1, ta,tb,(iray,),maxiter=1000)
            else: # there are more than 1 intersection points we are assuming 2
                #tm=fsolve(self.__f1, 0,(iray,),warning=False)
                #In new scipy version the warning kw is not supported
                tm=fsolve(self.__f1, 0,(iray,))
                
                if (tm<ta and tm<tb)or(tm>ta and tm>tb):
                    t=inf
                else:
                    dt=tb-ta
                    tta=tm-0.2*dt
                    ttb=tm+0.2*dt
                    t=brentq(self.__f1, tta,ttb,(iray,),maxiter=1000)
                
        else:
            fa=self.__f2(ta,iray)
            fb=self.__f2(tb,iray)
            
            if (fa<0 and fb>0) or (fa>0 and fb<0):
                t=brentq(self.__f2, ta,tb,(iray,),maxiter=1000)
            else: # there are more than 1 intersection points we are assuming 2
				#In new scipy version the warning kw is not supported
				#tm=fsolve(self.__f2, 0,(iray,),warning=False)
                tm=fsolve(self.__f2, 0,(iray,))
                
                if (tm<ta and tm<tb)or(tm>ta and tm>tb):
                    t=inf
                else:
                    dt=tb-ta
                    tta=tm-0.1*dt
                    ttb=tm+0.1*dt
                    t=brentq(self.__f2, tta,ttb,(iray,),maxiter=1000)
                
        ret_val= iray.pos+t*iray.dir
        
        return ret_val
               
    def _repr_(self):
        '''Return an string with the representation of an aspherical surface.
        '''
        return "Aspherical(shape="+str(self.shape)+",reflectivity="+\
                          str(self.reflectivity)+",Kx="+str(self.Kx)+\
                          ",Ky="+str(self.Ky)+",Ax="+str(self.Ax)+",Ay="+\
                          str(self.Ay)+",poly="+str(self.poly)+")"

