#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#cython: profile=True

#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amezquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amezquita Orozco
# Description:     Surface definition module
# Symbols Defined: surface
#------------------------------------------------------------------------------
#


cdef extern from "math.h":
    double sqrt(double) nogil
    double acos(double) nogil
    double cos(double) nogil
    bint isnan(double x) nogil 

cdef double pi= 3.1415926535897931

from warnings import warn
from sys import exit
from numpy import absolute, arccos, array, cos as npcos, dot, isinf as npisinf,inf, isnan as npisnan, \
                power, sqrt as npsqrt, sometrue,  zeros_like, ones_like,  histogram2d, \
                linspace, abs,  indices, argwhere, tan, polyfit,arange,where
from numpy.fft import fft2,  ifft2
from numpy.ma import array as ma_array

from scipy.signal import medfilt2d
from scipy import interpolate
from pylab import griddata, meshgrid

from inspect import getmembers, isroutine

from pyoptools.misc.pmisc import hitlist2int_list,hitlist2int, interpolate_g
from pyoptools.misc.cmisc.cmisc cimport norm_vect, empty_vec, vector_length
from pyoptools.misc.lsq import polyfit2d
from pyoptools.misc.Poly2D import ord2i
from pyoptools.raytrace.ray.ray cimport Ray, Rayf

from pyoptools.raytrace.shape.shape cimport Shape
from pyoptools.raytrace.shape.circular cimport Circular

from pyoptools.misc.definitions import inf_vect
from pyoptools.misc.picklable.picklable cimport Picklable


cimport numpy as np
np.import_array()
import cython

cdef extern from "math.h":
    bint isinf(double)

#Local declarations to improve speed. Need to check how can I make this global
cdef double infty=float("inf")


cdef class Surface(Picklable):
    '''
    :class:`Surface` is an abstract superclass for all optical surface objects.
    This class defines an API that all subclasses must support.
    
    All Surface subclasses share the attributes reflectivity, and shape,
    that are defined as follows:
    
    ============  =========================================================
    reflectivity  Float point number between 0 and 1. 0 indicates a 
                  completely transparent surface, 1 a completely reflective
                  surface. A value between are beam spliters.
    shape         Instance of the :class:`Shape`, that indicates the
                  surface's aperture shape.
    ============  =========================================================
    
    Each surface stores a hit_list, so the ray impact information can be 
    obtained. This list records the point of impact in the surface reference
    system, and a pointer to the hitting ray so intensity, wavelength, and other
    information about the ray can be retrieved. The reset method, clear this
    lists.
    
    '''
    
    #A reflectivity=0 indicates a transparent surface and a reflectivity =1 
    #'''
    
    #    indicates a perfect mirror. A value in between, indicates a beam splitter. The shape
    #    attribute is used to define the perimeter of the surface, and it is an
    #    instance of the Shape class.

    #~ '''

    def __init__(self,double reflectivity=0., shape=Circular(radius=10.)):
        self.reflectivity=reflectivity
        self.shape=shape
        self._hit_list=[]
        self.id =[] # The id of each surface gets registered when the component
                    #is created, and when the system is created
        Picklable.__init__(self,"reflectivity", "shape", "_hit_list","id")
        
            
        
    property hit_list:
        def __get__(self):
            # The list is converted to a tuple, so the internal list does not
            # get exposed, and it can not be modified by accident using an 
            # append
            return self._hit_list
        
        
    
    cpdef topo(self, x, y):
        '''Method that returns the topography of the surface
        
        The matrix returned is :math:`z=f(x,y)`. This method mus be overloaded in all 
        subclases of Surface.
        '''
        warn("Method topo, from class Surface, should be overloaded"+
             " in class "+self.__class__.__name__)

    cpdef int_nor(self,Ray iray):
        '''Point of intersection between a ray and a surface and the normal.
        
        Method that returns the point of intersection between a surface and
        a ray, and the vector normal to this point. I must return a tuple
        of the form ``((ix,iyi,z),(nx,ny,nz))``, where ``ix,iy,iz`` are the 
        coordinates of the intersection point, and nx,ny,nz are the components 
        of the vector normal to the surface.
        
        **Arguments:**
        
             iray -- incident ray in the coordinate system of the surface.
        
        This method uses surface.intersection, surface.normal, and 
        surface.shape.hit
        '''

        int_p= self.intersection(iray)
        N_= self.normal(int_p)
        
        # Check if the ray crosses the aperture
                    
        return int_p, N_


    cpdef np.ndarray  normal(self, int_p):
        '''Normal vector at the point ip.
        
        This method returns the normal vector at a specific intersection point, 
        given by ip. The normal vector must be normalized.
         
        This method must be overloaded in all Surface subclases. It contains
        the geometric specific code.
        
        **Arguments:**
        
        iray -- incident ray in the coordinate system of the surface.

        
        **Returns:**
        
        An array [dx,dy,dz] indicating the vector normal to the surface
        '''
            
        warn("Method normal from class Surface, should be overloaded"+
             " in class "+self.__class__.__name__)
        
    cpdef np.ndarray intersection(self,Ray iray):
        '''Point of intersection between a ray and a surface.
        
        This method returns the point of intersection  between the surface
        and the ray. This intersection point is calculated in the coordinate
        system of the surface.

        If there is not a point of intersection ie: ray is outside of the
        element aperture, it must return (numpy.inf,numpy.inf,numpy.inf)

        **Arguments:**
        
        iray -- incident ray
            
        iray must be in the coordinate system of the surface
        
        **Returns:**
        
        A vector [x,y,z] containing the coordinates of the point of 
        intersection between the ray and the surface. If there is no
        intersection returns [Inf, Inf, Inf].
        
        .. note::
            This method does not need to be overloaded. Overload 
            :meth:`surface._intersection` instead.
        '''
        
        int_p= self._intersection(iray)

        if self.shape.fhit(int_p[0],int_p[1],int_p[2])==False:
            return inf_vect

        return int_p

    cpdef _intersection(self,Ray iray):
        '''Point of intersection between a ray and a surface.

        This method returns the point of intersection  between the surface
        and the ray. This intersection point is calculated in surface coordinate
        system.
        
        This method should not check for the aperture  
     
        **Arguments:**
                
        ``iray`` -- incident ray

        iray must be in the coordinate system of the surface

        This function must be overloaded in all :class:`Surface` subclases. 
        It contains the geometric specific code. It must not  check for 
        the aperture.
        
        This function must not be called directly. You should call 
        :meth:`Surface.intersection` instead.
        '''
        warn("Method _intersection, from class Surface, should be overloaded"+
             " in class "+self.__class__.__name__)
             
    @cython.boundscheck(False)
    @cython.wraparound(False)
    
    cpdef distance(self,Ray iray):
        '''Distance propagated by a ray until intersection with a surface.
        
        This method returns the distance propagated by a ray since its origin
        and its point of intersection with a surface.
        
        **Args:**
        
        ``iray`` 
            Incident ray, iray must be in the coordinate system of the 
            surface
        
        **Returns:**
        
        A tuple ``(distance,point_of_intersection, surface)``, where:
        
        ``distance``
            The distance from the ray origin to the point of intersection 
            with the surface
        ``point_of_intersection``
            The point of intersection using the coordinate system of the 
            surface 
        ``surface``
            a pointer to the surface 
        
        '''
#~ 
        
        cdef np.ndarray[np.double_t, ndim=1] PI=self.intersection(iray)
        cdef double Dist
        # Dist is possitive if the current surface is ahead of the ray. 
        # If the surface is behind the ray, Dist becomes negative
        # If there is no intersection, Dist becomes inf
        if isinf(PI[0]) or isinf(PI[1]) or isinf(PI[2]) :#sometrue(isinf(PI)):
        #if PI[0]==inf or PI[0]==inf or PI[0]==inf:
            Dist=infty 
        else:
            Dist= dot(PI-iray.pos,iray.dir)
        
        if Dist>1e10:Dist=infty

        # Because of rounding errors, some times when the distance should be 0
        # It gives a very small number. A distance is too small, or 0 means 
        # that the ray was already refracted by the surface, and the surface
        # is already behind the ray.
        
        if Dist<1.e-10 : Dist= infty 

        return Dist,PI
        
    cpdef distance_s(self,Ray iray):
        '''Distance propagated by a ray until intersection with a surface.
        
        This method returns the distance propagated by a ray since its origin
        and its point of intersection with a surface. This method returns 
        positive and negative distances, and does not check for apertures.
        
        Arguments:
                
            iray -- incident ray

        iray must be in the coordinate system of the surface
        
        Return value:
        
            A tuple with the distance, the point of intersection using the
            coordinate system of the surface, and a pointer to the surface 
            (distance,point of intersection, surface)
        '''

        
        PI=self._intersection(iray)
        
        # Dist is possitive if the current surface is ahead of the ray. 
        # If the surface is behind the ray, Dist becomes negative
        # If there is no intersection, Dist becomes inf
        if isinf(PI[0]) or isinf(PI[1]) or isinf(PI[2]):
            Dist=inf
        else:
            Dist= dot(PI-iray.pos,iray.dir)
        
        if abs(Dist)>1e10:Dist=inf

        return Dist,PI# ,self
        
    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef propagate(self,Ray ri,double ni,double nr):
        '''Method to calculate the ray refracted on a surface.

        This method calculates the ray refracted (or rays refracted and 
        reflected)on a surface.
        
        Arguments:
        

             ri -- incident ray
             
             ni -- refraction index in the incident media n.
             
             nr -- refraction index in the refracted media

        ri must be in the coordinate system of the surface
        '''
        cdef double I,IA,gamma,gamma1,R
        cdef np.ndarray S1,PI,P,S2,A,A1,S3
        
        cdef np.float64_t *S1p  
        cdef np.float64_t *PIp  
        cdef np.float64_t *Pp  
        cdef np.float64_t *S2p 
        cdef np.float64_t *Ap  
        cdef np.float64_t *A1p 
        cdef np.float64_t *S3p 
        
        # Calculate the point of intersection of the ray with the surface, and
        # the normal vector to this point
        # The int_nor method should be overriden by each Surface subclasses

        #PI,P=self.int_nor(ri)
        PI= self.intersection(ri)
        P=  self.normal(PI)
        
        PIp = <np.float64_t *>np.PyArray_DATA(PI)
        Pp = <np.float64_t *>np.PyArray_DATA(P)
        
        
        # The information is added in the  system.propagate_ray method.
        # Add the information to the hitlist
        # self._hit_list.append((PI, ri))
        
        S1=ri._dir*ni #Use the numpy array instead the property. it is faster.
        S1p = <np.float64_t *>np.PyArray_DATA(S1)
        
        

        # Calculate the incident angle
        cdef double ddot=dot(S1,P)
        I=(acos(ddot/ni))
        
        
        #some times because rounding errors |ri.dir|>1. In this cases
        #I=nan
        if isnan(I): I=0.
        #IA=I
        # Take the correct normal
        #PA=P
        if I>pi/2:
            #P=-P
            Pp[0]=-Pp[0]
            Pp[1]=-Pp[1]
            Pp[2]=-Pp[2]
            ddot=dot(S1,P)
            I=(acos(ddot/ni))
        
        
        
       
        
        
        gamma= nr*sqrt((ni/nr*cos(I))**2 -(ni/nr)**2+1.)- ni*cos(I)
        
       
        #A=gamma*P
        A=empty_vec(3)
        Ap = <np.float64_t *>np.PyArray_DATA(A)
        Ap[0]=gamma*Pp[0]
        Ap[1]=gamma*Pp[1]
        Ap[2]=gamma*Pp[2]



        #S2=S1+A
        S2=empty_vec(3)
        S2p = <np.float64_t *>np.PyArray_DATA(S2)
        S2p[0]=S1p[0]+Ap[0]
        S2p[1]=S1p[1]+Ap[1]
        S2p[2]=S1p[2]+Ap[2]
        
        #S2=nr*(norm_vect(S2))
        norm_vect(S2)
        S2p[0]=nr*S2p[0]
        S2p[1]=nr*S2p[1]
        S2p[2]=nr*S2p[2]
       
        #R=acos(dot(S2/nr,P))

        # If the refraction index is negative, the surface should behave
        # as a mirror

        # If reflectivity==0 the optical surface is not a BeamSplitter

        # If there is total internal reflection, the surface should not 
        # behave as a BeamSplitter

        # If the optical system is a BeamSplitter this method should return
        # a list [Transmited ray, Reflected ray]
        
        #Note the fast ray creation function do not normalize the dir of the ray
        norm_vect(S2) #S2/sqrt(dot(S2,S2))
        
        if(ni/nr)<0:
            # This case should never happen
            # For a mirror use reflectivity=1

            warn("For a mirror use reflectivity=1, not n<0")
            return [Ray(pos=PI,dir=-S2,intensity=ri.intensity,
                        wavelength=ri.wavelength,n=absolute(ni),
                        label=ri.label, orig_surf=self.id)]
            

        elif (self.reflectivity ==0) and not(sometrue(npisnan(S2))):
            # Normal refraction case
            #return [Ray(pos=PI,dir=S2,intensity=ri.intensity,
            #           wavelength=ri.wavelength,n=nr,
            #           label=ri.label, orig_surf=self.id)]
            return [Rayf(PI,S2,ri.intensity, ri.wavelength, nr, ri.label, None ,0,self.id,0)]
        elif sometrue(npisnan(S2)):
            # Total internal refraction case
            gamma1= -2.*ni*cos(I)
        
            #A1=gamma1*P
            A1=empty_vec(3)
            A1p = <np.float64_t *>np.PyArray_DATA(A1)
            A1p[0]=gamma1*Pp[0]
            A1p[1]=gamma1*Pp[1]
            A1p[2]=gamma1*Pp[2]
            
            #S3=S1+A1
            S3=empty_vec(3)
            S3p = <np.float64_t *>np.PyArray_DATA(S3)
            S3p[0]=S1p[0]+A1p[0]
            S3p[1]=S1p[1]+A1p[1]
            S3p[2]=S1p[2]+A1p[2]

        
            norm_vect(S3) #S3=S3/sqrt(dot(S3,S3))
        
        
            #return [Ray(pos=PI,dir=S3,
            #            intensity=ri.intensity,
            #            wavelength=ri.wavelength,n=ni,
            #            label=ri.label, orig_surf=self.id)]
            return [Rayf(PI,S3,ri.intensity, ri.wavelength, ni, ri.label, None ,0,self.id,0)]
             
        else:
            # BeamSplitter case
            gamma1= -2.*ni*cos(I)
            #A1=gamma1*P
            A1=empty_vec(3)
            A1p = <np.float64_t *>np.PyArray_DATA(A1)
            A1p[0]=gamma1*Pp[0]
            A1p[1]=gamma1*Pp[1]
            A1p[2]=gamma1*Pp[2]
            
            #S3=S1+A1
            S3=empty_vec(3)
            S3p = <np.float64_t *>np.PyArray_DATA(S3)
            S3p[0]=S1p[0]+A1p[0]
            S3p[1]=S1p[1]+A1p[1]
            S3p[2]=S1p[2]+A1p[2]

            norm_vect(S3) #S3=S3/sqrt(dot(S3,S3))
            if self.reflectivity!=1.:
                return [#Ray(pos=PI,dir=S2,
                        #    intensity=ri.intensity*(1.-self.reflectivity),
                        #    wavelength=ri.wavelength,n=nr, label=ri.label, orig_surf=self.id),
                        Rayf(PI,S2,ri.intensity*(1.-self.reflectivity), ri.wavelength, nr, ri.label, None ,0,self.id,0),
                        #Ray(pos=PI,dir=S3,
                        #    intensity=ri.intensity*self.reflectivity,
                        #    wavelength=ri.wavelength,n=ni,label=ri.label, orig_surf=self.id)
                        Rayf(PI,S3,ri.intensity*self.reflectivity, ri.wavelength, ni, ri.label, None ,0,self.id,0)
                        ]
            else:
                #return [Ray(pos=PI,dir=S3,
                #            intensity=ri.intensity*self.reflectivity,
                #            wavelength=ri.wavelength,n=ni,label=ri.label, orig_surf=self.id)]
                return [Rayf(PI,S3,ri.intensity, ri.wavelength, ni, ri.label, None ,0,self.id,0)]

    cpdef pw_propagate1(self, Ray ri,ni,nr, rsamples, isamples, knots):
        '''Method to calculate wavefront emerging from the surface
        
        .. Note::
            This method needs to be checked, because the incoming plane 
            wave has not a constant intensity. The ray distribution is 
            affected by the surface topography. it is better to use the 
            slow version pw_propagate.

        This method calculates the wavefront emerging from the optical surface
        when illuminated by an unity amplitude plane wave. 
        The k vector of the incoming PW is given by the direction of ri. 
        The wavelength of the PW is given by the wavelength of ri.
        The origin of ri is not used at all.
       
        The input and output planes are located at z=0 in the surface coordinated 
        system. 
        
        Note: this method needs to be checked, because the incoming plane wave
        has not a constant intensity. The ray distribution is affected by the 
        surface topography. it is better to use the slow version pw_propagate1.
        
        Note: The ray comes from the negative side. Need to change this
        
        Arguments:
        

             ri -- isurface.pyncident ray
             
             ni -- refraction index in the incident media n.
             
             nr -- refraction index in the refracted media

        ri must be in the coordinate system of the surface
        '''
         
        # Calculate the position maximum and minimum z values for the surface
        X, Y, Z=self.shape.mesh(topo=self.topo)
        
        #Calculate a mesh where the test rays will be shot from
        X, Y, H=self.shape.mesh(ndat=rsamples)#,  size=(xmin, xmax, ymin, ymax))
        X, Y, Z=self.shape.mesh(ndat=rsamples,  topo=self.topo)

        
        
        # The ecuation where the rays are shooted to is Z=0
        N_=array((0, 0, 1))
        
        # The ecuation where the rays are shooted from is
        # X*dx+Y*dy+Z*dz=0, where dx,dy,dx are the plane normal.
        N_p=ri.dir
        
        ix, iy=indices(X.shape)
        
        # Create a list of indexes to be able to select the points that are going 
        # to be used as spline generators, and as control points
        idx=where((ix+iy) %2 ==0, False, True)
        
        X=X.flatten()
        Y=Y.flatten()
        Z=Z.flatten()
        H=H.flatten()
        idx=idx.flatten()
        # Get the rays that hit the surface
        li=H.nonzero()[0]
        #rin=Ray( dir=L1, wavelength=wavelength)
        dir=ri.dir
        S1=dir*ni
        #Intersection point and optical path for the incident rays
        xi=[]
        yi=[]
        zi=[]
        ii=[]
        
        for i in li:
            #Destination of the ray
            P1=array((X[i], Y[i], Z[i]))
            
            #Normal to the destination point
            Np=self.normal(P1)
            
            # Calculate incident and refracted angle. To optimize the routine, 
            # only normal diffracted ray is taken into account
            
            # Calculate the incident angle
            I=(arccos(dot(S1,Np)/ni))
            
            # Take the correct normal
            if I>pi/2:
                Np=-Np
                I=(arccos(dot(S1,Np)/ni))
            # Calculate the diffracted direction
            
            gamma= nr*sqrt(power(ni/nr*cos(I),2)-power(ni/nr,2)+1.)- ni*cos(I)
            A=gamma*Np
            S2=nr*(S1+A)
            L2=(S2/sqrt(dot(S2,S2)))
            ###  
            
            # Calculate the distance between P1 ant the plane X*dx+Y*dy+Z*dz=0, 
            # along the N_p direction
            din=-dot(N_p,-P1)
            
            # Calculate the distance between P1 Z=0 along the L2 line
            ddi=-P1[2]/L2[2]
            
            #Calculate the intersection point
            PIR=P1+ddi*L2

            
            # Calculate optical path
            d=din*ni+ddi*nr
            if d!=inf:
                x, y, z=PIR
                xi.append(x)
                yi.append(y)
                zi.append(d)
                ii.append(idx[i])    
        xi=array(xi)
        yi=array(yi)
        zi=array(zi)
        ii=array(ii)
        
        xmax=xi.max()#=abs(array(xi+yi)).max()
        ymax=xmax
        xmin=xi.min()#-ymax
        ymin=xmin#-xmax
        
        xx=linspace(xmin, xmax, isamples[0])
        yy=linspace(ymin, ymax, isamples[1])
        
        # Use only half of the samples to create the Spline, 
        isp=argwhere(ii==True)
        ich=argwhere(ii==False)
        
        xsp=xi[isp]
        ysp=yi[isp]
        zsp=zi[isp]
        
        xch=xi[ich]
        ych=yi[ich]
        zch=zi[ich]
        
        #Distribute homogeneously the knots
        xk=linspace(xsp.min(), xsp.max(),knots)
        yk=linspace(ysp.min(), ysp.max(),knots)
        
        # LSQBivariateSpline using some knots gives smaller error than 
        # SmoothBivariateSpline
        di=interpolate.LSQBivariateSpline(xsp, ysp, zsp, xk[1:-1],  yk[1:-1])
        #di=interpolate.SmoothBivariateSpline(xsp, ysp, zsp)
        
        # Evaluate error
        zch1=di.ev(xch, ych)
        er=(zch.flatten()-zch1).std()
        
        #I, xe, ye=histogram2d(yi, xi, (xx, yy))
        I=hitlist2int(xi, yi, xi,  xx, yy)
        d=ma_array(di(xx,yy).transpose(), mask=I.mask)
        
        #XD, YD=meshgrid(xx, yy)
        
        #d1=griddata(xi,  yi,  zi,  xx, yy)
        return I, d, er
        
        
    
    cpdef pw_propagate(self, Ray ri,ni,nr, rsamples, shape, order,z):
        '''Method to calculate wavefront emerging from the surface


        This method calculates the wavefront emerging from the optical surface
        when illuminated by an unity amplitude plane wave. 
        The k vector of the incoming PW is given by the direction of ri. 
        The wavelength of the PW is given by the wavelength of ri.
        The origin of ri is not used at all.
               
        Arguments:
        

             ri -- surface incident ray
             
             ni -- refraction index in the incident media n.
             
             nr -- refraction index in the refracted media
             
             rsamples -- number of rays used to sample the surface
             
             shape -- shape of the output wavefront
             
             order -- order of the polynomial fit
             
             z -- Z position of the input and output plane. The origin is the surface vertex

        ri must be in the coordinate system of the surface
        Note: The ray comes from the negative side. Need to change this
        '''
        from ray_trace.surface import Plane
        
        xi,yi,zi=self.pw_propagate_list(ri,ni,nr,rsamples,z)
    
        xmin,xmax,ymin,ymax=self.shape.limits()
        xx=linspace(xmin, xmax, shape[0])
        yy=linspace(ymin, ymax, shape[1])
        
        #I, xe, ye=histogram2d(yi, xi, (xx, yy))
        I=hitlist2int(xi, yi, xi,  xx, yy)

        
        p,er=polyfit2d(xi, yi, zi,order=order)
        
        d=ma_array(p.evalvv(xx,yy),mask=I.mask)
        #d,er=interpolate_g(xi,yi,zi,xx,yy,knots=knots, error=True,mask=I.mask)
        
        #XD, YD=meshgrid(xx, yy)

        #d1=griddata(xi,  yi,  zi,  xx, yy)
        return I, d, er
        
    cpdef pw_propagate_list(self, Ray ri,ni,nr, rsamples,z):
        '''Method to calculate wavefront emerging from the surface

        This method calculates samples of the  wavefront emerging from 
        the optical surface when illuminated by an unity amplitude plane 
        wave. 
        The k vector of the incoming PW is given by the direction of ri. 
        The wavelength of the PW is given by the wavelength of ri.
        The origin of ri is not used at all.
        
        The returned value is a list containing the x,y coordinates of the ray list
        in the output surface, and the optical path at such point.
        
        To create an output matrix, this values must be interpolated.
        
        
               
        Arguments:
        

             ri -- incident ray
             
             ni -- refraction index in the incident media n.
             
             nr -- refraction index in the refracted media
             
             rsamples -- number of rays used to sample the surface (Tuple)
             
             z -- Z position of the input and output plane. The origin is the surface vertex

        ri must be in the coordinate system of the surface
        Note: The ray comes from the negative side. Need to change this
        '''
        from pyoptools.raytrace.surface import Plane
        
        #plane where the aperture is located 
        Za=z
        
        # Create an array of rays to simulate the plane wave
        dir=ri.dir
        xmin,xmax,ymin,ymax=self.shape.limits()
        
        # Calculate the position maximum and minimum z values for the surface
        #X, Y, H=self.shape.mesh(ndat=rsamples)
        #X, Y, H=self.shape.mesh(ndat=rsamples)

        
        #X;Y;Z coordinates of the points the rays are aiming to  in the 
        # aperture plane. The 20% increase in size, is to assure that
        # all the surface is sampled even if the rays are tilted
        dx=float(xmax-xmin)*.1
        dy=float(ymax-ymin)*.1
        
        Xa=linspace(xmin-dx, xmax+dx,int(rsamples[0]*1.2))
        Ya=linspace(ymin-dy, ymax+dy,int(rsamples[1]*1.2))
        #print Xa,Ya
        Xa,Ya=meshgrid(Xa,Ya)
      
        #
        Xa=Xa.flatten()
        Ya=Ya.flatten()  
        
        # The ecuation where the rays are shooted from is
        # X*dx+Y*dy+Z*dz=0, where dx,dy,dx are the plane normal.
        Z=-(Xa*dir[0]+Ya*dir[1])/dir[2]+Za
        #H=H.flatten()

        #li=H.nonzero()[0]
        pl=Plane()
        xi=[]
        yi=[]
        zi=[]
        for i in range(len(Xa)):
                     
            P0=array((Xa[i],Ya[i],Z[i]))
            
            ri=Ray(pos=P0,dir=dir,wavelength=ri.wavelength)
            #if the ray do not intersect the surface, break current iteration
            if sometrue(npisinf(self.intersection(ri))): continue
            
            #print self.intersection(ri)
            rd=self.propagate(ri,ni,nr)
          
            #take only the transmited ray
            rd=rd[0]
            
            #Translate rd, to put it in the apperture reference plane
            rd.pos[2]=rd.pos[2]-Za
            di=self.distance_s(ri)[0]
            dr=pl.distance_s(rd)[0]
            
            PI=pl._intersection(rd)
            
            # Calculate optical path
            d=di*ni+dr*nr
            #print d,di,dr,PI
            if d!=inf:
                x, y, z=PI
                xi.append(x)
                yi.append(y)
                zi.append(d)
        xi=array(xi)
        yi=array(yi)
        zi=array(zi)
        
        return xi,yi,zi
        
            
        
    cpdef wf_propagate(self, wf,ni,nr, samples,shape, knots):
        '''Method to calculate wavefront emerging from the surface


        This method calculates the wavefront emerging from the optical surface
        when illuminated by an arbitrary wavefront.
              
        The input and output planes are located at z=0 in the surface coordinated 
        system. 
        
        Arguments:
        
             wf -- Field instance containing the incoming wavefront
             
             ni -- refraction index in the incident media n.
             
             nr -- refraction index in the refracted media.
             
             samples -- Tuple containing the number of rays used to sample the field.
             
             shape -- Tuple containing the shape of the output field
        '''
        from ray_trace.surface import Plane
        # Get the wavefront ray representation
        # TODO: This representation only takes into account the phase, but not the intensity.
        #       This has to be fixed, because in practice this is not OK
        rays=wf.rayrep(samples[0],samples[1])
        
        #TODO: check which rays pass inside the apperture
        
        #rin=Ray( dir=L1, wavelength=wavelength)
        
        #Intersection point and optical path for the incident rays

        xi=[]
        yi=[]
        zi=[]
    
        pl=Plane()
        for ri in rays:
            #Calculate the intersection point
            rd=self.propagate(ri,ni,nr)
            #Take only de transmited ray
            rd=rd[0]
            
            # Incident ray propagation distance until the optical surface
            di=self.distance(ri)[0]
            
            
            
            # Refracted ray propagation until the output surface (plane Z=0)
            # The distance methos is not used, because it eliminates the negative
            # propagation values
            
            PI=pl._intersection(rd)
            dr=dot(PI-rd.pos,rd.dir)

            d=di*ni+dr*nr+ri.optical_path_parent()

            if d!=inf:
                x, y, z=PI
                xi.append(x)
                yi.append(y)
                zi.append(d)
            else:
                print ri
                
        d=interpolate_g(xi,yi,zi,samples=shape,knots=knots, error=False,mask=None)
        
        return d
        
    def __repr__(self):
        #~ '''Return an string with the representation of the optical surface
#~ 
        #~ It must be overloaded in all subclasses
        #~ '''

        return "OptSurf(reflectivity="+str(self.reflectivity)+")"

    

    
    cpdef reset(self):
        ''' Remove information from the hit_list
        '''
        self._hit_list=[]
        
        
    cpdef pw_cohef(self,ni,nr,ilimit, slimit, step, order, rsamples,zb):
        '''Method to generate the taylor polinomial coheficients that can be 
        used to obtain the phase and intensity for diferent plane wave 
        inclinations.
     
        Notes: 
               - The pupil is normalized to the radius of the lens
               - This method asumes a circular shaped pupil
               - The phase is returned as an optical path
               - The intensity is normalized to 1
               
        Arguments:
        
        =========  ======================================================        
        ni,nr      Refraction index from the incident and refracted sides
        ilimit     Inferior limit for incidence angle of the plane wave 
                   in radians
        slimit     Superior limit for incidence angle of the plane wave 
                   in radians
        step       Step to be used to generate the interpolation data
        order      Order of the Taylor interpolation used
        rsamples   Tuple containing the number of ray samples to be used 
                   in each direction
        zb         Z position of the plane where the measurementas are 
                   made. The origin is the vertex of the surface. If 
                   None, an estimate position is taken.
        =========  ======================================================  
        
        '''
        xm=self.shape.limits()[1]
        #print xm
        #Get the z position of the border of the surface
        
        if zb == None:
            zm=self.topo(xm,0)
        else: zm=zb
        
        # Optical path coheficient list 
        opcl=[]
        
        # Intensity coheficient list
        icl=[]
        
        xd=[]
        for i in arange(ilimit,slimit,step):
            xd.append(i)
            
            x,y,d =self.pw_propagate_list(Ray(dir=(tan(i), 0, 1)),ni,nr, rsamples=rsamples,z=zm)
            
            #Normalize the pupil
            x=x/xm; y=y/xm
            
            # Get the optical path polynomial coheficients
            pf,ef=polyfit2d(x, y, d,order=order)
            
            #get the intensity data
            xi,yi,I= hitlist2int_list(x, y)
            
            # Get the intensity polynomial coheficients
            pi,ei=polyfit2d(xi, yi, I,order=order)
            
            #TODO: Print something if error is too big 
            #if ei>0.001: print ""
            #if ef>1e-6: print ""
            
            opcl.append(pf.cohef.flatten())

            icl.append(pi.cohef.flatten())    

        dph=array(opcl)
        di=array(icl)

        # Lista de los coheficientes de los polinomios para generar los coheficientes
        phcohef=[] 
        
        #get the number of coheficients the polinomial expansion has
        ncohef=ord2i(order)
        
        for i in range(ncohef):
            #figure()
            #plot(xd,d[:,i])
            pof=polyfit(xd,dph[:,i],15)
            phcohef.append(pof)
        
        icohef=[]
        
        for i in range(ncohef):
            #figure()
            #plot(xd,d[:,i])
            pof=polyfit(xd,di[:,i],15)
            icohef.append(pof)
        
        return phcohef, icohef, zm 
        
    
    def polylist(self):
        
        points=[]
        
        X,Y=self.shape.pointlist()
        Z=self.topo(X,Y)
        
        for i in range(len(X)):
            points.append((X[i],Y[i],Z[i]))
        
        from matplotlib.delaunay import delaunay
        
        #Need to find a beter way to do this not using delaunay# or maybe to generate all using triangulations????
        cs,e,trip,trin=delaunay(X,Y)
        return points, trip
        
