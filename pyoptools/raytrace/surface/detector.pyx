#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# cython: profile=True
#------------------------------------------------------------------------------
# Copyright (c) 2007, Ricardo Amézquita Orozco
# All rights reserved.
#
# This software is provided without warranty under the terms of the GPLv3
# license included in LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.
#
#
# Author:          Ricardo Amézquita Orozco
# Description:     Detector surface definitión module
# Symbols Defined: CCD
#------------------------------------------------------------------------------
''' Module that defines optical detector surfaces
'''



from numpy import zeros, asarray, float64
#from enthought.traits.api import Bool, Property, List,  Tuple,  Float

from pyoptools.raytrace.surface.plane cimport Plane
from pyoptools.raytrace.ray.ray cimport Ray

from pyoptools.misc.pmisc import wavelength2RGB
from pyoptools.raytrace.shape.rectangular cimport Rectangular
#from gui.plotutils import plot, figure, cm
#from gui.plotutils import *

cdef class ArrayDetector(Plane):
    '''**CCD like detector surface.**

    Description the ArrayDetector class acts similar as a real CCD device. It holds
    a list of the coordinates of the rays that hits it.

     Example of a CCD:

        >>> cs=ArrayDetector(size=(10,10))
    '''



    # CCD phisical size
    #size=Tuple(Float(10), Float(10))
    cdef public tuple size

    def __init__(self,size=(10,10),*args, **kwargs):
        
        # Create a detector with a dummy size|
        Plane.__init__(self,shape=Rectangular(size=size), *args, **kwargs)
        self.size=size
        self.addkey("size")
        #Adjust the size
        #self.shape.size=self.size
        
        
        #~ #Add attributes to the state list        
        #~ self.state.append(self.size)
        #~ #print"Warning the transparent attribute in ArrayDetector is not working"
    
    #~ def __reduce__(self):
       #~ 
        #~ args=(self.size,  self.reflectivity, self.shape)
        #~ return(type(self),args,self.__getstate__())
    
     
    
    def get_histogram(self,size=(256,256)):
        """Method that returns the number of ray impacts per unit of area of the
        detector.
        
        Parameter:
    
        *size*
            Size of the detector in pixels. The phisical size of the detector is
            given at the surface creation.
        """
        px,py=size
        sx,sy=self.size
        dx,dy=sx/px,sy/py

        retval=zeros((px,py))

        # hit_list[0] holds the impact coordinates
        for i in self._hit_list:
            x,y,z=i[0]#z should always be 0
            nx=int(px*(x+sx/2.)/sx)
            ny=int(py*(y+sy/2.)/sy)
            retval[nx,ny]+=1
        return retval

    def get_color_histogram(self,size=(256,256)):
        """Method that returns the number of ray impacts per unit of area of the
        detector. It simulates a color image.
        
        Parameter:
    
        *size*
            Size of the detector in pixels. The phisical size of the detector is
            given at the surface creation.
        """
        px,py=size
        sx,sy=self.size
        dx,dy=sx/px,sy/py
        retval=zeros((px,py,3))
        for i in self._hit_list:
            x,y,z=i[0]
            # z must be always 0
            nx=int(px*(x+sx/2.)/sx)
            ny=int(py*(y+sy/2.)/sy)

            #TODO: The intensities are not taken into account. This should be fixed
            r,g,b=wavelength2RGB(i[1].wavelength)
            retval[nx,ny,0]+=r
            retval[nx,ny,1]+=g
            retval[nx,ny,2]+=b
        return retval

#    def spot_diagram(self,fig=None, title='Spot Diagram',style='o',color=False):
#        '''Plot a spot diagram in a pylab figure
        
#        Method that plots a spot diagram of the rays hitting the CCD.
        
#        *Attributes:*
        
#        *fig*
#            Pylab figure where the plot will be made. If set to None
#            a new figure will be created.
        
#        *style*
#            Symbol to be used to represent the spot. See the pylab plot 
#            documentation for more information.
        
#        *label*
#            String containing the label to show in the figure for this spot diagram.
#            Can be used to identify different spot diagrams on the same figure.
#        '''
        
#        if fig == None:
#            fig=figure()
#        X=[]
#        Y=[]
#        COL=[]
#        if len(self.__d_surf._hit_list) >0:
#            for i in self.__d_surf._hit_list:
#                p=i[0]
                # Hitlist[1] points to the incident ray
#                col=wavelength2RGB(i[1].wavelength)
#                X.append(p[0])
#                Y.append(p[1])
#                COL.append(col)
#        if label== None:
#            plot(X, Y, style,  figure=fig)
#        else:
#            plot(X, Y, style,label=label,figure=fig)
#            legend()
#        return fig
        
#    def im_show(self,size=(256,256),cmap=cm.gray,title='Image',color=False):
#        """Method to simulate the image detected in the CCD.
        
        
#        Parameters:
#        ___________
    
#        *title*
#            Title of the figure
            
#        *color*
#            When set to true, in simulates a color image.
            
#        """
        
#        if color==False:
#            data= self.get_histogram(size)
#            imshow(data, cmap=cmap)
#        else:
#            rdata=self.get_color_histogram(size)
#            im=imshow(rdata)
        #pf.Show()

#    def vtk_actor(self):
#        actor=PlaneSurf.vtkActor(self)
#        actor.property.color=(0,0,1)
#        return actor
