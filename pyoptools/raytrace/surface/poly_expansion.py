 #!/usr/bin/env python
# -*- coding: UTF-8 -*-

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

The aspherical optical surfaces defined in this module are defined by the following equation

Z(s):=(c*s^2)/(1+sqrt(1-(1+k)*c^2*s^2))+TaylorPoly(x,y)

The taylor polynomial is defined as in the TaylorPoly Class. 
The aspherical optical surface, is modeled as a taylor polynomial.

This module is not working and should not be used
'''

from numpy import  array, asarray, arange, polyadd, polymul, polysub, polyval,\
     dot, inf, roots, zeros, meshgrid, sqrt,where, abs
#from enthought.traits.api import Tuple, Float, Array, Int
import sympy

from pyoptools.raytrace.surface.taylor_poly import TaylorPoly

class AsphericalHO(TaylorPoly):
    """Class that defines a high order aspherical surface
    
    The aspherical surface is defined as:
    Z=(Ax*x**2+Ay*y**2)/(1+sqrt(1-(1+Kx)*Ax**2*x**2-(1+Ky)*Ay**2*y**2))+TaylorPoly()

    The TaylorPoly is defined by a array in the same way as it is defined in the
    TaylorPoly Class
"""

    
    # Curvature of the surface
    Ax=Float(.1)
    
    Ay=Float(.1)
    # Conic constant
    Kx=Float(0)
    
    Ky=Float(0)
    
    # Aspherical coheficient definition
    ho_cohef=Array('d') #,shape=(3,),value=(0.,0.,0.))
    #  [[ x0y0, x1y0, x2y0,.....], 
    #  [[ x0y1, x1y1, x2y1,.....], 
    #  [[ x0y2, x1y2, x2y2,.....], 
    #  [[ x0y3, x1y3, x2y3,.....],
    #  [[ ... , ... , ... ,.....], 
    #  [[ x0y.., x1y..,x2y., ...], 
    
    def __init__(self, **traits):
        TaylorPoly.__init__(self, **traits)
        #Declare the analytical function
  
        Z=sympy.Function("Z")
        Ax,Ay,Kx,Ky=sympy.symbols(("Ax","Ay","Kx","Ky"))
        x, y =sympy.symbols('xy')
        Z=(Ax*x**2+Ay*y**2)/(1+sympy.sqrt(1-(1+Kx)*Ax**2*x**2-(1+Ky)*Ay**2*y**2));
        
        #Calculate taylor polynomial coheficients
        cohef=[[Z, ],]
        order=self.n
        for i in range(0, order+1, 2):
            if i!=0:
                cohef.append([sympy.diff(cohef[i/2-1][0], y, 2), ])
            for j in range(2, order-i+1, 2):
                cohef[i/2].append(sympy.diff(cohef[i/2][j/2 -1], x, 2))
        

        A_x=self.Ax
        A_y=self.Ay
        K_x=self.Kx
        K_y=self.Ky
        
        c=zeros((self.n+1, self.n+1))
        for i in range(0, order/2+1):
            for j in range(0,order/2- i+1):
                cohef[j][i]=cohef[j][i].subs(x, 0).subs(y, 0).subs(Ax, A_x).subs(Ay, A_y).subs(Kx, K_x).subs(Ky, K_y)/(sympy.factorial(2*i)*sympy.factorial(2*j))
                c[2*j, 2*i]=cohef[j][i].evalf()
        
        # Add the high order corrections
        if len(self.ho_cohef.shape)==2:
            cx, cy = c.shape 
            dx, dy =self.ho_cohef.shape
            mx=array((cx, dx)).max()
            my=array((cy, dy)).max()
            self.cohef=zeros((mx, my))
            self.cohef[0:cx, 0:cy]=c
            self.cohef[0:dy, 0:dy]=self.cohef[0:dy, 0:dy]+self.ho_cohef
        else:
            self.cohef=c
    def error(self, nsteps=100):
        """Method that returns the error in the polynomial expansion
        
        
        Returned Value
        ==============
        
        (PolyValue-ExactValue).std()
        """
        
        X, Y, Z= self.shape.mesh()
        
        Z0=self.eval_poly(X, Y)
        
        Ax,Ay,Kx,Ky=self.Ax, self.Ay, self.Kx, self.Ky
        
        Z1=(Ax*X**2+Ay*Y**2)/(1+sqrt(1-(1+Kx)*Ax**2*X**2-(1+Ky)*Ay**2*Y**2))
        
        return Z0*Z, Z1*Z
