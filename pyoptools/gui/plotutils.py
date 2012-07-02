#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Module that imports some pylab functions so it can be used in the system, 
without interfering with the wx threads in the wxRayTra.py
It can be used also in pylab
"""
from matplotlib import use
use('WXAgg') 
import matplotlib

matplotlib.interactive(True)

import matplotlib.pyplot as pl
from types import FunctionType
from wx import CallAfter
import time
import traceback
from pyoptools.misc.pmisc import wavelength2RGB

#Check if this is called from a standard python, ipython console, or from a wxipython console
#using the enviroment variable shell
try: 
    import os
    os.environ['shell']
    __WXRAYTRA=True
except KeyError:
    __WXRAYTRA=False

#Define the colormaps
cm=matplotlib.cm

class AsyncCall:
    #from threading import Event
    ''' Queues a func to run in thread of MainLoop.
    Code may wait() on self.complete for self.result to contain
    the result of func(*ar,**kwar).  It is set upon completion.
    Wait() does this.
    When no return value is needed, use wx.CallAfter'''
    def __init__( self, func, *ar, **kwar ):
        self.ready=False
        self.error=None
        self.func, self.ar, self.kwar= func, ar, kwar
        CallAfter( self.TimeToRun )
    def TimeToRun( self ):
        try:
            self.result=self.func( *self.ar, **self.kwar )
        except Exception as detail:
            # Print original error in console
            traceback.print_exc()
            # capture exception text. What can I do with this?
            self.txt=traceback.format_exc()

            self.error=detail
        self.ready=True
    def Wait( self, timeout= None, failval= None ):
        it=time.time()
        while not self.ready:
            time.sleep(0.1)
            if timeout!=None:
                if timeout>(time.time-it):
                    return failval
            if self.error!=None:
                print "Error redirecting the function ",self.func.__name__
                error=self.error
                self.error=None
                raise error
        return self.result  




## Redefine the plot functions if needed

if __WXRAYTRA:
    
    # Define all the pylab graphic functions
    a=dir(pl)
    for f in a:
        if (type(pl.__getattribute__(f)) is FunctionType) and f[0]!="_":
            func=" lambda *arguments, **keywords:AsyncCall(pl.%s,*arguments, **keywords).Wait()"%(f)
            ef=eval(func)
            ef.__doc__=pl.__getattribute__(f).__doc__
            globals()[f]=ef
    
    #Define other graphic functions
    import glplotframe2 as glp
    def glPlotFrame(os=None):
        pf= AsyncCall(glp.glPlotFrame, os)
        return pf.Wait()
else:
    from matplotlib.pyplot import *
    from glplotframe2 import glPlotFrame
    

## Definition of some help ploting functions. This must be moved to another module
from numpy import array
def spot_diagram(s):
    """Plot the spot diagram for the given surface, or element
    """
    hl=s.hit_list
    X=[]
    Y=[]
    COL=[]
    if len(hl) >0:
        for i in hl:
            p=i[0]
            # Hitlist[1] points to the incident ray
            col=wavelength2RGB(i[1].wavelength)
            X.append(p[0])
            Y.append(p[1])
            COL.append(col)
    max=array(X+Y).max
    min=array(X+Y).min
    plot(X,Y,"o",)
    axis("equal")
