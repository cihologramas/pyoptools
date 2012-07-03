#!/usr/bin/env python
# -*- coding: UTF-8 -*-

#from enthought.tvtk.api import tvtk
#from enthought.tvtk.tools import ivtk
#from mayavi import ivtk
#from gui.plotutils import AsyncCall

#from enthought.traits.api import Instance

#from ray_trace.system import System
#from wx import CallAfter

#class _PlotFrame(ivtk.IVTK):
#    opsys=Instance(System)    
#    def __init__(self, **traits):
#        ivtk.IVTK.__init__(self,**traits)
#        #self.open()
#        #self.plot()
#        
#    def plot(self):
#        AsyncCall(self.open).Wait()
#        actor=self.opsys.tvtk_actor()
#        AsyncCall(self.scene.add_actor, actor).Wait()
#        self.reset()
#    def clear(self):
#        AsyncCall(self.open).Wait()
#        AsyncCall(self.scene.renderer.remove_all_view_props).Wait()
#        self.reset()#

#    def reset(self):
#        AsyncCall(self.scene.reset_zoom)

#    def close(self):
#        AsyncCall(self._close)
#        #AsyncCall(ivtk.IVTK.close, self)
#    def _close(self):
#        self.scene.renderer.remove_all_view_props()
#        ivtk.IVTK.close(self)
        
#def PlotFrame(**traits):
#    OB= AsyncCall(_PlotFrame, **traits).Wait()
#    OB.plot()
#    return OB

#_pf= Instance(tvtk_frame,()) #this initializes the tvtk_frame

#    def __init__(self,opsys):
        #HasTraits.__init__(self,**traits)
#        self._pf=AsyncCall(tvtk_frame, opsys).Wait()
#        AsyncCall(self._pf.open).Wait()
        #self._pf.scene.background=(0,0,0)
        
#        if self.opsys!=None:
#            print self.opsys
            #self.actor=self.opsys.tvtk_actor()
#            self._pf.plot(self.opsys)
        #self._pf.reset()
        
#    def get_scene(self):
#        """ Returns the tvtk scene
#        """
#        return self._pf.scene

#    def add_text(self,text="",scale=(1.,1.,1.), position=(0.,0.,0.),follower=True):
        #atext = tvtk.TextSource()
#        atext = tvtk.VectorText()
#        atext.text=(text)
#        textMapper = tvtk.PolyDataMapper()
#        textMapper.input=atext.output
        #textActor = tvtk.Actor()
#        textActor = tvtk.Follower()
#        textActor.mapper=textMapper
#        textActor.scale=scale
#        textActor.position=position
#        self._pf.scene.add_actor(textActor)
#        if follower==True:
#            textActor.camera=self._pf.scene.camera
#        self._pf.reset()
#        return textActor
#    def add_caption(self,caption="",attachment_point=(0,0,0),border=False):
#        
#        cap=tvtk.CaptionActor2D(caption=caption,
#                                attachment_point=attachment_point,
#                                border=border, position2=(.12,.05))
#        #position2 is to make the caption smaller, but I dont really
#        #understand how it is working
#        self._pf.scene.add_actor(cap)
#        return cap



# Create the axes and the associated mapper and actor.
#def axes_actor(origin=(0,0,0), scale=(1.,1.,1.)):
#    axes = tvtk.Axes()
#    axes.origin=(0, 0, 0)
#    axesMapper = tvtk.PolyDataMapper()
#    axesMapper.input=axes.output 
#    #SetInputConnection(axes.GetOutputPort())
#    axesActor = tvtk.Actor()
#    axesActor.mapper=axesMapper
#    axesActor
#    return axesActor
