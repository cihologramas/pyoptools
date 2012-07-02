import wx
from StringIO import StringIO
import Image as PIL
from IPython.core.display import Image

from OpenGL.GL import *
from OpenGL.GLU import *

from pyoptools.raytrace.system import System
from oglframe import OGLFrame
from types import NoneType
from math import sqrt
from pyoptools.gui.glwindow import v3distsq,glTranslateScene,glRotateScene

#~ import Image

class glPlotFrame(OGLFrame):
    def __init__(self, os=None):
        if not isinstance(os,(System,NoneType)):
            raise TypeError
        OGLFrame.__init__(self, None, -1, '3D System')
        self.glcanvas.os=os
        self.Show()
    def translate(self,x,y):
        # Scale mouse translations to object viewplane so object tracks with mouse
        win_height = max( 1,self.glcanvas.w)
        obj_c      = (self.glcanvas.xcenter, self.glcanvas.ycenter, self.glcanvas.zcenter)
        
        win        = gluProject( obj_c[0], obj_c[1], obj_c[2] )
        obj        = gluUnProject( win[0], win[1] + 0.5 * win_height, win[2] )
        dist       = sqrt( v3distsq( obj, obj_c ) )
        scale      = abs( dist / ( 0.5 * win_height ) )

        glTranslateScene(scale, x, y, 0, 0)
        self.glcanvas.wxRedraw()
    
    def OnUp(self, event):
        self.translate(0,-10)
        event.Skip()
    
    def OnDown(self, event):
        self.translate(0,10)
        event.Skip()
    
    
    def OnLeft(self, event):
        self.translate(-10,0)
        event.Skip()

    def OnRight(self, event):
        self.translate(10,0)
        event.Skip()

    def scale(self, inout):
        """Scale the scene.  Achieved by moving the eye position."""
        scale = 1 - 0.01 * (inout)
        self.glcanvas.distance = self.glcanvas.distance * scale
        self.glcanvas.wxRedraw()
    def OnZmOut(self, event):
        self.scale(-5)
        event.Skip()
        
    def OnZmIn(self, event):
        self.scale(5)
        event.Skip()
        
    def rotate(self, rx, ry):
        glRotateScene(0.5,
                    self.glcanvas.xcenter, self.glcanvas.ycenter, self.glcanvas.zcenter,
                    rx, ry, 0, 0)
        self.glcanvas.wxRedraw()
    def OnSpUp(self, event):
        self.rotate(0,-2)
        event.Skip()
    def OnSpDown(self, event):
        self.rotate(0,2)
        event.Skip()
    def OnSpLeft(self, event):
        self.rotate(-2,0)
        event.Skip()
    def OnSpRight(self, event):
        self.rotate(2,0)
        event.Skip()
    
    def spin(self,d):
        # rotate about z
        sz = self.glcanvas.GetClientSizeTuple()
        sz = (sz[0]/2, sz[1]/2)
        #~ xp = event.GetX()
        #~ yp = event.GetY()
        #~ dy = (self.ymouse-yp)
        #~ dx = (self.xmouse-xp)
        #~ if yp > sz[1]:
        #~ dx = dx * -1
        #~ if xp < sz[0]:
        #~ dy = dy * -1
        #~ d = dx + dy
        glMatrixMode(GL_MODELVIEW);
        m = glGetDouble(GL_MODELVIEW_MATRIX)
        glLoadIdentity()
        glTranslatef(self.glcanvas.xcenter,self.glcanvas.ycenter,self.glcanvas.zcenter)
        glRotatef(.5*d,0,0,1.)
        glTranslatef(-self.glcanvas.xcenter,-self.glcanvas.ycenter,-self.glcanvas.zcenter)
        glMultMatrixd(m)
        self.glcanvas.wxRedraw()
    
    def OnSpCCW(self,event):
        self.spin(2)
        event.Skip()
    def OnSpCW(self,event):
        self.spin(-2)
        event.Skip()
    def GetImageData(self):
        x,y,width,height = glGetDoublev(GL_VIEWPORT)
        x=int(x)
        y=int(y)
        width=int(width)
        height=int(height)
        glPixelStorei(GL_PACK_ALIGNMENT, 1)
        data = glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE)
        image = PIL.fromstring( "RGB", (width, height), data )
        image = image.transpose( PIL.FLIP_TOP_BOTTOM)
        temppng=StringIO()
        image.save(temppng, "PNG" )
        data=temppng.getvalue()
        temppng.close()
        #~ print 'Saved image to %s'% (os.path.abspath( filename))
        #~ return image
        return Image(data,embed=True)
        
if __name__ == "__main__":
    from pyoptools.all import *
    # Blue
    r_b= parallel_beam_c(size=(10,10),num_rays=(10,10), wavelength=.470)
    N_BK7=schott['BK7']
    bs=BeamSplitingCube(size=50,material=N_BK7,reflectivity=0.5)
    #Definition of a detector plane
    ccd=CCD(size=(20, 20))
    os=System(complist=[(bs,(0,0,50),(0,0,0)),
                    (ccd,(0,0,100),(0,0,0)),
                    ],n=1)

    #Add the ray sources
    os.ray_add(r_b)
    #Propagate the rays
    os.propagate()

    
    
    app = wx.PySimpleApp(0)
    wx.InitAllImageHandlers()
    oglFrame = glPlotFrame2(os)
    app.SetTopWindow(oglFrame)
    oglFrame.Show()
    app.MainLoop()

