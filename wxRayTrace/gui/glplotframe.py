from wxRayTrace.gui.glwindow import wxAdvancedGLWindow
from wx.glcanvas import WX_GL_DOUBLEBUFFER, WX_GL_RGBA
#~ from wxRayTrace.gui.glplotframe2 import glPlotFrame2 as glPlotFrame

import sys,math

from numpy import array, sqrt, dot, pi
from wx import Frame, DefaultPosition, Size

from OpenGL.GL import *
from OpenGL.GLU import *

from pyoptools.raytrace.system import System
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Surface
from pyoptools.misc.pmisc import wavelength2RGB

    
#def glPlotFrame(os=None):
#    pf= AsyncCall(glFrame, os)
#    return pf.Wait()
    

#~ class glPlotFrame(Frame):
    #~ def __init__(self, os=None):
        #~ if not isinstance(os,System):
            #~ raise TypeError
        #~ Frame.__init__(self, None, -1, '3D System', DefaultPosition, Size(400,400))
        #~ canvas=glCanvas(self, os)
        #~ self.Show()

class glCanvas(wxAdvancedGLWindow):    
    def __init__(self,parent,  os=None):
        wxAdvancedGLWindow.__init__(self, parent, attribList=[WX_GL_DOUBLEBUFFER, WX_GL_RGBA])
        self.os=os
        self.scene=None

    def InitGL(self):
        #self.set_base_distance(500)
        #self.set_distance(500)
        #set up lighting
        glLightfv(GL_LIGHT0, GL_AMBIENT, [1.0, 1.0, 1.0, 1.0])
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glClearColor(0.7,0.7,0.7,0.0)
        glShadeModel(GL_SMOOTH)
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_COLOR_MATERIAL)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        self.scene = glGenLists(1)
        glNewList(self.scene,GL_COMPILE);
        self.DrawGLL()
        glEndList();

    def DrawGL(self):
        if self.scene!=None:
            glCallList(self.scene)

    def DrawGLL(self):
        #Draw Rays
        # "RRRRRRRRRRRRRRRRRRRRRRxR", self.os
        if self.os != None:
			for i in self.os.prop_ray:
				self.DrawRay(i)
			#Draw Components
			for comp in self.os.complist:
				self.DrawComp(comp)
                  
    def DrawComp(self, comp):
        C,P,D = comp
        #glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glTranslatef(P[0],P[1],P[2])
        glRotatef(180*D[2]/pi,0.0,0.0,1.0)
        glRotatef(180*D[1]/pi,0.0,1.0,0.0)
        glRotatef(180*D[0]/pi,1.0,0.0,0.0)
        if isinstance(C, Component):
            for surf in C.surflist:
                S, P, D=surf
                self.DrawSurf(S, P, D)
        elif isinstance(C, System):
            for comp in C.complist:
                self.DrawComp(comp)
        glPopMatrix()
 
    def DrawSurf(self, surf, P, D):
        if isinstance(surf, Surface):
            #points, polylist =surf.shape.polylist(surf.topo)
            points, polylist =surf.polylist()
            glPushMatrix()
            glTranslatef(P[0],P[1],P[2])    
            glRotatef(180*D[2]/pi,0.0,0.0,1.0)
            glRotatef(180*D[1]/pi,0.0,1.0,0.0)
            glRotatef(180*D[0]/pi,1.0,0.0,0.0)
            glColor4f(.7,.7,.7, 0.5)
            for p in polylist:
                if len(p)==3:
                    p0=points[p[0]]
                    p1=points[p[1]]
                    p2=points[p[2]]
                    glBegin(GL_TRIANGLES) #Drawing Using Triangles
                    glVertex3f( p0[0], p0[1], p0[2])
                    glVertex3f( p1[0], p1[1], p1[2])
                    glVertex3f( p2[0], p2[1], p2[2])
                    glEnd()                
                elif len(p)==4:
                    p0=points[p[0]]
                    p1=points[p[1]]
                    p2=points[p[2]]
                    p3=points[p[3]]
                    glBegin(GL_QUADS)           # Start Drawing The Cube
                    glVertex3f( p0[0], p0[1], p0[2])
                    glVertex3f( p1[0], p1[1], p1[2])
                    glVertex3f( p2[0], p2[1], p2[2])
                    glVertex3f( p3[0], p3[1], p3[2])
                    glEnd()                
            glPopMatrix()
 
    def DrawRay(self, ray):
        P1=ray.pos
        w=ray.wavelength
        rc,gc,bc=wavelength2RGB(w)
        if len(ray.childs)>0:
            P2=ray.childs[0].pos
        else:
            P2=P1+10.*ray.dir

        if ray.intensity!=0:
            glBegin(GL_LINES)
            glColor4f(rc,gc,bc, 1.)
            glVertex3f( P1[0], P1[1], P1[2])
            glVertex3f( P2[0], P2[1], P2[2])
            glEnd()
        for i in ray.childs:
            self.DrawRay(i)
         
