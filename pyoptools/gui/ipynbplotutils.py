#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Ipython Notebook specific plotting utilities
Se hizo una correccion siguiendo
http://g.sweyla.com/blog/2014/osmesa-pyopengl-310/
"""

import six
from numpy import array, pi, sqrt, degrees
import PIL.Image as PImage

from OpenGL.GL import *
from OpenGL.GLU import *

from OpenGL import arrays
try:
    from OpenGL.platform import CurrentContextIsValid
    from OpenGL.raw.osmesa.mesa import (OSMesaCreateContext,
                                        OSMesaMakeCurrent,
                                        OSMesaDestroyContext)

except:
    print("need OSMesa installed, and the following environment variable")
    print("export PYOPENGL_PLATFORM=osmesa")

from pyoptools.raytrace.system import System
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Surface
from pyoptools.misc.pmisc import wavelength2RGB, cross

# Nota toca exportar la veriable de ambiente
# export PYOPENGL_PLATFORM=osmesa


def draw_sys(os):
    if os is not None:
        for i in os.prop_ray:
            draw_ray(i)
        # Draw Components
        for comp in os.complist:
            C, P, D = comp
            draw_comp(C, P, D)


def draw_comp(C, P, D):
    # C,P,D = comp
    # glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glTranslatef(P[0], P[1], P[2])
    glRotatef(180 * D[2] / pi, 0.0, 0.0, 1.0)
    glRotatef(180 * D[1] / pi, 0.0, 1.0, 0.0)
    glRotatef(180 * D[0] / pi, 1.0, 0.0, 0.0)

    if isinstance(C, Component):
        for surf in C.surflist:
            S, P, D = surf
            draw_surf(S, P, D)
    elif isinstance(C, System):
        for comp in C.complist:
            C, P, D = comp
            draw_comp(C, P, D)
    glPopMatrix()


def draw_surf(surf, P, D):
    if isinstance(surf, Surface):
        # points, polylist =surf.shape.polylist(surf.topo)
        points, polylist = surf.polylist()
        glPushMatrix()
        glTranslatef(P[0], P[1], P[2])
        glRotatef(180 * D[2] / pi, 0.0, 0.0, 1.0)
        glRotatef(180 * D[1] / pi, 0.0, 1.0, 0.0)
        glRotatef(180 * D[0] / pi, 1.0, 0.0, 0.0)
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                     [1., 1., 0, 0.7])
        # glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, [1.,1.,0.,1.])
        # glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, [1.,0.,0.,1.])
        for p in polylist:
            if len(p) == 3:
                p0 = points[p[0]]
                p1 = points[p[1]]
                p2 = points[p[2]]
                v0 = array(p1)-array(p0)
                v1 = array(p2)-array(p0)
                v3 = cross(v0, v1)
                v3 = v3 / sqrt(v3[0]**2 + v3[1]**2 + v3[2]**2)
                glBegin(GL_TRIANGLES)  # Drawing Using Triangles
                glNormal3f(v3[0], v3[1], v3[2])
                glVertex3f(p0[0], p0[1], p0[2])
                glVertex3f(p1[0], p1[1], p1[2])
                glVertex3f(p2[0], p2[1], p2[2])
                glEnd()
            elif len(p) == 4:
                p0 = points[p[0]]
                p1 = points[p[1]]
                p2 = points[p[2]]
                p3 = points[p[3]]
                v0 = array(p1)-array(p0)
                v1 = array(p2)-array(p0)
                v3 = cross(v0, v1)
                v3 = v3/sqrt(v3[0]**2+v3[1]**2+v3[2]**2)
                glBegin(GL_QUADS)  # Start Drawing The Cube
                glNormal3f(v3[0], v3[1], v3[2])
                glVertex3f(p0[0], p0[1], p0[2])
                glVertex3f(p1[0], p1[1], p1[2])
                glVertex3f(p2[0], p2[1], p2[2])
                glVertex3f(p3[0], p3[1], p3[2])
                glEnd()
        glPopMatrix()


def draw_ray(ray):
    P1 = ray.pos
    w = ray.wavelength
    rc, gc, bc = wavelength2RGB(w)
    if len(ray.childs) > 0:
        P2 = ray.childs[0].pos
    else:
        P2 = P1 + 10. * ray.dir

    if ray.intensity != 0:
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
                     [rc, gc, bc, 1.])
        glBegin(GL_LINES)
        glColor4f(rc, gc, bc, 1.)
        glVertex3f(P1[0], P1[1], P1[2])
        glVertex3f(P2[0], P2[1], P2[2])
        glEnd()
    for i in ray.childs:
        draw_ray(i)


class Plot3D(object):
    """
    Generate a 3D image of the optical element under study

    Parameters:
    ===========

    os     Optical element (system, component or surface) to be drawn
    center Tuple (x,y,z) with the coordinates of the center of the drawing
    size   Tuple (width, height) of the requested image
    rot    list of (rx,ry,rz) tuples containing a series of rotation angles
    scale  scale for the image

    The rotations are applied first rx, then ry and then rz

    Attributes:
    ===========
    buffer    3D buffer of the system
    image     PIL image of the buffer
    """
    def __init__(self, os,
                 center=(0, 0, 0),
                 size=(400, 400),
                 rot=[(0, 0, 0)],
                 scale=1.):

        self.buffer = self.buffer_3d(os, center=center,
                                     size=size, rot=rot, scale=scale)
        self.image = self._buffer_to_image()

    def show(self):
        '''Show the PIL image with the default viewer'''
        self.image.show()

    def _repr_jpeg_(self):
        '''Ipython can embed the jpg image directly
        from http://ipython.org/ipython-doc/dev/config/integrating.html
        '''
        temppng = six.BytesIO()
        self.image.save(temppng, "JPEG")
        data = temppng.getvalue()
        temppng.close()
        return data

    def _buffer_to_image(self):
        h, w, _c = self.buffer.shape
        image = PImage.frombytes("RGBA", (w, h), self.buffer)
        image = image.transpose(PImage.FLIP_TOP_BOTTOM)
        return image

    def buffer_3d(self, os, center=(0, 0, 0),
                  size=(400, 400), rot=[(0, 0, 0)],
                  scale=1.):

        left = -size[0] / 2
        right = size[0] / 2
        top = size[1] / 2
        bottom = -size[1] / 2

        ctx = OSMesaCreateContext(GL_RGBA, None)

        width, height = int(size[0] * scale), int(size[1] * scale)

        buf = arrays.GLubyteArray.zeros((height, width, 4))
        assert(OSMesaMakeCurrent(ctx, buf, GL_UNSIGNED_BYTE, width, height))
        assert(CurrentContextIsValid())

        light_ambient = [.5, 0.5, 0.5, 1.0]
        light_diffuse = [1.0, 1.0, 1.0, 1.0]
        light_specular = [1.0, 1.0, 1.0, 1.0]
        light_position = [0.0, 1.0, -1.0, 1.]

        glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient)
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse)
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular)

        glShadeModel(GL_SMOOTH)  # No parece estar funcionando
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_BLEND)
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)

        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()

        glOrtho(left, right, top, bottom, -1000.0, 10000.0)

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

        # Las rotaciones se aplican en orden inverso. No se por que, pero
        # esta funcionando
        for rx, ry, rz in rot[::-1]:
            glRotatef(degrees(rz), 0, 0, 1)
            glRotatef(degrees(ry), 0, 1, 0)
            glRotatef(degrees(rx), 1, 0, 0)

        glTranslatef(-center[0], -center[1], -center[2])
        glLightfv(GL_LIGHT0, GL_POSITION, light_position)

        if isinstance(os, System):
            draw_sys(os)
        elif isinstance(os, Component):
            draw_comp(os, (0, 0, 0), (0, 0, 0))
        elif isinstance(os, Surface):
            draw_surf(os, (0, 0, 0), (0, 0, 0))

        glFinish()
        return buf
