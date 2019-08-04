#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pyoptools.raytrace.system import System
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Surface
from pyoptools.misc.pmisc import wavelength2RGB, cross
try:
    import pythreejs as py3js
except ModuleNotFoundError:
    print("need py3js installed to be able to plot systems in Jupyter notebooks")
from pyoptools.misc.pmisc import wavelength2RGB, cross
from numpy import pi
def surf2mesh(S,P=(0,0,0),D=(0,0,0),wire=False):
    
    color="#ffff00"
        
    points,polylist = S.polylist()

    #Conversion para quethreejs la entienda

    polylist=list(polylist)

    lpoly=[]
    lpoints=[]

    for l in points:
        lpoints.append(list(l))

    for l in polylist:
        lpoly.append(list(map(int,l)))
    
    vertices = lpoints

    faces = lpoly

    # Map the vertex colors into the 'color' slot of the faces
    # Map the normals
    nfaces=[]
    for f in faces:
        p0 = points[f[0]]
        p1 = points[f[1]]
        p2 = points[f[2]]
        v0 = array(p1)-array(p0)
        v1 = array(p2)-array(p0)
        v3 = cross(v0, v1)
        v3 = tuple(v3 / sqrt(v3[0]**2 + v3[1]**2 + v3[2]**2))
        
        nfaces.append(f + [v3, color, None])
        
    # Create the geometry:
    
    surfaceGeometry = py3js.Geometry(vertices=vertices,
        faces=nfaces,
        #colors=vertexcolors
                           )
    
    
    #surfaceGeometry = py3js.SphereGeometry(radius=300, widthSegments=32, heightSegments=24)
    
    if wire:
        surfaceGeometry = py3js.WireframeGeometry(surfaceGeometry)
        
    # Calculate normals per face, for nice crisp edges:
    surfaceGeometry.exec_three_obj_method('computeFaceNormals')

    surfaceMaterial=py3js.MeshPhongMaterial( color=color,
                                             ambient="#050505", 
                                             specular="#ffffff",
                                             shininess= 15,
                                             emissive="#000000",
                                             side='DoubleSide',
                                             transparent = True,
                                             opacity=.8)
    #surfaceMaterial = py3js.MeshLambertMaterial(color='red',side='DoubleSide')
    
    # Create a mesh. Note that the material need to be told to use the vertex colors.
    surfaceMesh = py3js.Mesh(
        geometry=surfaceGeometry,
        material= surfaceMaterial,)
    
    surfaceMesh.rotation=*D,"ZYX"
    surfaceMesh.position=tuple(P)
    return surfaceMesh

def comp2mesh(C, P, D):
    c=py3js.Group()
    if isinstance(C, Component):
        for surf in C.surflist:
            sS, sP, sD = surf
            s=surf2mesh(sS,sP,sD)
            c.add(s)
    
    elif isinstance(C, System):
       for comp in C.complist:
            sC, sP, sD = comp
            c.add(comp2mesh(sC, sP, sD))
    #glPopMatrix()
    c.rotation=*D,"ZYX"
    c.position=tuple(P)
    return c

def ray2mesh(ray):
    rays=py3js.Group()

    P1 = ray.pos
    w = ray.wavelength
    rc, gc, bc = wavelength2RGB(w)
    rc=int(255*rc)
    gc=int(255*gc)
    bc=int(255*bc)
    material = py3js.LineBasicMaterial(color = "#{:02X}{:02X}{:02X}".format(rc,gc,bc))

    if len(ray.childs) > 0:
        P2 = ray.childs[0].pos
    else:
        P2 = P1 + 10. * ray.dir
    
    if ray.intensity != 0:
        
        geometry = py3js.Geometry()
    
        geometry.vertices =  [list(P1),list(P2)]
        
        line = py3js.Line( geometry, material)
        
        rays.add(line)
    
    for i in ray.childs:
        rays.add(ray2mesh(i))
    return rays


def sys2mesh(os):
    s=py3js.Group()
    if os is not None:
        for i in os.prop_ray:
            s.add(ray2mesh(i))
        # Draw Components
        n=0
        for comp in os.complist:    
            C, P, D = comp
            c=comp2mesh(C, P, D)
            s.add(c)
    return s
    
def Plot3D(S,size=(800,200)):
    width,height=size
    
    
    light =  py3js.DirectionalLight(color='#ffffff',
                                    intensity=.7,
                                    position=[0, 1000,0])
    alight =  py3js.AmbientLight(color='#777777',)
    
    # Set up a scene and render it:
    #cam = py3js.PerspectiveCamera(position=[0, 0, 500], fov=70, children=[light], aspect=width / height)


    cam = py3js.OrthographicCamera(-width/2,width/2, height/2, -height/2,children=[light],position=[-500, 00,0])

    c=sys2mesh(S)

    scene = py3js.Scene(children=[c, alight,cam],background="#000000")

    renderer = py3js.Renderer(camera=cam, background='black', background_opacity=1,
                          scene=scene, controls=[py3js.OrbitControls(controlling=cam)],width=width, height=height)

    return(renderer)

    
    
