#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module with functions and classes to represent the pyoptools objects
in `jupyter notebooks <http://jupyter.org>`_.
"""
from IPython.display import display
from pyoptools.raytrace.system import System
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.surface import Surface
from pyoptools.misc.pmisc import wavelength2RGB, cross, rot_x, rot_y, rot_z

try:
    import pythreejs as py3js
except ModuleNotFoundError:
    print("need py3js installed to be able to plot systems in Jupyter notebooks")
from pyoptools.misc.pmisc import wavelength2RGB, cross
from numpy import pi, array, dot
from math import sqrt
from matplotlib import colors

__all__ = ["Plot3D"]


import pythreejs as pjs
import numpy as np

def create_transformed_scene():
    # Create a mesh with a basic geometry
    mesh = pjs.BoxBufferGeometry()

    # Create a group and add the mesh
    group = pjs.Group()
    group.add(mesh)

    # Disable automatic updates of the transformation matrix
    group.matrixAutoUpdate = False

    # Manually compute and set the transformation matrix
    D = [45, 30, 60]  # Rotation angles in degrees for X, Y, Z
    P = [1, 2, 3]     # Translation vectors for X, Y, Z
    transformation_matrix = create_transformation_matrix(D, P)
    group.matrix = transformation_matrix

    # Setup the scene and renderer
    scene = pjs.Scene(children=[group])
    camera = pjs.PerspectiveCamera(position=[2, 2, 2], fov=50)
    renderer = pjs.Renderer(camera=camera, scene=scene, controls=[pjs.OrbitControls(controlling=camera)])

    return renderer

def create_transformation_matrix(D, P):
    # Convert degrees to radians for rotation
    #theta, phi, psi = np.radians(D)  # Assuming D is given in degrees
    
    #theta, phi, psi = D
    
    psi,phi,theta = D
 
    # Rotation matrices
    Rz = np.array([
        [np.cos(theta), -np.sin(theta), 0, 0],
        [np.sin(theta), np.cos(theta), 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    Ry = np.array([
        [np.cos(phi), 0, np.sin(phi), 0],
        [0, 1, 0, 0],
        [-np.sin(phi), 0, np.cos(phi), 0],
        [0, 0, 0, 1]
    ])
    Rx = np.array([
        [1, 0, 0, 0],
        [0, np.cos(psi), -np.sin(psi), 0],
        [0, np.sin(psi), np.cos(psi), 0],
        [0, 0, 0, 1]
    ])

    # Translation matrix
    T = np.array([
        [1, 0, 0, P[0]],
        [0, 1, 0, P[1]],
        [0, 0, 1, P[2]],
        [0, 0, 0, 1]
    ])

    # Composite transformation matrix
    return T @ Rz @ Ry @ Rx


def surf2mesh(S, P=(0, 0, 0), D=(0, 0, 0), wire=False):

    color = "#ffff00"

    points, polylist = S.polylist()

    # Conversion para quethreejs la entienda

    polylist = list(polylist)

    lpoly = []
    lpoints = []

    for l in points:
        lpoints.append(list(l))

    for l in polylist:
        lpoly.append(list(map(int, l)))

    vertices = lpoints

    faces = lpoly

    # Map the vertex colors into the 'color' slot of the faces
    # Map the normals
    nfaces = []
    for f in faces:
        p0 = points[f[0]]
        p1 = points[f[1]]
        p2 = points[f[2]]
        v0 = array(p1) - array(p0)
        v1 = array(p2) - array(p0)
        v3 = cross(v0, v1)
        v3 = tuple(v3 / sqrt(v3[0] ** 2 + v3[1] ** 2 + v3[2] ** 2))

        nfaces.append(f + [v3, color, None])

    # Create the geometry:

    surfaceGeometry = py3js.Geometry(
        vertices=vertices,
        faces=nfaces,
        # colors=vertexcolors
    )

    # surfaceGeometry = py3js.SphereGeometry(radius=300, widthSegments=32, heightSegments=24)

    if wire:
        surfaceGeometry = py3js.WireframeGeometry(surfaceGeometry)

    # Calculate normals per face, for nice crisp edges:
    surfaceGeometry.exec_three_obj_method("computeFaceNormals")

    surfaceMaterial = py3js.MeshPhongMaterial(
        color=color,
        # ambient="#050505",
        specular="#ffffff",
        shininess=15,
        emissive="#000000",
        side="DoubleSide",
        transparent=True,
        opacity=0.8,
    )
    # surfaceMaterial = py3js.MeshLambertMaterial(color='red',side='DoubleSide')

    # Create a mesh. Note that the material need to be told to use the vertex colors.
    surfaceMesh = py3js.Mesh(
        geometry=surfaceGeometry,
        material=surfaceMaterial,
    )

    group = py3js.Group()

    group.add(surfaceMesh)
    group.matrixAutoUpdate = False

    transformation_matrix = create_transformation_matrix(D, P)

    group.matrix = tuple(transformation_matrix.flatten(order="F"))


    #surfaceMesh.position = tuple(P)
    #surfaceMesh.rotateZ(D[2])
    #surfaceMesh.rotateY(D[1])
    #surfaceMesh.rotateX(D[0])

    #group.position = tuple(P)
    #group.rotateZ(D[2])
    #group.rotateY(D[1])
    #group.rotateX(D[0])



    #return surfaceMesh
    return group


def comp2mesh(C, P, D):
    c = py3js.Group()
    if isinstance(C, Component):
        for surf in C.surflist:
            sS, sP, sD = surf
            s = surf2mesh(sS, sP, sD)
            c.add(s)

    elif isinstance(C, System):
        for comp in C.complist:
            sC, sP, sD = comp
            c.add(comp2mesh(sC, sP, sD))

    c.matrixAutoUpdate = False

    transformation_matrix = create_transformation_matrix(D, P)

    c.matrix = tuple(transformation_matrix.flatten(order="F"))

    #c.position = tuple(P)
    #c.rotateZ(D[2])
    #c.rotateY(D[1])
    #c.rotateX(D[0])

    return c

def ray2list(ray):
    rays = []

    P1 = ray.pos
    if len(ray.childs) > 0:
        P2 = ray.childs[0].pos
    else:
        P2 = P1 + 10.0 * ray.dir

    if ray.intensity != 0:

        line = [list(P1), list(P2)]
        rays.append(line)

    for i in ray.childs:
        rays.extend(ray2list(i))
    return rays


def ray2mesh(ray):
    rays = py3js.Group()

    if ray.draw_color is None:
        color = wavelength2RGB(ray.wavelength)
    else:
        color = colors.to_rgb(ray.draw_color)

    int_colors = [int(255 * c) for c in color]
    material = py3js.LineBasicMaterial(color="#{:02X}{:02X}{:02X}".format(*int_colors))

    rl = ray2list(ray)

    for r in rl:
        geometry = py3js.Geometry()
        geometry.vertices = r
        line = py3js.Line(geometry, material)
        rays.add(line)

    return rays


# def ray2mesh(ray):
#    rays=py3js.Group()

#    P1 = ray.pos
#    w = ray.wavelength
#    rc, gc, bc = wavelength2RGB(w)
#    rc=int(255*rc)
#    gc=int(255*gc)
#    bc=int(255*bc)
#    material = py3js.LineBasicMaterial(color = "#{:02X}{:02X}{:02X}".format(rc,gc,bc))

#    if len(ray.childs) > 0:
#        P2 = ray.childs[0].pos
#    else:
#        P2 = P1 + 10. * ray.dir

#    if ray.intensity != 0:

#        geometry = py3js.Geometry()

#        geometry.vertices =  [list(P1),list(P2)]

#        line = py3js.Line( geometry, material)

#        rays.add(line)

#    for i in ray.childs:
#        rays.add(ray2mesh(i))
#    return rays


def sys2mesh(os):
    s = py3js.Group()
    if os is not None:
        for i in os.prop_ray:
            s.add(ray2mesh(i))
        # Draw Components
        n = 0
        for comp in os.complist:
            C, P, D = comp
            c = comp2mesh(C, P, D)
            s.add(c)
    return s


def Plot3D(
    S, size=(800, 200), center=(0, 0, 0), rot=[(pi / 3.0, pi / 6.0, 0)],
    scale=1):
    """Function to create 3D interactive visualization widgets in a jupyter
    notebook

    Args:
        S: (:class:`~pyoptools.raytrace.system.System`,
            :class:`~pyoptools.raytrace.component.Component` or
            :class:`~pyoptools.raytrace.component.Component`) Object to plot
        size: (Tuple(float,float)) Field of view in X and Y for the window
            shown in the notebook.
        center: (Tuple(float,float,float) Coordinate of the center of the
            visualization window given in the coordinate system of the object
            to plot.
        rot:   List of tuples. Each tuple describe an (Rx, Ry, Rz) rotation and
               are applied in order to generate the first view of the window.
        scale: (float)  Scale factor applied to the rendered window
    Returns:
        pyjs renderer needed to show the image in the jupiter notebook.

    """
    width, height = size

    light = py3js.DirectionalLight(
        color="#ffffff", intensity=0.7, position=[0, 1000, 0]
    )
    alight = py3js.AmbientLight(
        color="#777777",
    )

    # Set up a scene and render it:
    #cam = py3js.PerspectiveCamera(position=[0, 0, 500], fov=70, children=[light], aspect=width / height)

    pos = array((0, 0, 500))

    for r in rot:
        pos = dot(rot_z(r[2]), pos)
        pos = dot(rot_y(r[1]), pos)
        pos = dot(rot_x(r[0]), pos)

    #cam = py3js.OrthographicCamera(
    #    -width / 2 * scale,
    #    width / 2 * scale,
    #    height / 2 * scale,
    #    -height / 2 * scale,
    #    children=[light],
    #    position=list(pos),
    #    zoom=scale,
    #)

    cam = py3js.OrthographicCamera(
        -width / 2 * scale,
        width / 2 * scale,
        height / 2 * scale,
        -height / 2 * scale,
        children=[light],
        position=list(pos),
        zoom=scale,
    )

    if isinstance(S, System):
        c = sys2mesh(S)
    elif isinstance(S, Component):
        c = comp2mesh(S, (0, 0, 0), (0, 0, 0))
    else:
        c = surf2mesh(S, (0, 0, 0), (0, 0, 0))

    scene = py3js.Scene(children=[c, alight, cam], background="#000000")
    oc = py3js.OrbitControls(controlling=cam)
    oc.target = center
    renderer = py3js.Renderer(
        camera=cam,
        background="black",
        background_opacity=1,
        scene=scene,
        controls=[oc],
        width=width * scale,
        height=height * scale,
    )
    
    return renderer
