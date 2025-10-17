
"""Module with functions and classes to represent the pyoptools objects
in `jupyter notebooks <http://jupyter.org>`_.
"""

from pyoptools.raytrace.system import System
from pyoptools.raytrace.component import Component
from pyoptools.misc.pmisc import wavelength2RGB, cross, rot_x, rot_y, rot_z
from numpy import array

try:
    import pythreejs as py3js
except ModuleNotFoundError:
    print("need pythreejs installed to be able to plot systems in Jupyter notebooks")

from numpy import pi, array, dot, sin, cos
from math import sqrt
from matplotlib import colors

__all__ = ["Plot3D"]


def create_transformation_matrix(P, D):
    """
    Create a composite transformation matrix from translation and rotation parameters.

    This function generates a 4x4 transformation matrix that combines both rotation and translation.
    The rotation is applied first, followed by the translation. The rotation is defined by Euler angles
    psi (rotation around x-axis), phi (rotation around y-axis), and theta (rotation around z-axis).

    Parameters
    ----------
    P : tuple
        A 3-element tuple representing the translation vector (Tx, Ty, Tz).
    R : tuple
        A 3-element tuple representing the Euler angles (psi, phi, theta) in radians.

    Returns
    -------
    numpy.ndarray
        A 4x4 numpy array representing the composite transformation matrix.

    Notes
    -----
    The rotation matrices are defined in the order XYZ, which means the rotation around the x-axis is applied first,
    followed by the rotation around the y-axis, and finally the rotation around the z-axis. The translation is applied
    after all rotations have been applied.
    """

    psi, phi, theta = D

    # Rotation matrices
    Rz = array(
        [
            [cos(theta), -sin(theta), 0, 0],
            [sin(theta), cos(theta), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ]
    )
    Ry = array(
        [
            [cos(phi), 0, sin(phi), 0],
            [0, 1, 0, 0],
            [-sin(phi), 0, cos(phi), 0],
            [0, 0, 0, 1],
        ]
    )
    Rx = array(
        [
            [1, 0, 0, 0],
            [0, cos(psi), -sin(psi), 0],
            [0, sin(psi), cos(psi), 0],
            [0, 0, 0, 1],
        ]
    )

    # Translation matrix
    T = array([[1, 0, 0, P[0]], [0, 1, 0, P[1]], [0, 0, 1, P[2]], [0, 0, 0, 1]])

    # Composite transformation matrix
    return T @ Rz @ Ry @ Rx


def surf2mesh(S, P=(0, 0, 0), D=(0, 0, 0), wire=False):
    """
    Convert a pyOpTools surface to a mesh representation suitable for rendering
    with pythreejs.

    This function takes a surface object, which is expected to have a method
    `polylist` that returns a list of points and polygons defining the surface.
    It then creates a mesh representation of this surface, which can be
    displayed in a Jupyter notebook using pythreejs. The function allows for the
    mesh to be translated and rotated according to the provided parameters, and
    optionally displayed as a wireframe.

    Parameters
    ----------
    S : object
        The surface object to be converted into a mesh. This object must have a
        `polylist` method that returns two elements: a list of points and a list
         of polygons.
    P : tuple of float, optional
        A 3-element tuple representing the translation vector (Tx, Ty, Tz) to be
        applied to the mesh.
        Default is (0, 0, 0), meaning no translation.
    D : tuple of float, optional
        A 3-element tuple representing the Euler angles (psi, phi, theta) in
        radians for rotation of the mesh. Default is (0, 0, 0), meaning no
        rotation.
    wire : bool, optional
        If True, the mesh is created as a wireframe. Otherwise, a solid mesh is
        created. Default is False.

    Returns
    -------
    pythreejs.Group
        A pythreejs Group object containing the mesh representation of the
        surface, which can be directly added to a pythreejs scene for rendering.

    Notes
    -----
    The function computes normals for the faces of the mesh to ensure proper
    lighting and shading when rendered. The color of the mesh is hardcoded to
    yellow, and the material properties are set to create a somewhat shiny and
    semi-transparent appearance.
    """

    color = "#ffff00"

    points, polylist = S.polylist()

    polylist = list(polylist)

    lpoly = []
    lpoints = []

    for point in points:
        lpoints.append(list(point))

    for poly in polylist:
        lpoly.append(list(map(int, poly)))

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
    )

    if wire:
        surfaceGeometry = py3js.WireframeGeometry(surfaceGeometry)

    # Calculate normals per face, for nice crisp edges:
    surfaceGeometry.exec_three_obj_method("computeFaceNormals")

    surfaceMaterial = py3js.MeshPhongMaterial(
        color=color,
        specular="#ffffff",
        shininess=15,
        emissive="#000000",
        side="DoubleSide",
        transparent=True,
        opacity=0.8,
    )

    # Create a mesh. Note that the material need to be told to use the vertex colors.
    surfaceMesh = py3js.Mesh(
        geometry=surfaceGeometry,
        material=surfaceMaterial,
    )

    group = py3js.Group()

    group.add(surfaceMesh)
    group.matrixAutoUpdate = False

    transformation_matrix = create_transformation_matrix(P, D)

    group.matrix = tuple(transformation_matrix.flatten(order="F"))

    return group


def comp2mesh(C, P, D):
    """
    Convert a pyOpTools Component or (Sub) System into a mesh representation for
    rendering with pythreejs.

    This function takes a Component or System object from pyOpTools and converts
    it into a mesh representation using pythreejs. This enables the visualization
    of the object in a Jupyter notebook. The function is capable of handling both
    individual components and systems composed of multiple components. It applies
    a transformation to the mesh based on the provided translation and rotation
    parameters.

    Parameters
    ----------
    C : Component or System
        The pyOpTools Component or System object to be converted into a mesh. If
        a Component is provided, the function will convert its surfaces into
        meshes. If a System is provided, the function will recursively convert all
        components within the system into meshes.
    P : tuple of float
        A 3-element tuple representing the translation vector (Tx, Ty, Tz) to be
        applied to the mesh. This moves the mesh to the specified position in the
        3D space.
    D : tuple of float
        A 3-element tuple representing the Euler angles (psi, phi, theta) in
        radians for rotation of the mesh. This rotates the mesh according to the
        specified angles.

    Returns
    -------
    pythreejs.Group
        A pythreejs Group object containing the mesh representation of the
        Component or System. This group can be directly added to a pythreejs scene
        for rendering in a Jupyter notebook.

    Notes
    -----
    The function disables automatic matrix updates for the returned group to
    ensure that the manually set transformation matrix remains unchanged. This is
    crucial for maintaining the correct position and orientation of the mesh in
    the scene.
    """
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

    transformation_matrix = create_transformation_matrix(P, D)

    c.matrix = tuple(transformation_matrix.flatten(order="F"))

    return c


def ray2list(ray):
    """
    Convert a ray and its child rays into a list of line segments.

    This function recursively processes a ray and its children (if any) to extract
    their positions and directions, converting them into a list of line segments.
    Each line segment is represented by a pair of points (start and end), where
    each point is a list of its coordinates. The function is designed to work with
    rays that have a tree-like structure, where each ray can spawn multiple child rays.

    Parameters
    ----------
    ray : object
        The ray object to be converted. This object must have attributes `pos` (position),
        `dir` (direction), `childs` (list of child rays), and `intensity` (indicating
        the ray's intensity).

    Returns
    -------
    list
        A list of line segments, where each line segment is represented as a list
        containing two lists, each of which represents the coordinates of the start
        and end points of the segment.

    Notes
    -----
    - The function assumes that the ray object and its children have the necessary
      attributes (`pos`, `dir`, `childs`, `intensity`) correctly defined.
    - If a ray has no children, its direction is extended by a factor of 10.0 units
      to create the end point of the line segment.
    - Only rays with a non-zero intensity are considered for conversion into line segments.
    """
    rays = []

    P1 = ray.origin
    if len(ray.childs) > 0:
        P2 = array(ray.childs[0].origin)
    else:
        P2 = P1 + 10.0 * array(ray.direction)

    if ray.intensity != 0:

        line = [list(P1), list(P2)]
        rays.append(line)

    for i in ray.childs:
        rays.extend(ray2list(i))
    return rays


def ray2mesh(ray):
    """
    Convert a ray object into a pythreejs Group containing line representations
    for visualization.

    This function takes a ray object, which may include child rays, and
    converts it into a collection of line segments that can be rendered using
    pythreejs. Each line segment represents the path of the ray or one of its
    child rays. The color of the lines is determined by the wavelength of the
    ray, if `draw_color` is not specified, or by the specified `draw_color`.

    Parameters
    ----------
    ray : object
        The ray object to be converted into a mesh. This object must have
        attributes such as `wavelength` and optionally `draw_color` to
        determine the color of the line. It should also support a structure
        for child rays, allowing recursive processing.

    Returns
    -------
    pythreejs.Group
        A pythreejs Group object containing line representations of the ray and
        its child rays. Each line is created with a material color based on the
        ray's wavelength or its specified `draw_color`, allowing for visual
        differentiation of rays in a rendered scene.

    Notes
    -----
    - The function uses the `ray2list` function to convert the ray and its
      children into a list of line segments.
    - The color of each line is determined by converting the ray's wavelength
      to RGB values if `draw_color` is not provided. If `draw_color` is
      provided, it is converted to RGB values directly.
    - The function creates a `LineBasicMaterial` with the calculated color and
      applies it to each line segment, adding them to a pythreejs Group for
      easy rendering.
    """
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


def sys2mesh(os):
    """
    Convert an optical system into a mesh representation for rendering with pythreejs.

    This function takes an optical system object, which contains properties of rays and components,
    and converts it into a mesh representation. This allows for the visualization of the entire
    optical system, including its rays and components, in a 3D space using pythreejs. The function
    iterates over the rays and components of the system, converting each into its mesh representation
    and adding them to a pythreejs Group object.

    Parameters
    ----------
    os : System
        The optical system object to be converted into a mesh. This object should contain properties
        of rays (`prop_ray`) and a list of components (`complist`), where each component is described
        by its geometry and transformation parameters.

    Returns
    -------
    pythreejs.Group
        A pythreejs Group object containing the mesh representations of the optical system's rays
        and components. This group can be directly added to a pythreejs scene for rendering in a
        Jupyter notebook.

    Notes
    -----
    The function checks if the optical system object is not None before proceeding with the conversion.
    It first converts the rays of the system using the `ray2mesh` function and then iterates over the
    components of the system, converting each component into a mesh using the `comp2mesh` function.
    """
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
    S, size=(800, 200), center=(0, 0, 0), rot=[(pi / 3.0, pi / 6.0, 0)], scale=1
):
    """
    Creates a 3D interactive visualization of an optical system, component,
    or surface within a Jupyter notebook using pythreejs.

    This function sets up a 3D scene with the specified object at its center,
    applying rotation and scaling as defined. It utilizes pythreejs to render
    the scene, allowing for interactive exploration of the object in 3D space.
    The visualization includes directional and ambient lighting to enhance the
    appearance of the object.

    Parameters
    ----------
    S : System or Component
        The optical system, component, or surface object to be visualized. This
        object must be compatible with the pyOpTools framework and capable of
        being converted into a mesh representation for 3D rendering.
    size : tuple of float, optional
        The size of the visualization window in the notebook, specified as
        (width, height) in pixels. Default is (800, 200).
    center : tuple of float, optional
        The coordinates of the center of the visualization window in the
        object's coordinate system, specified as (x, y, z). Default is (0, 0, 0),
        which centers the view on the origin of the object's coordinate system.
    rot : list of tuples, optional
        A list of rotation tuples, each representing a rotation around the x, y,
        and z axes (in radians). These rotations are applied in sequence to the
        initial view of the object. Default is [(pi / 3.0, pi / 6.0, 0)], which
        applies a specific initial rotation to the object.
    scale : float, optional
        A scale factor applied to the rendered scene, affecting the zoom level
        of the visualization. Default is 1, which renders the scene at a 1:1 scale.

    Returns
    -------
    pythreejs.Renderer
        A pythreejs Renderer object that displays the 3D visualization of the
        specified object within a Jupyter notebook. This renderer is interactive,
        allowing users to rotate, zoom, and pan the view within the notebook.

    Examples
    --------
    >>> Plot3D(my_system, size=(800, 600), center=(0, 0, 0),
    ... rot=[(pi/3, pi/6, 0)], scale=1)
        This example creates a 3D visualization of 'my_system' with a window
        size of 800x600 pixels, centered at the origin, with an initial rotation
        and a scale factor of 1.
    """
    width, height = size

    light = py3js.DirectionalLight(
        color="#ffffff", intensity=0.7, position=[0, 1000, 0]
    )
    alight = py3js.AmbientLight(
        color="#777777",
    )

    pos = array((0, 0, 500))

    for r in rot:
        pos = dot(rot_z(r[2]), pos)
        pos = dot(rot_y(r[1]), pos)
        pos = dot(rot_x(r[0]), pos)

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
