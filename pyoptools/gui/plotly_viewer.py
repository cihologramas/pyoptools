"""Interactive 3D visualization for pyOpTools optical systems using Plotly.

Provides a modern, universal WebGL 3D viewer compatible with Marimo,
JupyterLab 4, VS Code Interactive Window, and standalone HTML exports.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from matplotlib import colors

from pyoptools.misc.pmisc import wavelength2RGB

__all__ = ["plot_system_plotly"]

try:
    import plotly.graph_objects as go
except ImportError:
    go = None


def _transformation_matrix(
    translation: tuple[float, float, float] | list[float] | np.ndarray,
    rotation: tuple[float, float, float] | list[float] | np.ndarray,
) -> np.ndarray:
    """Creates a 4x4 homogeneous transformation matrix from translation and Euler rotations.

    Rotations are applied in the order: Rx(psi) -> Ry(phi) -> Rz(theta), followed by translation.
    Matches the exact coordinate transform convention in pyOpTools.
    """
    psi, phi, theta = rotation
    Rz = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0.0, 0.0],
            [np.sin(theta), np.cos(theta), 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    Ry = np.array(
        [
            [np.cos(phi), 0.0, np.sin(phi), 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [-np.sin(phi), 0.0, np.cos(phi), 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    Rx = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, np.cos(psi), -np.sin(psi), 0.0],
            [0.0, np.sin(psi), np.cos(psi), 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    T = np.array(
        [
            [1.0, 0.0, 0.0, translation[0]],
            [0.0, 1.0, 0.0, translation[1]],
            [0.0, 0.0, 1.0, translation[2]],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=float,
    )
    return T @ Rz @ Ry @ Rx


def _extract_ray_segments(ray) -> list[tuple[np.ndarray, np.ndarray, str]]:
    """Recursively extract visible ray segments and their display colors.
    Ray colors are computed from the wavelength, unless they have been manually
    specified, in which case matplotlib compatible color definitions are
    resolved.
    """
    segments = []
    p1 = np.array(ray.origin, dtype=float)
    if len(ray.childs) > 0:
        p2 = np.array(ray.childs[0].origin, dtype=float)
    else:
        p2 = p1 + 15.0 * np.array(ray.direction, dtype=float)

    if ray.draw_color is None:
        color = wavelength2RGB(ray.wavelength)
    else:
        color = colors.to_rgb(ray.draw_color)
    color_hex = "#{:02X}{:02X}{:02X}".format(*(int(255 * c) for c in color))
    if getattr(ray, "intensity", 1.0) != 0:
        segments.append((p1, p2, color_hex))

    for child in ray.childs:
        segments.extend(_extract_ray_segments(child))
    return segments


def _resolve_render_style(component, surface) -> tuple[str, float, bool, float]:
    """Resolve surface overrides against its owning component's render style.
    Render styles are owned by the component and can be resolved consistantly
    by different renderers.
    """
    color = surface.render_color
    if color is None:
        color = getattr(component, "render_color", "#CBD5E1")

    opacity = surface.render_opacity
    if opacity is None:
        opacity = getattr(component, "render_opacity", 0.5)

    smooth = surface.render_smooth
    if smooth is None:
        smooth = getattr(component, "render_smooth", False)

    thickness = surface.render_thickness
    if thickness is None:
        thickness = getattr(component, "render_thickness", 0.0)

    return color, opacity, smooth, thickness


def _extrude_planar_mesh(
    points: np.ndarray, triangles: np.ndarray, thickness: float
) -> tuple[np.ndarray, np.ndarray]:
    """Return a closed display-only slab for a planar triangular mesh.
    This is used when a cosmetic render thickness is specified, e.g. for
    screens and stops, so they are visible when viewed edge on.
    """
    if thickness == 0.0 or triangles.ndim != 2 or triangles.shape[1] != 3:
        return points, triangles

    normal = None
    for triangle in triangles:
        p0, p1, p2 = points[triangle]
        candidate = np.cross(p1 - p0, p2 - p0)
        magnitude = np.linalg.norm(candidate)
        if magnitude > 0.0:
            normal = candidate / magnitude
            break
    if normal is None:
        return points, triangles

    scale = max(np.ptp(points, axis=0).max(), 1.0)
    if not np.allclose((points - points[0]) @ normal, 0.0, atol=scale * 1e-9):
        return points, triangles

    oriented_triangles = triangles.copy()
    for index, triangle in enumerate(oriented_triangles):
        p0, p1, p2 = points[triangle]
        if np.dot(np.cross(p1 - p0, p2 - p0), normal) < 0.0:
            oriented_triangles[index, [1, 2]] = oriented_triangles[index, [2, 1]]

    edges: dict[tuple[int, int], list[tuple[int, int]]] = {}
    for triangle in oriented_triangles:
        for start, end in zip(triangle, np.roll(triangle, -1), strict=True):
            edges.setdefault(tuple(sorted((start, end))), []).append((start, end))
    if any(len(edge) > 2 for edge in edges.values()):
        return points, triangles

    vertex_count = len(points)
    offset = normal * (thickness / 2.0)
    slab_points = np.vstack((points + offset, points - offset))
    front = oriented_triangles
    back = oriented_triangles[:, [0, 2, 1]] + vertex_count
    sides = []
    for edge in edges.values():
        if len(edge) != 1:
            continue
        start, end = edge[0]
        sides.extend(
            (
                (start, end + vertex_count, end),
                (start, start + vertex_count, end + vertex_count),
            )
        )

    if sides:
        slab_triangles = np.vstack((front, back, np.asarray(sides, dtype=int)))
    else:
        slab_triangles = np.vstack((front, back))
    return slab_points, slab_triangles


def _iter_leaf_components(complist, parent_transform: np.ndarray | None = None):
    """Yield renderable components with transforms composed through subsystems."""
    if parent_transform is None:
        parent_transform = np.eye(4)

    for component, position, rotation in complist:
        transform = parent_transform @ _transformation_matrix(position, rotation)
        if hasattr(component, "complist"):
            yield from _iter_leaf_components(component.complist, transform)
        else:
            yield component, transform


def plot_system_plotly(
    obj: Any,
    title: str = "Optical System (3D View)",
    dark_mode: bool = False,
    width: int = 900,
    height: int = 550,
    camera_eye: dict[str, float] | None = None,
    show_rays: bool = True,
    ray_width: int = 4,
) -> Any:
    """Creates an interactive 3D WebGL visualization of an optical system using Plotly.

    Parameters:
        obj: pyOpTools System, Component, Surface, Ray, or list of Rays to visualize.
        title: Title of the visualization plot.
        dark_mode: If True, uses dark background theme.
        width: Plot width in pixels.
        height: Plot height in pixels.
        camera_eye: Optional dictionary with 'x', 'y', 'z' camera viewpoint coordinates.
        show_rays: Whether to render rays.
        ray_width: Width of the ray lines in pixels.

    Returns:
        plotly.graph_objects.Figure ready to display in Marimo, Jupyter, or export to HTML.
    """
    if go is None:
        raise ImportError(
            "Plotly is required for plot_system_plotly. Install it via 'pip install plotly'."
        )

    traces = []

    # Identify object type: System, Component, Surface, or Rays
    complist = []
    rays_to_render = []

    if hasattr(obj, "complist"):
        # System
        complist = obj.complist
        if show_rays and hasattr(obj, "prop_ray"):
            rays_to_render = obj.prop_ray
    elif hasattr(obj, "surflist"):
        # Single Component
        complist = [(obj, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))]
    elif hasattr(obj, "polylist"):
        # Single Surface wrapped in a dummy component structure
        class _DummyComp:
            def __init__(self, s):
                self.surflist = {"S": (s, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))}

        complist = [(_DummyComp(obj), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))]
    elif isinstance(obj, (list, tuple)) and len(obj) > 0 and hasattr(obj[0], "origin"):
        # List of rays
        if show_rays:
            rays_to_render = obj
    elif hasattr(obj, "origin") and hasattr(obj, "direction"):
        # Single ray
        if show_rays:
            rays_to_render = [obj]

    # 1. Render Optical Components with Hierarchical Transformations
    for component, T_comp in _iter_leaf_components(complist):
        comp_name = type(component).__name__

        surflist = getattr(component, "surflist", [])
        if isinstance(surflist, dict):
            surf_items = list(surflist.items())
        else:
            surf_items = list(enumerate(surflist))

        for surf_key, surf_item in surf_items:
            if isinstance(surf_item, (list, tuple)) and len(surf_item) == 3:
                surf_obj, surf_pos, surf_rot = surf_item
            else:
                surf_obj, surf_pos, surf_rot = (
                    surf_item,
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                )

            if not hasattr(surf_obj, "polylist"):
                continue

            pts, polys = surf_obj.polylist()
            pts = np.array(pts, dtype=float)
            polys = np.array(polys, dtype=int)
            if len(pts) == 0 or len(polys) == 0:
                continue
            color, opacity, smooth, thickness = _resolve_render_style(
                component, surf_obj
            )
            pts, polys = _extrude_planar_mesh(pts, polys, thickness)

            # Composite 4x4 matrix: World = T_comp @ T_surface
            T_surf = _transformation_matrix(surf_pos, surf_rot)
            T_total = T_comp @ T_surf

            pts_4d = np.hstack([pts, np.ones((len(pts), 1), dtype=float)])
            world_pts = (T_total @ pts_4d.T).T[:, :3]

            mesh = go.Mesh3d(
                x=world_pts[:, 0],
                y=world_pts[:, 1],
                z=world_pts[:, 2],
                i=polys[:, 0],
                j=polys[:, 1],
                k=polys[:, 2],
                name=f"{comp_name} ({surf_key})",
                color=color,
                opacity=opacity,
                flatshading=not smooth,
                lighting={
                    "ambient": 0.7,
                    "diffuse": 0.9,
                    "roughness": 0.1,
                    "specular": 0.8,
                    "fresnel": 0.4,
                },
                lightposition={"x": 100, "y": 200, "z": 500},
                hoverinfo="name",
                showlegend=False,
            )
            traces.append(mesh)

    # 2. Render propagated rays grouped by their effective display color.
    ray_groups: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {}
    for ray in rays_to_render:
        for p1, p2, color in _extract_ray_segments(ray):
            ray_groups.setdefault(color, []).append((p1, p2))

    for color, segments in ray_groups.items():
        rx, ry, rz = [], [], []
        for p1, p2 in segments:
            rx.extend([p1[0], p2[0], None])
            ry.extend([p1[1], p2[1], None])
            rz.extend([p1[2], p2[2], None])

        ray_trace = go.Scatter3d(
            x=rx,
            y=ry,
            z=rz,
            mode="lines",
            line={"color": color, "width": ray_width},
            name="Beam",
            hoverinfo="name",
            showlegend=True,
        )
        traces.append(ray_trace)

    # 3. Configure Camera & True 1:1:1 Physical Aspect Ratio
    template = "plotly_dark" if dark_mode else "plotly_white"
    fig = go.Figure(data=traces)
    fig.update_layout(
        title={"text": title, "x": 0.05, "y": 0.95},
        template=template,
        width=width,
        height=height,
        scene={
            "aspectmode": "data",  # Guaranteed true 1:1:1 optical scale
            "xaxis": {
                "title": "X (mm)",
                "gridcolor": "#E5E7EB" if not dark_mode else "#333333",
            },
            "yaxis": {
                "title": "Y (mm)",
                "gridcolor": "#E5E7EB" if not dark_mode else "#333333",
            },
            "zaxis": {
                "title": "Z (mm)",
                "gridcolor": "#E5E7EB" if not dark_mode else "#333333",
            },
            "camera": {
                "eye": camera_eye or {"x": 1.0, "y": 0, "z": 0},
                "up": {"x": 0, "y": -1, "z": 0},
                "center": {"x": 0, "y": 0, "z": 0},
                "projection": {"type": "orthographic"},
            },
        },
        margin={"l": 10, "r": 10, "b": 10, "t": 50},
        legend={"x": 0.02, "y": 0.98},
    )
    return fig
