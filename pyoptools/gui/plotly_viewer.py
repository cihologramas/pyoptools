"""Interactive 3D visualization for pyOpTools optical systems using Plotly.

Provides a modern, universal WebGL 3D viewer compatible with Marimo,
JupyterLab 4, VS Code Interactive Window, and standalone HTML exports.
"""

from __future__ import annotations

from typing import Any

import numpy as np

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


def _extract_ray_segments(ray) -> list[tuple[np.ndarray, np.ndarray, float]]:
    """Recursively extracts line segments and wavelengths from a pyOpTools ray and its children."""
    segments = []
    p1 = np.array(ray.origin, dtype=float)
    if len(ray.childs) > 0:
        p2 = np.array(ray.childs[0].origin, dtype=float)
    else:
        p2 = p1 + 15.0 * np.array(ray.direction, dtype=float)

    wl = float(getattr(ray, "wavelength", 0.589))
    if getattr(ray, "intensity", 1.0) != 0:
        segments.append((p1, p2, wl))

    for child in ray.childs:
        segments.extend(_extract_ray_segments(child))
    return segments


def _wavelength_to_color(wl: float) -> str:
    """Map optical wavelength to an RGB hex color with accurate laser beam rendering.

    Wavelength can be provided in micrometers (e.g. 0.405) or nanometers (e.g. 405.0).
    405 nm GaN laser diode light is rendered as crisp Royal Blue.
    """
    wl_nm = wl * 1000.0 if wl < 10.0 else wl
    if 380 <= wl_nm <= 440:
        return "#2563EB"  # 405 nm Royal Blue Laser Diode
    elif 440 < wl_nm <= 495:
        return "#0284C7"  # Cyan / Sky blue (488 nm Argon)
    elif 495 < wl_nm <= 570:
        return "#16A34A"  # Emerald Green (532 nm DPSS)
    elif 570 < wl_nm <= 590:
        return "#EAB308"  # Sodium Yellow (589 nm)
    elif 590 < wl_nm <= 620:
        return "#EA580C"  # Orange (604 nm)
    elif 620 < wl_nm <= 750:
        return "#DC2626"  # Helium-Neon Red (633 nm)
    elif wl_nm > 750:
        return "#7F1D1D"  # Near-Infrared
    return "#2563EB"


def _generate_lens_side_mesh(
    component, T_comp: np.ndarray
) -> tuple[np.ndarray, np.ndarray] | tuple[None, None]:
    """Generate watertight side walls connecting S1 and S2 perimeter points."""
    surflist = getattr(component, "surflist", {})
    if isinstance(surflist, (list, tuple)):
        surflist = dict(enumerate(surflist))

    s1_item = surflist.get("S1") or surflist.get(0)
    s2_item = surflist.get("S2") or surflist.get(1)
    if not s1_item or not s2_item:
        return None, None

    surf1, pos1, rot1 = (
        s1_item
        if isinstance(s1_item, (list, tuple)) and len(s1_item) == 3
        else (s1_item, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    )
    surf2, pos2, rot2 = (
        s2_item
        if isinstance(s2_item, (list, tuple)) and len(s2_item) == 3
        else (s2_item, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    )

    if not hasattr(surf1, "polylist") or not hasattr(surf2, "polylist"):
        return None, None

    pts1, _ = surf1.polylist()
    pts2, _ = surf2.polylist()
    if len(pts1) == 0 or len(pts2) == 0:
        return None, None

    T_surf1 = _transformation_matrix(pos1, rot1)
    T_surf2 = _transformation_matrix(pos2, rot2)
    T_tot1 = T_comp @ T_surf1
    T_tot2 = T_comp @ T_surf2

    pts1_4d = np.hstack([pts1, np.ones((len(pts1), 1), dtype=float)])
    pts2_4d = np.hstack([pts2, np.ones((len(pts2), 1), dtype=float)])

    world1 = (T_tot1 @ pts1_4d.T).T[:, :3]
    world2 = (T_tot2 @ pts2_4d.T).T[:, :3]

    nx, ny = getattr(getattr(surf1, "shape", None), "samples", (30, 30))
    if len(world1) != nx * ny or len(world2) != nx * ny:
        return None, None

    grid1 = world1.reshape((ny, nx, 3))
    grid2 = world2.reshape((ny, nx, 3))

    perimeter_coords = []
    for c in range(nx - 1):
        perimeter_coords.append((0, c))
    for r in range(ny - 1):
        perimeter_coords.append((r, nx - 1))
    for c in range(nx - 1, 0, -1):
        perimeter_coords.append((ny - 1, c))
    for r in range(ny - 1, 0, -1):
        perimeter_coords.append((r, 0))

    N = len(perimeter_coords)
    side_pts = np.empty((2 * N, 3), dtype=float)
    for idx, (r, c) in enumerate(perimeter_coords):
        side_pts[idx] = grid1[r, c]
        side_pts[N + idx] = grid2[r, c]

    side_polys = np.empty((2 * N, 3), dtype=int)
    for k in range(N):
        k_next = (k + 1) % N
        p1, p2, p3, p4 = k, k_next, N + k_next, N + k
        side_polys[2 * k] = [p1, p2, p3]
        side_polys[2 * k + 1] = [p1, p3, p4]

    return side_pts, side_polys


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
        rays_to_render = obj
    elif hasattr(obj, "origin") and hasattr(obj, "direction"):
        # Single ray
        rays_to_render = [obj]

    # 1. Render Optical Components with Hierarchical Transformations
    for comp in complist:
        component, comp_pos, comp_rot = comp
        T_comp = _transformation_matrix(comp_pos, comp_rot)
        comp_name = type(component).__name__

        # Material styling: clean, distinct optical materials
        if "Polygon" in comp_name or "Prism" in comp_name:
            color = "#38BDF8"  # Refractive optical glass
            opacity = 0.40
            flatshading = True
        elif (
            "Cylindrical" in comp_name
            or "Lens" in comp_name
            or "Spherical" in comp_name
        ):
            color = "#BAE6FD"  # High-transparency optical glass
            opacity = 0.45
            flatshading = False
        elif "Mirror" in comp_name:
            color = "#E2E8F0"  # Silver/aluminum mirror
            opacity = 0.90
            flatshading = True
        elif "CCD" in comp_name or "Detector" in comp_name or "pd" in str(comp).lower():
            color = "#475569"  # Silicon detector sensor
            opacity = 0.90
            flatshading = True
        elif "Stop" in comp_name or "Aperture" in comp_name:
            color = "#1E293B"  # Matte dark chassis
            opacity = 0.95
            flatshading = True
        else:
            color = "#CBD5E1"
            opacity = 0.50
            flatshading = True

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

            # Composite 4x4 matrix: World = T_comp @ T_surface
            T_surf = _transformation_matrix(surf_pos, surf_rot)
            T_total = T_comp @ T_surf

            pts_4d = np.hstack([pts, np.ones((len(pts), 1), dtype=float)])
            world_pts = (T_total @ pts_4d.T).T[:, :3]

            facet_color = (
                "#F59E0B" if getattr(surf_obj, "reflectivity", 0) == 1 else color
            )
            surf_opacity = opacity

            is_planar = "Plane" in type(surf_obj).__name__
            mesh = go.Mesh3d(
                x=world_pts[:, 0],
                y=world_pts[:, 1],
                z=world_pts[:, 2],
                i=polys[:, 0],
                j=polys[:, 1],
                k=polys[:, 2],
                name=f"{comp_name} ({surf_key})",
                color=facet_color,
                opacity=surf_opacity,
                flatshading=True if is_planar else flatshading,
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

    # 2. Render Propagated Rays grouped by wavelength
    ray_groups: dict[float, list[tuple[np.ndarray, np.ndarray]]] = {}
    for ray in rays_to_render:
        for p1, p2, wl in _extract_ray_segments(ray):
            ray_groups.setdefault(round(wl, 5), []).append((p1, p2))

    for wl, segments in ray_groups.items():
        rx, ry, rz = [], [], []
        for p1, p2 in segments:
            rx.extend([p1[0], p2[0], None])
            ry.extend([p1[1], p2[1], None])
            rz.extend([p1[2], p2[2], None])

        wl_nm = int(wl * 1000.0 if wl < 10.0 else wl)
        beam_color = _wavelength_to_color(wl)
        ray_trace = go.Scatter3d(
            x=rx,
            y=ry,
            z=rz,
            mode="lines",
            line={"color": beam_color, "width": ray_width},
            name=f"Beam ({wl_nm} nm)",
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
                "eye": camera_eye or {"x": -1.6, "y": -1.6, "z": 1.3},
                "up": {"x": 0, "y": 0, "z": 1},
            },
        },
        margin={"l": 10, "r": 10, "b": 10, "t": 50},
        legend={"x": 0.02, "y": 0.98},
    )
    return fig
