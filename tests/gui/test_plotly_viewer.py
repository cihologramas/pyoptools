import pickle

import numpy as np
import pytest

from pyoptools.gui.ipywidgets import Plot3D
from pyoptools.gui.plotly_viewer import _extrude_planar_mesh, plot_system_plotly
from pyoptools.raytrace.comp_lib import (
    CCD,
    Block,
    CylindricalLens,
    Doublet,
    RectMirror,
    RoundMirror,
    Stop,
)
from pyoptools.raytrace.component import Component
from pyoptools.raytrace.ray import Ray
from pyoptools.raytrace.shape import Rectangular
from pyoptools.raytrace.surface import ArrayDetector, Plane
from pyoptools.raytrace.system import System

plotly = pytest.importorskip("plotly")
import plotly.graph_objects as go


def test_plot_system_plotly_uses_ray_draw_color():
    ray = Ray(origin=(0, 0, 0), direction=(0, 0, 1), draw_color="rebeccapurple")

    fig = plot_system_plotly(ray)

    assert len(fig.data) == 1
    assert fig.data[0].line.color == "#663399"


def test_plot_system_plotly_hides_all_ray_inputs_when_disabled():
    ray = Ray()

    assert len(plot_system_plotly(ray, show_rays=False).data) == 0
    assert len(plot_system_plotly([ray], show_rays=False).data) == 0


def test_plot_system_plotly_uses_optical_axis_camera():
    fig = plot_system_plotly(Ray())

    camera = fig.layout.scene.camera
    assert camera.eye.to_plotly_json() == {"x": 1.0, "y": 0, "z": 0}
    assert camera.up.to_plotly_json() == {"x": 0, "y": -1, "z": 0}
    assert camera.center.to_plotly_json() == {"x": 0, "y": 0, "z": 0}
    assert camera.projection.type == "orthographic"


def test_plot_system_plotly():
    lens = CylindricalLens(
        size=(20, 20),
        thickness=6,
        curvature_s1=1.0 / 50.0,
        curvature_s2=0.0,
        material=1.5,
    )
    sys = System(complist=[(lens, (0, 0, 20), (0, 0, 0))])
    ray = Ray(origin=(0, 0, 0), direction=(0, 0, 1), wavelength=0.405)
    sys.ray_add(ray)
    sys.propagate()

    fig = plot_system_plotly(sys, title="Test System")
    assert isinstance(fig, go.Figure)
    assert len(fig.data) >= 7  # 6 lens surfaces + 1 ray trace
    trace_names = [t.name for t in fig.data]
    assert "Beam" in trace_names
    assert any("CylindricalLens" in name for name in trace_names)


def test_plot_system_plotly_renders_nested_doublet():
    doublet = Doublet(
        radius=5,
        curvature_s1=1.0 / 50.0,
        curvature_s2=-1.0 / 50.0,
        curvature_s3=-1.0 / 75.0,
        thickness_l1=3,
        thickness_l2=2,
        material_l1=1.5,
        material_l2=1.6,
    )
    sys = System(complist=[(doublet, (100, 20, 30), (0, np.pi / 2, 0))])

    fig = plot_system_plotly(sys)

    lens_traces = [trace for trace in fig.data if "SphericalLens" in trace.name]
    assert len(lens_traces) == 8
    assert min(np.min(trace.x) for trace in lens_traces) > 80
    assert min(np.min(trace.y) for trace in lens_traces) > 10


def test_plot3d_auto_and_component():
    cube = Block(size=(10, 10, 10), material=1.5)
    fig = Plot3D(cube, backend="auto")
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 6  # 6 cube faces


def test_plot_single_surface():
    surf = Plane(shape=Rectangular(size=(15, 15)))
    fig = plot_system_plotly(surf)
    assert isinstance(fig, go.Figure)
    assert len(fig.data) == 1


def test_plot_system_plotly_uses_component_and_surface_styles():
    block = Block(
        size=(10, 10, 10),
        material=1.5,
        render_color="#112233",
        render_opacity=0.25,
        render_smooth=True,
    )
    block.surflist["S1"][0].render_color = "#445566"
    block.surflist["S1"][0].render_opacity = 0.75
    block.surflist["S1"][0].render_smooth = False
    block.surflist["S1"][0].render_thickness = 0.4

    fig = plot_system_plotly(block)

    traces = {trace.name: trace for trace in fig.data}
    assert traces["Block (S1)"].color == "#445566"
    assert traces["Block (S1)"].opacity == 0.75
    assert traces["Block (S1)"].flatshading is True
    assert np.ptp(traces["Block (S1)"].z) == pytest.approx(0.4)
    assert traces["Block (S2)"].color == "#112233"
    assert traces["Block (S2)"].opacity == 0.25
    assert traces["Block (S2)"].flatshading is False


def test_plot_single_surface_uses_surface_style():
    detector = ArrayDetector(
        size=(15, 15),
        render_color="#102030",
        render_opacity=0.6,
        render_smooth=True,
    )

    fig = plot_system_plotly(detector)

    assert fig.data[0].color == "#102030"
    assert fig.data[0].opacity == 0.6
    assert fig.data[0].flatshading is False


def test_render_thickness_defaults_and_component_inheritance():
    detector = ArrayDetector(size=(15, 15))

    assert np.ptp(plot_system_plotly(detector).data[0].z) == pytest.approx(0.1)
    assert np.ptp(plot_system_plotly(CCD()).data[0].z) == pytest.approx(0.1)
    assert np.ptp(plot_system_plotly(Stop()).data[0].z) == pytest.approx(0.1)


def test_render_thickness_extrusion_and_zero_thickness_compatibility():
    surface = Plane(shape=Rectangular(size=(10, 10)))
    points, triangles = surface.polylist()
    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=int)

    slab_points, slab_triangles = _extrude_planar_mesh(points, triangles, 0.2)
    unchanged_points, unchanged_triangles = _extrude_planar_mesh(points, triangles, 0.0)

    assert len(slab_points) == 2 * len(points)
    assert len(slab_triangles) > 2 * len(triangles)
    assert np.ptp(slab_points[:, 2]) == pytest.approx(0.2)
    assert np.allclose(unchanged_points, points)
    assert np.array_equal(unchanged_triangles, triangles)


def test_render_thickness_follows_surface_transform():
    surface = Plane(shape=Rectangular(size=(10, 10)))
    component = Component(
        surflist=[(surface, (0, 0, 0), (np.pi / 2, 0, 0))],
        render_thickness=0.2,
    )
    system = System(complist=[(component, (1, 2, 3), (0, 0, 0))])

    trace = plot_system_plotly(system).data[0]

    assert np.ptp(trace.y) == pytest.approx(0.2)
    assert np.mean(trace.y) == pytest.approx(2.0)


def test_mirror_surface_style_is_independent_of_reflectivity():
    mirror = RoundMirror(reflectivity=0)

    fig = plot_system_plotly(mirror)

    front_surface = next(
        trace for trace in fig.data if trace.name == "RoundMirror (S1)"
    )
    assert front_surface.color == "#8b03a6"
    assert front_surface.opacity == 1.0


def test_mirror_surface_render_properties_update_the_reflective_surface():
    for mirror in (RoundMirror(), RectMirror()):
        assert mirror.surface_render_color == "#8b03a6"
        assert mirror.surface_render_opacity == 1.0

        mirror.surface_render_color = "#123456"
        mirror.surface_render_opacity = 0.75

        assert mirror.surflist["S1"][0].render_color == "#123456"
        assert mirror.surflist["S1"][0].render_opacity == 0.75


def test_mirror_constructor_sets_reflective_surface_render_properties():
    for mirror in (
        RoundMirror(surface_render_color="#123456", surface_render_opacity=0.75),
        RectMirror(surface_render_color="#123456", surface_render_opacity=0.75),
    ):
        assert mirror.surface_render_color == "#123456"
        assert mirror.surface_render_opacity == 0.75


def test_render_style_survives_pickle_round_trip():
    component = Component(
        render_color="#123456",
        render_opacity=0.4,
        render_smooth=True,
        render_thickness=0.3,
    )
    surface = Plane(
        render_color="#654321",
        render_opacity=0.6,
        render_smooth=True,
        render_thickness=0.2,
    )

    restored_component = pickle.loads(pickle.dumps(component))
    restored_surface = pickle.loads(pickle.dumps(surface))

    assert restored_component.render_color == "#123456"
    assert restored_component.render_opacity == 0.4
    assert restored_component.render_smooth is True
    assert restored_component.render_thickness == 0.3
    assert restored_surface.render_color == "#654321"
    assert restored_surface.render_opacity == 0.6
    assert restored_surface.render_smooth is True
    assert restored_surface.render_thickness == 0.2
