import pytest

from pyoptools.gui.ipywidgets import Plot3D
from pyoptools.gui.plotly_viewer import _wavelength_to_color, plot_system_plotly
from pyoptools.raytrace.comp_lib import (
    Block,
    CylindricalLens,
)
from pyoptools.raytrace.ray import Ray
from pyoptools.raytrace.shape import Rectangular
from pyoptools.raytrace.surface import Plane
from pyoptools.raytrace.system import System

plotly = pytest.importorskip("plotly")
import plotly.graph_objects as go


def test_wavelength_to_color():
    # 405 nm laser diode
    assert _wavelength_to_color(0.405) == "#2563EB"
    assert _wavelength_to_color(405.0) == "#2563EB"
    # 532 nm green
    assert _wavelength_to_color(0.532) == "#16A34A"
    # 633 nm red
    assert _wavelength_to_color(0.633) == "#DC2626"


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
    assert any("Beam (405 nm)" in name for name in trace_names)
    assert any("CylindricalLens" in name for name in trace_names)


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
