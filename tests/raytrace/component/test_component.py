# standard imports

# third-party imports

# local imports
from pyoptools.raytrace import component, shape, surface
from pyoptools.raytrace.mat_lib import material


def test_component_surflist():
    S0 = surface.Spherical(shape=shape.Circular(radius=50), curvature=1.0 / 200.0)
    S1 = surface.Spherical(shape=shape.Circular(radius=50), curvature=1.0 / 200.0)
    S2 = surface.Cylinder(radius=50, length=10)
    L1 = component.Component(
        surflist=[
            (S0, (0, 0, -5), (0, 0, 0)),
            (S1, (0, 0, 5), (0, 0, 0)),
            (S2, (0, 0, 6.5), (0, 0, 0)),
        ],
        material=material.schott["BK7"],
    )
    print(type(L1.surflist))
    for surf in L1.surflist:
        assert len(surf) == 3, "Wrong surface length"


def test_cylindrical_lens():
    import pytest

    from pyoptools.raytrace.comp_lib import CylindricalLens
    from pyoptools.raytrace.ray import Ray

    # Plano-convex cylindrical lens
    lens = CylindricalLens(
        size=(20, 20),
        thickness=10,
        curvature_s1=1.0 / 50.0,
        curvature_s2=0.0,
        material=1.5,
    )
    surfs = lens.surflist
    assert len(surfs) == 6, f"Expected 6 surfaces (S1..S6), got {len(surfs)}"
    for key in ["S1", "S2", "S3", "S4", "S5", "S6"]:
        assert key in surfs, f"Missing surface {key}"

    # Ray hitting through front surface S1
    r_front = Ray(origin=(0, 0, -20), direction=(0, 0, 1))
    res_front = lens.propagate(r_front, 1.0)
    assert len(res_front) > 0

    # Ray hitting lateral surface S5
    r_side = Ray(origin=(-20, 0, 0.5), direction=(1, 0, 0))
    res_side = lens.propagate(r_side, 1.0)
    assert len(res_side) > 0

    # Check non-physical curvature error
    with pytest.raises(ValueError):
        CylindricalLens(size=(50, 20), thickness=10, curvature_s1=1.0 / 20.0)


def test_pentaprism():
    import numpy as np

    from pyoptools.raytrace.comp_lib import PentaPrism
    from pyoptools.raytrace.ray import Ray
    from pyoptools.raytrace.system import System

    penta = PentaPrism(s=20, material=1.5)
    assert len(penta.surflist) == 7, (
        f"Expected 7 surfaces (S1..S7), got {len(penta.surflist)}"
    )
    for key in ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]:
        assert key in penta.surflist, f"Missing surface {key}"

    # Ray entering entrance face S1 and exiting 90 degrees turned through S2
    sys_p = System(complist=[(penta, (0, 0, 0), (0, 0, 0))])
    r_in = Ray(origin=(0, 0, -20), direction=(0, 0, 1), wavelength=0.589)
    sys_p.ray_add(r_in)
    sys_p.propagate()
    final_rays = r_in.get_final_rays()
    assert len(final_rays) == 1
    # Should exit pointing along +X
    np.testing.assert_allclose(final_rays[0].direction, [1.0, 0.0, 0.0], atol=1e-5)

    # Ray hitting top closing surface S6 (y = 10)
    r_top = Ray(origin=(0, 20, 0), direction=(0, -1, 0))
    res_top = penta.propagate(r_top, 1.0)
    assert len(res_top) > 0


def test_doveprism():

    from pyoptools.raytrace.comp_lib import DovePrism
    from pyoptools.raytrace.ray import Ray

    dove = DovePrism(s=15, length=60, material=1.5)
    assert len(dove.surflist) == 6, (
        f"Expected 6 surfaces (S1..S6), got {len(dove.surflist)}"
    )
    for key in ["S1", "S2", "S3", "S4", "S5", "S6"]:
        assert key in dove.surflist, f"Missing surface {key}"

    # Ray hitting top surface S5 (y = 7.5)
    r_top = Ray(origin=(0, 20, 0), direction=(0, -1, 0))
    res_top = dove.propagate(r_top, 1.0)
    assert len(res_top) > 0


def test_component_hit_list():
    """Verify that Component.hit_list aggregates surface hits in component coordinates."""
    import numpy as np

    from pyoptools.raytrace.comp_lib import SphericalLens
    from pyoptools.raytrace.ray import Ray
    from pyoptools.raytrace.system import System

    lens = SphericalLens(
        radius=25.0, thickness=5.0, curvature_s1=1.0 / 50.0, curvature_s2=-1.0 / 50.0
    )
    s = System(complist=[(lens, (0, 0, 30), (0, 0, 0))])
    r = Ray(origin=(0, 0, 0), direction=(0, 0, 1))
    s.propagate_ray(r)

    hl = lens.hit_list
    assert len(hl) == 2, f"Expected 2 hits (front & back surface), got {len(hl)}"
    p0, r0 = hl[0]
    p1, r1 = hl[1]
    # In component coordinates, front surface vertex is at z = -2.5, rear at z = +2.5
    np.testing.assert_allclose(p0, [0.0, 0.0, -2.5], atol=1e-5)
    np.testing.assert_allclose(p1, [0.0, 0.0, 2.5], atol=1e-5)
