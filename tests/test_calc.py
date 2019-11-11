# standard imports

# third-party imports
import numpy as np
import pytest

# local imports
from pyoptools.raytrace._comp_lib.ccd import CCD
from pyoptools.raytrace._comp_lib.spherical_lens import SphericalLens
from pyoptools.raytrace._comp_lib.stop import Stop
import pyoptools.raytrace.calc.calc as calc
from pyoptools.raytrace.mat_lib import material
from pyoptools.raytrace.ray import Ray
from pyoptools.raytrace.shape.circular import Circular
from pyoptools.raytrace.system.system import System


def test_intersection():
    # TODO: choose better example
    # Intersecting
    pos = (0, 0, 0)
    r1 = Ray(pos=pos, dir=(0, 0.2, 1))
    r2 = Ray(pos=pos, dir=(1, 2, 1))
    ip, rv = calc.intersection(r1, r2)
    np.testing.assert_almost_equal(ip, pos)
    assert(rv == True)

    # Not intersecting
    r1 = Ray(pos=(0, 0, 0), dir=(0, 0.2, 1))
    r2 = Ray(pos=(1, 0, 0), dir=(0, 0.2, 1))
    ip, rv = calc.intersection(r1, r2)
    np.testing.assert_almost_equal(ip, [np.nan, np.nan, np.nan])
    assert(rv == False)


def test_nearest_points():
    # Real closest point
    r1 = Ray(pos=(0, 0, 0), dir=(1, 1, 0))
    r2 = Ray(pos=(1, 0, 1), dir=(0, 1, 0))
    p1, p2, d, rv = calc.nearest_points(r1, r2)
    np.testing.assert_almost_equal(p1, [1.0, 1.0, 0.0])
    np.testing.assert_almost_equal(p2, [1.0, 1.0, 1.0])
    assert(d == 1)
    assert(rv == True)

    # Virtual closest point
    r1 = Ray(pos=(0, 0, 0), dir=(1, 1, 0))
    r2 = Ray(pos=(1, 10, 1), dir=(0, 1, 0))
    p1, p2, d, rv = calc.nearest_points(r1, r2)
    np.testing.assert_almost_equal(p1, [1.0, 1.0, 0.0])
    np.testing.assert_almost_equal(p2, [1.0, 1.0, 1.0])
    assert(d == 1)
    assert(rv == False)


def test_chief_ray_search():
    l1 = SphericalLens(
        radius=25,
        curvature_s1=1.0 / 100.0,
        curvature_s2=-1.0 / 100,
        thickness=10,
        material=material.schott["BK7"],
    )
    l2 = SphericalLens(
        radius=25,
        curvature_s1=1.0 / 100.0,
        curvature_s2=-1.0 / 100,
        thickness=10,
        material=material.schott["BK7"],
    )
    c = CCD()
    ap = Stop(shape=Circular(radius=30), ap_shape=Circular(radius=25))
    s = System(
        complist=[
            (l1, (0, 0, 100), (0, 0, 0)),
            (l2, (0, 0, 120), (0, 0, 0)),
            (ap, (0, 0, 110), (0, 0, 0)),
            (c, (0, 0, 150), (0, 0, 0)),
        ],
        n=1,
    )
    chief_ray = calc.chief_ray_search(s, ap, (0, 10, 0), (0, -1, 1))

    # TODO: should implement Ray.__eq__() and use it here
    np.testing.assert_almost_equal(chief_ray.pos, [0, 10, 0])
    np.testing.assert_almost_equal(
        chief_ray.dir, [3.58848263e-04, -9.26093228e-02, 9.95702458e-01], decimal=3
    )
    np.testing.assert_almost_equal(chief_ray.intensity, 1)
    np.testing.assert_almost_equal(chief_ray.wavelength, 0.58929)
    np.testing.assert_almost_equal(chief_ray.n, 1)
    np.testing.assert_equal(chief_ray.label, "")
    np.testing.assert_equal(chief_ray.orig_surf, None)
    np.testing.assert_almost_equal(chief_ray.order, 0)