# standard imports

# third-party imports
from nose.tools import assert_equal
import numpy as np

# local imports
import pyoptools.raytrace.calc.calc as calc
from pyoptools.raytrace.ray import Ray


def test_intersection():
    # TODO: choose better example
    # Intersecting
    pos = (0, 0, 0)
    r1 = Ray(pos=pos, dir=(0, 0.2, 1))
    r2 = Ray(pos=pos, dir=(1, 2, 1))
    ip, rv = calc.intersection(r1, r2)
    np.testing.assert_almost_equal(ip, pos)
    assert_equal(rv, True)

    # Not intersecting
    r1 = Ray(pos=(0, 0, 0), dir=(0, 0.2, 1))
    r2 = Ray(pos=(1, 0, 0), dir=(0, 0.2, 1))
    ip, rv = calc.intersection(r1, r2)
    np.testing.assert_almost_equal(ip, [np.nan, np.nan, np.nan])
    assert_equal(rv, False)


def test_nearest_points():
    # Real closest point
    r1 = Ray(pos=(0, 0, 0), dir=(1, 1, 0))
    r2 = Ray(pos=(1, 0, 1), dir=(0, 1, 0))
    p1, p2, d, rv = calc.nearest_points(r1, r2)
    np.testing.assert_almost_equal(p1, [1., 1., 0.])
    np.testing.assert_almost_equal(p2, [1., 1., 1.])
    assert_equal(d, 1)
    assert_equal(rv, True)

    # Virtual closest point
    r1 = Ray(pos=(0, 0, 0), dir=(1, 1, 0))
    r2 = Ray(pos=(1, 10, 1), dir=(0, 1, 0))
    p1, p2, d, rv = calc.nearest_points(r1, r2)
    np.testing.assert_almost_equal(p1, [1., 1., 0.])
    np.testing.assert_almost_equal(p2, [1., 1., 1.])
    assert_equal(d, 1)
    assert_equal(rv, False)


if __name__ == "__main__":
    import nose

    nose.runmodule()
