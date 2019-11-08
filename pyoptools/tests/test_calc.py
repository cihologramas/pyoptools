from nose.tools import assert_equal
from numpy import array, nan
import numpy as np

import pyoptools.raytrace.calc.calc as calc
from pyoptools.raytrace.ray import Ray


def test_intersection():
    # Intersecting
    pos = (0, 0, 0)
    r1 = Ray(pos=pos, dir=(0, 0.2, 1))
    r2 = Ray(pos=pos, dir=(1, 2, 1))
    ip, rv = calc.intersection(r1, r2)
    np.testing.assert_array_equal(ip, pos)
    assert_equal(rv, True)

    # Not intersecting
    r1 = Ray(pos=(0, 0, 0), dir=(0, 0.2, 1))
    r2 = Ray(pos=(1, 0, 0), dir=(0, 0.2, 1))
    ip, rv = calc.intersection(r1, r2)
    np.testing.assert_array_equal(ip, [nan, nan, np.nan])
    assert_equal(rv, False)


if __name__ == "__main__":
    import nose

    nose.runmodule()
