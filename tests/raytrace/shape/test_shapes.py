import pickle

import pytest
from pyoptools.raytrace.shape.circular import Circular
from pyoptools.raytrace.shape.polygon import Polygon


def test_circular():
    c = Circular()
    assert c.radius == 1.0, "Unexpected radius"


def test_polygon_validation():
    with pytest.raises(ValueError):
        Polygon(coord=((0, 0), (1, 1)))  # Less than 3 vertices


def test_polygon_triangle():
    poly = Polygon(coord=((0, 0), (10, 0), (0, 10)))
    assert poly.coord == ((0.0, 0.0), (10.0, 0.0), (0.0, 10.0))
    assert poly.limits() == (0.0, 10.0, 0.0, 10.0)
    assert poly.hit((2, 2, 0)) is True
    assert poly.hit((8, 8, 0)) is False
    assert poly.hit((-1, 2, 0)) is False


def test_polygon_trapezoid():
    # Isosceles trapezoid
    poly = Polygon(coord=((-10, -5), (10, -5), (5, 5), (-5, 5)))
    assert poly.limits() == (-10.0, 10.0, -5.0, 5.0)
    assert poly.hit((0, 0, 0)) is True
    assert poly.hit((0, -4.5, 0)) is True
    assert poly.hit((0, 4.5, 0)) is True
    assert poly.hit((8, 3, 0)) is False  # outside slanted side


def test_polygon_pentagon():
    poly = Polygon(coord=((10, -10), (-10, -10), (-18, 10), (-10, 18), (10, 10)))
    assert poly.limits() == (-18.0, 10.0, -10.0, 18.0)
    assert poly.hit((0, 0, 0)) is True
    assert poly.hit((-12, 12, 0)) is True
    assert poly.hit((15, 0, 0)) is False


def test_polygon_pointlist():
    poly = Polygon(coord=((0, 0), (10, 0), (10, 10), (0, 10)), samples=5)
    X, Y = poly.pointlist()
    assert len(X) == len(Y)
    assert len(X) >= 4  # Includes vertices plus interior samples
    for x, y in zip(X, Y):
        assert poly.hit((x, y, 0)) or (x in (0, 10) and y in (0, 10))


def test_polygon_pickling():
    poly = Polygon(coord=((-5, -5), (5, -5), (0, 5)), samples=15)
    data = pickle.dumps(poly)
    restored = pickle.loads(data)
    assert restored.coord == poly.coord
    assert restored.samples == poly.samples
    assert restored.hit((0, 0, 0)) is True
