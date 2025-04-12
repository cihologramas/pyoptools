from math import isnan

from pyoptools.raytrace.ray.ray import Ray
from pyoptools.raytrace.shape.circular import Circular

# local imports
from pyoptools.raytrace.surface.spherical import Spherical


def test_surface_ray_intersection_backpropagation():
    """Test ray-surface intersection calculation with rays propagating in
    different directions.

    This test verifies two cases:
    1. When a ray propagates in a direction that doesn't intersect the surface
       (forward +z), the intersection point should contain NaN values.
    2. When a ray propagates toward the surface (backward -z),
       a valid intersection point should be calculated.
    """
    # Create a test spherical surface with 10mm radius of curvature and 10mm
    # aperture
    s = Spherical(curvature=1 / 10.0, shape=Circular(radius=10.0))

    # Case 1: Ray pointing away from the surface (in +z direction)
    r = Ray([0.0, 0.0, 10.0], [0.0, 0.0, 1.0])
    intersection_point = s.intersection(r)

    # Verify no intersection occurs (values should be NaN)
    assert all(isnan(coord) for coord in intersection_point)

    # Case 2: Ray pointing toward the surface (in -z direction)
    r = Ray([0.0, 0.0, 10.0], [0.0, 0.0, -1.0])
    intersection_point = s.intersection(r)

    # Verify a valid intersection is found (no NaN values)
    assert all(not isnan(coord) for coord in intersection_point)
