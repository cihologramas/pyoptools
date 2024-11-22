import pytest
from math import isclose, sin, cos, radians

from pyoptools.raytrace.ray.ray import Ray
from pyoptools.raytrace.system.idealcomponent import IdealThickLens
from pyoptools.raytrace.system.system import System
from pyoptools.raytrace.shape.circular import Circular
from pyoptools.raytrace.calc.calc import intersection


def test_parallel_rays():
    focal_length = 100
    lens_thickness = 10
    a = IdealThickLens(
        lens_shape=Circular(radius=10.001),
        lens_thickness=lens_thickness,
        focal_length=focal_length,
    )
    s = System(complist=[(a, (0, 0, 0), (0, 0, 0))])

    # Test 1: Rays through center and edge
    r1 = Ray((0, 0, -20), (0, 0, 1))
    r2 = Ray((10, 0, -20), (0, 0, 1))
    s.ray_add([r1, r2])
    s.propagate()
    r1f, r2f = r1.get_final_rays()[0], r2.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r1f, r2f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)

    # Test 2: Rays at different Y positions
    r3 = Ray((0, 5, -20), (0, 0, 1))
    r4 = Ray((0, -5, -20), (0, 0, 1))
    s.ray_add([r3, r4])
    s.propagate()
    r3f, r4f = r3.get_final_rays()[0], r4.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r3f, r4f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)

    # Test 3: Rays at diagonal positions
    r5 = Ray((5, 5, -20), (0, 0, 1))
    r6 = Ray((-5, -5, -20), (0, 0, 1))
    s.ray_add([r5, r6])
    s.propagate()
    r5f, r6f = r5.get_final_rays()[0], r6.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r5f, r6f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)

    # Test 4: Rays at different X positions
    r7 = Ray((7, 0, -20), (0, 0, 1))
    r8 = Ray((-7, 0, -20), (0, 0, 1))
    s.ray_add([r7, r8])
    s.propagate()
    r7f, r8f = r7.get_final_rays()[0], r8.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r7f, r8f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)

    # Test 5: Rays near the edge
    r9 = Ray((9, 0, -20), (0, 0, 1))
    r10 = Ray((-9, 0, -20), (0, 0, 1))
    s.ray_add([r9, r10])
    s.propagate()
    r9f, r10f = r9.get_final_rays()[0], r10.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r9f, r10f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)


def test_skewed_rays():
    focal_length = 100
    lens_thickness = 10
    a = IdealThickLens(
        lens_shape=Circular(radius=30.001),
        lens_thickness=lens_thickness,
        focal_length=focal_length,
    )
    s = System(complist=[(a, (0, 0, 0), (0, 0, 0))])

    # Test 1: Parallel rays at 45 degrees in XZ plane
    angle = radians(45)
    r1 = Ray((5, 0, -20), (sin(angle), 0, cos(angle)))
    r2 = Ray((-5, 0, -20), (sin(angle), 0, cos(angle)))
    s.ray_add([r1, r2])
    s.propagate()
    r1f, r2f = r1.get_final_rays()[0], r2.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r1f, r2f)[0]
    print(int_x, int_y, int_z)
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)

    # Test 2: Rays at 30 degrees in YZ plane
    angle = radians(30)
    r3 = Ray((0, 5, -20), (0, sin(angle), cos(angle)))
    r4 = Ray((0, -5, -20), (0, sin(angle), cos(angle)))
    s.ray_add([r3, r4])
    s.propagate()
    r3f, r4f = r3.get_final_rays()[0], r4.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r3f, r4f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)

    # Test 3: Rays at 60 degrees in both XZ and YZ planes
    angle = radians(60)
    r5 = Ray((5, 5, -10), (0, -sin(angle), cos(angle)))
    r6 = Ray((-5, -5, -10), (0, -sin(angle), cos(angle)))
    s.ray_add([r5, r6])
    s.propagate()
    r5f, r6f = r5.get_final_rays()[0], r6.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r5f, r6f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)

    # Test 4: Rays at varying angles
    r7 = Ray((3, -2, -20), (-0.2, 0.1, 0.9))
    r8 = Ray((-3, 2, -20), (-0.2, 0.1, 0.9))
    s.ray_add([r7, r8])
    s.propagate()
    r7f, r8f = r7.get_final_rays()[0], r8.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r7f, r8f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)

    # Test 5: Rays at shallow angle
    angle = radians(1)
    r9 = Ray((9, 0, -20), (-sin(angle), 0, cos(angle)))
    r10 = Ray((-9, 0, -20), (-sin(angle), 0, cos(angle)))
    s.ray_add([r9, r10])
    s.propagate()
    r9f, r10f = r9.get_final_rays()[0], r10.get_final_rays()[0]
    int_x, int_y, int_z = intersection(r9f, r10f)[0]
    assert isclose(int_z, focal_length + lens_thickness / 2, rel_tol=1e-9)
