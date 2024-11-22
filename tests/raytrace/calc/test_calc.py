import pytest
import numpy as np
from itertools import permutations

from pyoptools.raytrace._comp_lib.ccd import CCD
from pyoptools.raytrace._comp_lib.spherical_lens import SphericalLens
from pyoptools.raytrace._comp_lib.stop import Stop
import pyoptools.raytrace.calc.calc as calc
from pyoptools.raytrace.library import library
from pyoptools.raytrace.mat_lib import material
from pyoptools.raytrace.ray import Ray, parallel_beam_c
from pyoptools.raytrace.shape.circular import Circular
from pyoptools.raytrace.system.system import System


def test_intersection():
    # Intersecting exactly
    expected_intersection_point = (3, 1, 1)
    d1 = np.array((0.1, 0.2, 0.3))
    d2 = np.array((1, 0.5, 0))
    t1 = 2
    t2 = 3
    ray1 = Ray(origin=expected_intersection_point - t1 * d1, direction=d1)
    ray2 = Ray(origin=expected_intersection_point - t2 * d2, direction=d2)
    intersection_point, real_ = calc.intersection(ray1, ray2)
    np.testing.assert_almost_equal(intersection_point, expected_intersection_point)
    assert real_

    # Approximately intersecting
    ray1 = Ray(
        origin=np.array((-1.18414888e-15, 0.0, 5.603)),
        direction=np.array((-5.22664407e-18, 0.0, 1.0)),
    )
    ray2 = Ray(
        origin=np.array((0.01068222, 0.0, 5.60299917)),
        direction=np.array((-5.56560819e-05, 0.0, 9.99999998e-01)),
    )
    intersection_point, real_ = calc.intersection(ray1, ray2)
    expected_intersection_point = np.array((0.0, 0.0, 1.97535672e02))
    np.testing.assert_almost_equal(intersection_point, expected_intersection_point)
    assert real_

    # Not intersecting
    ray1 = Ray(origin=(0, 0, 0), direction=(0, 0.2, 1))
    ray2 = Ray(origin=(1, 0, 0), direction=(0, 0.2, 1))
    intersection_point, real_ = calc.intersection(ray1, ray2)
    np.testing.assert_almost_equal(intersection_point, [np.nan, np.nan, np.nan])
    assert not real_


def test_nearest_points():
    # Real closest point
    ray1 = Ray(origin=(0, 0, 0), direction=(1, 1, 0))
    ray2 = Ray(origin=(1, 0, 1), direction=(0, 1, 0))
    closest_point_on_ray1, closest_point_on_ray2, distance, real_ = calc.nearest_points(
        ray1, ray2
    )
    np.testing.assert_almost_equal(closest_point_on_ray1, [1.0, 1.0, 0.0])
    np.testing.assert_almost_equal(closest_point_on_ray2, [1.0, 1.0, 1.0])
    assert distance == 1
    assert real_

    # Virtual closest point
    ray1 = Ray(origin=(0, 0, 0), direction=(1, 1, 0))
    ray2 = Ray(origin=(1, 10, 1), direction=(0, 1, 0))
    closest_point_on_ray1, closest_point_on_ray2, distance, real_ = calc.nearest_points(
        ray1, ray2
    )
    np.testing.assert_almost_equal(closest_point_on_ray1, [1.0, 1.0, 0.0])
    np.testing.assert_almost_equal(closest_point_on_ray2, [1.0, 1.0, 1.0])
    assert distance == 1
    assert not real_


@pytest.mark.skip(
    reason="This test is failing, need to check the chief ray search implementation"
)
def test_chief_ray_search():
    lens1 = SphericalLens(
        radius=25,
        curvature_s1=1.0 / 100.0,
        curvature_s2=-1.0 / 100,
        thickness=10,
        material=material.schott["N-BK7"],
    )
    lens2 = SphericalLens(
        radius=25,
        curvature_s1=1.0 / 100.0,
        curvature_s2=-1.0 / 100,
        thickness=10,
        material=material.schott["N-BK7"],
    )
    ccd = CCD()
    aperture_place = Stop(shape=Circular(radius=30), ap_shape=Circular(radius=25))
    s = System(
        complist=[
            (lens1, (0, 0, 100), (0, 0, 0)),
            (lens2, (0, 0, 120), (0, 0, 0)),
            (aperture_place, (0, 0, 110), (0, 0, 0)),
            (ccd, (0, 0, 150), (0, 0, 0)),
        ],
        n=1,
    )
    chief_ray = calc.chief_ray_search(s, aperture_place, (0, 10, 0), (0, -1, 1))

    np.testing.assert_almost_equal(chief_ray.origin, [0, 10, 0])
    np.testing.assert_almost_equal(
        chief_ray.direction, [3.58848263e-04, -9.26093228e-02, 9.95702458e-01], decimal=3
    )
    np.testing.assert_almost_equal(chief_ray.intensity, 1)
    np.testing.assert_almost_equal(chief_ray.wavelength, 0.58929)
    np.testing.assert_almost_equal(chief_ray.n, 1)
    np.testing.assert_equal(chief_ray.label, "")
    np.testing.assert_equal(chief_ray.orig_surf, None)
    np.testing.assert_almost_equal(chief_ray.order, 0)


@pytest.mark.skip(reason="Please write a proper test for pupil_location")
def test_pupil_location():
    assert False


def test_paraxial_location():
    lens1 = library.Edmund["45179"]  # f=200 r= 25
    optical_axis = Ray(origin=(0, 0, -10000), direction=(0, 0, 1), wavelength=0.55)
    ccd = CCD(size=(10, 10))
    s = System(
        complist=[
            (lens1, (0, 0, 100), (0, np.pi, 0)),
            (ccd, (0, 0, 320.053), (0, 0, 0)),
        ],
        n=1,
    )
    PB = parallel_beam_c(
        origin=(0, 0, 50),
        direction=(0, 0, 0),
        size=(15, 15),
        num_rays=(15, 15),
        wavelength=0.55,
    )
    s.ray_add(PB)
    s.propagate()

    image_location, real_ = calc.paraxial_location(s, optical_axis)
    np.testing.assert_almost_equal(image_location, [-5.59109334e-16, 0, 3.07249900e02])
    assert not real_


@pytest.mark.skip(reason="Test for find_aperture is pending.")
def test_find_aperture():
    for p in permutations([11, 13, 17, 19]):
        ccd_size = (p[0], p[1])
        aperture_size = (p[2], p[3])
        ccd = CCD(size=ccd_size)
        result = calc.find_aperture(ccd, size=aperture_size)
        expected = np.zeros_like(result)
        np.testing.assert_equal(result, expected)

        assert result.shape == aperture_size


def test_find_ppp():
    lens1 = library.Edmund["45179"]  # f=200 r= 25
    optical_axis = Ray(origin=(0, 0, -10), direction=(0, 0, 1), wavelength=0.55)
    s = System(
        complist=[
            (lens1, (0, 0, 0), (0, np.pi, 0)),
        ],
        n=1,
    )

    result = calc.find_ppp(s, optical_axis)
    print(result)
    np.testing.assert_almost_equal(result, [0, 0.0, 4.56699768])


@pytest.mark.skip(reason="Please write a proper test for get_optical_path_ep")
def test_get_optical_path_ep():
    assert False


@pytest.mark.skip(reason="Please write a proper test for find_reference_sphere_radius")
def test_find_reference_sphere_radius():
    assert False


@pytest.mark.skip(reason="Please write a proper test for aux_paral_f")
def test_aux_paral_f():
    assert False


@pytest.mark.skip(reason="Please write a proper test for parallel_propagate")
def test_parallel_propagate():
    assert False


@pytest.mark.skip(reason="Please write a proper test for aux_paral_f_ns")
def test_aux_paral_f_ns():
    assert False


@pytest.mark.skip(reason="Please write a proper test for parallel_propagate_ns")
def test_parallel_propagate_ns():
    assert False
