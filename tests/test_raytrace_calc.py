# standard imports
from itertools import permutations

# third-party imports
import numpy as np

# local imports
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
    # TODO: choose better example
    # Intersecting
    expected_intersection_point = (0, 0, 0)
    ray1 = Ray(pos=expected_intersection_point, dir=(0, 0.2, 1))
    ray2 = Ray(pos=expected_intersection_point, dir=(1, 2, 1))
    intersection_point, real_ = calc.intersection(ray1, ray2)
    np.testing.assert_almost_equal(intersection_point, expected_intersection_point)
    assert real_ == True

    # Not intersecting
    ray1 = Ray(pos=(0, 0, 0), dir=(0, 0.2, 1))
    ray2 = Ray(pos=(1, 0, 0), dir=(0, 0.2, 1))
    intersection_point, real_ = calc.intersection(ray1, ray2)
    np.testing.assert_almost_equal(intersection_point, [np.nan, np.nan, np.nan])
    assert real_ == False


def test_nearest_points():
    # Real closest point
    ray1 = Ray(pos=(0, 0, 0), dir=(1, 1, 0))
    ray2 = Ray(pos=(1, 0, 1), dir=(0, 1, 0))
    closest_point_on_ray1, closest_point_on_ray2, distance, real_ = calc.nearest_points(
        ray1, ray2
    )
    np.testing.assert_almost_equal(closest_point_on_ray1, [1.0, 1.0, 0.0])
    np.testing.assert_almost_equal(closest_point_on_ray2, [1.0, 1.0, 1.0])
    assert distance == 1
    assert real_ == True

    # Virtual closest point
    ray1 = Ray(pos=(0, 0, 0), dir=(1, 1, 0))
    ray2 = Ray(pos=(1, 10, 1), dir=(0, 1, 0))
    closest_point_on_ray1, closest_point_on_ray2, distance, real_ = calc.nearest_points(
        ray1, ray2
    )
    np.testing.assert_almost_equal(closest_point_on_ray1, [1.0, 1.0, 0.0])
    np.testing.assert_almost_equal(closest_point_on_ray2, [1.0, 1.0, 1.0])
    assert distance == 1
    assert real_ == False


def test_chief_ray_search():
    lens1 = SphericalLens(
        radius=25,
        curvature_s1=1.0 / 100.0,
        curvature_s2=-1.0 / 100,
        thickness=10,
        material=material.schott["BK7"],
    )
    lens2 = SphericalLens(
        radius=25,
        curvature_s1=1.0 / 100.0,
        curvature_s2=-1.0 / 100,
        thickness=10,
        material=material.schott["BK7"],
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


# def test_pupil_location():
#     # pupil_location(opsys,ccds,opaxis)
#     assert False, "Please write a proper test"


def test_paraxial_location():
    # image_location, real_ = paraxial_location(opsys, opaxis):
    lens1 = library.Edmund.get("45179")  # f=200 r= 25
    optical_axis = Ray(pos=(0, 0, -10000), dir=(0, 0, 1), wavelength=0.55)
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
    assert real_ == False


# TODO better test. Here is a tentative to better understand what does this function using a blackbox approach.
def test_find_apperture():
    # def find_apperture(ccd, size=(5,5)):

    # test that it always return a zeros array
    for p in permutations([11, 13, 17, 19]):
        ccd_size = (p[0], p[1])
        aperture_size = (p[2], p[3])
        ccd = CCD(size=ccd_size)
        result = calc.find_apperture(ccd, size=aperture_size)
        expected = np.zeros_like(result)
        np.testing.assert_equal(result, expected)

    # test the shape of the output array
    for p in permutations([11, 13, 17, 19]):
        ccd_size = (p[0], p[1])
        aperture_size = (p[2], p[3])
        ccd = CCD(size=ccd_size)
        result = calc.find_apperture(ccd, size=aperture_size)
        print(p, aperture_size, result.shape)
        if result.shape != aperture_size:
            print("---")
        assert result.shape == aperture_size
    # mostly work, but not all times. Looks like a bug.


def test_find_ppp():
    # find_ppp(opsys, opaxis)
    lens1 = library.Edmund.get("45179")  # f=200 r= 25
    optical_axis = Ray(pos=(0, 0, -10000), dir=(0, 0, 1), wavelength=0.55)
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

    result = calc.find_ppp(s, optical_axis)
    np.testing.assert_almost_equal(result, [0.0, 0.0, 104.5670357])


# def test_get_optical_path_ep():
#     get_optical_path_ep(opsys, opaxis, raylist, stop=None, r=None)
#     assert False, "Please write a proper test"


# def test_find_reference_sphere_radius():
#     find_reference_sphere_radius(ip, pl):
#     assert False, "Please write a proper test"


# def test_aux_paral_f():
#     aux_paral_f(x):
#     assert False, "Please write a proper test"


# def test_parallel_propagate():
#     parallel_propagate(os,r , np=None):
#     assert False, "Please write a proper test"


# def test_aux_paral_f_ns():
#     aux_paral_f_ns(x):
#     assert False, "Please write a proper test"


# def test_parallel_propagate_ns():
#     parallel_propagate_ns(os,rg, dp, r, np=None):
#     assert False, "Please write a proper test"
